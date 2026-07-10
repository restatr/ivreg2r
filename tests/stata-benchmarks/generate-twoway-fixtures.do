/*===========================================================================
  generate-twoway-fixtures.do
  ---------------------------
  Generates CSV benchmark fixtures for two-way cluster-robust VCE testing.
  Uses Cameron-Gelbach-Miller (2006) formula via Stata's ivreg2 cluster2().

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-twoway-fixtures.do

  ---------------------------------------------------------------------------
  M-14 RE-BASE (2026-07-05)
  ---------------------------------------------------------------------------
  The prior base for this family, `sim_twoway` (a simulated 25-firm x 20-obs
  panel with correlated firm/year effects), is RETIRED per Frank's ruling
  2026-07-05. The DGP was not
  reproducible: covariate draws were assigned after an unstable `sort
  year_id`, so re-running the generator remapped the seeded RNG draws to
  different firm/year cells and the checked-in fixtures could not be
  regenerated. The recorded known-dirty gate (Hansen J suppressed on the
  weighted cells under a rank-deficient S) was also shown to be
  draw-dependent rather than a stable property of the design: a re-run draw
  instead hit a ranktest degeneracy (missing KP Wald F / first-stage SE), and
  a stable-sort draw was dirtier still. A real-data alternative on this
  project's existing dataset stable was probed and rejected: abdata H88 +
  cluster(id year) is structurally rank-deficient on every cell (the d2.
  operators leave ~7 effective years against L=7 instruments).

  This family now re-bases onto the `cigar` panel (Baltagi & Levin 1992;
  Baltagi, Griffin & Xiong 2000; distributed with the R `plm` package as
  `Cigar`; also the empirical-example dataset in the GFIC focused-moment-
  selection paper, which motivates the price-endogeneity specification used
  below). The panel is 46 US states x 30 years (1963-1992), cached with
  verified provenance at ../validation/data/cigar.dta (source of record).
  No data CSV is exported here: the R
  side tests against the bundled `data(cigar)`, and the Stata side reads the
  same cached .dta -- both trace to one file of record, so there is nothing
  to duplicate into a fixture CSV (machine-readable-fixtures rule; the rule
  concerns computed *results*, not re-exporting input data that already has
  a canonical machine-readable copy).

  The canonical OLS two-way-cluster parity demo remains hf-owned (nlswork,
  help-file cells H104/H105) and is untouched by this file. The cells below
  are D5a: there is no upstream Stata/ivreg2 worked example of IV estimation
  with two-way clustering, so this file supplies the only such coverage,
  using the cigar/GFIC price-endogeneity specification as the running
  example. Gate: expected-clean -- confirmed by probes on 2026-07-05 that
  every cell below computes every requested diagnostic with zero warnings
  (no rank-deficiency, no suppressed statistics, unlike the retired DGP).

  Compositional retirements (documented, not generated): the wt x small and
  dof x small intersections are NOT generated as separate cells. The small-
  sample correction is applied uniformly in the shared R VCV assembly with
  no weight- or dofminus-specific branch; it is exercised for the two-way
  case by the cl2_small cell below, and the correction-factor property is
  covered independently of any fixture by a fixture-free test in
  test-vcov-cluster-twoway.R (M-15 review precedent). Crossing `small` with
  `wt` or `dof` again would test the same shared code path a third/fourth
  time without new coverage.
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"


/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results with two-way cluster metadata
---------------------------------------------------------------------------*/
capture program drop save_twoway_results
program define save_twoway_results
    syntax, prefix(string) suffix(string) outdir(string)

    local N = e(N)
    local K = e(rankxx)
    local L = e(rankzz)

    // --- Coefficients and SEs ---
    matrix b = e(b)
    matrix V = e(V)
    local names : colnames b
    local ncols = colsof(b)

    quietly {
        preserve
        clear
        set obs `ncols'
        gen str32 term = ""
        gen double estimate = .
        gen double std_error = .

        forvalues i = 1/`ncols' {
            local nm : word `i' of `names'
            replace term = "`nm'" in `i'
            replace estimate = b[1, `i'] in `i'
            replace std_error = sqrt(V[`i', `i']) in `i'
        }
        export delimited using "`outdir'/`prefix'_coef_`suffix'.csv", replace
        restore
    }

    // --- Full VCV matrix ---
    matrix V = e(V)
    quietly {
        preserve
        clear
        local dim = rowsof(V)
        svmat double V
        gen str32 term = ""
        forvalues i = 1/`dim' {
            local nm : word `i' of `names'
            replace term = "`nm'" in `i'
        }
        order term
        forvalues i = 1/`dim' {
            local nm : word `i' of `names'
            local cnm = subinstr("`nm'", ".", "_", .)
            rename V`i' vcov_`cnm'
        }
        export delimited using "`outdir'/`prefix'_vcov_`suffix'.csv", replace
        restore
    }

    // --- Diagnostics (with N_clust1 and N_clust2) ---
    quietly {
        preserve
        clear
        set obs 1

        gen double N = `N'
        gen double K = `K'
        gen double L = `L'
        gen double r2 = e(r2)
        gen double r2_a = e(r2_a)
        gen double r2u = e(r2u)
        gen double r2c = e(r2c)
        gen double rss = e(rss)
        gen double rmse = e(rmse)
        gen double sigmasq = e(rmse)^2
        gen double F_stat = e(F)
        gen double F_p = e(Fp)
        gen double F_df1 = e(Fdf1)
        gen double F_df2 = e(Fdf2)

        // Initialize all conditional columns
        gen double sargan = .
        gen double sarganp = .
        gen double sargandf = .
        gen double j = .
        gen double jp = .
        gen double jdf = .
        gen double idstat = .
        gen double idp = .
        gen double iddf = .
        gen double cdf = .
        gen double widstat = .
        gen double arf = .
        gen double arfp = .
        gen double archi2 = .
        gen double archi2p = .
        gen double ardf = .
        gen double ardf_r = .
        gen double cstat = .
        gen double cstatp = .
        gen double cstatdf = .
        gen double N_clust = .
        gen double N_clust1 = .
        gen double N_clust2 = .
        gen double sstat = .
        gen double sstatp = .
        gen double sstatdf = .

        // Overidentification
        capture replace sargan = e(sargan)
        capture replace sarganp = e(sarganp)
        capture replace sargandf = e(sargandf)
        capture replace j = e(j)
        capture replace jp = e(jp)
        capture replace jdf = e(jdf)

        // Underidentification
        capture replace idstat = e(idstat)
        capture replace idp = e(idp)
        capture replace iddf = e(iddf)

        // Weak identification
        capture replace cdf = e(cdf)
        capture replace widstat = e(widstat)

        // Anderson-Rubin
        capture replace arf = e(arf)
        capture replace arfp = e(arfp)
        capture replace archi2 = e(archi2)
        capture replace archi2p = e(archi2p)
        capture replace ardf = e(ardf)
        capture replace ardf_r = e(ardf_r)

        // Endogeneity / C-statistic
        capture replace cstat = e(cstat)
        capture replace cstatp = e(cstatp)
        capture replace cstatdf = e(cstatdf)

        // Endogeneity test (endog() option)
        gen double estat = .
        gen double estatp = .
        gen double estatdf = .
        capture replace estat = e(estat)
        capture replace estatp = e(estatp)
        capture replace estatdf = e(estatdf)

        // Stock-Wright S statistic
        capture replace sstat = e(sstat)
        capture replace sstatp = e(sstatp)
        capture replace sstatdf = e(sstatdf)

        // Cluster counts
        capture replace N_clust = e(N_clust)
        capture replace N_clust1 = e(N_clust1)
        capture replace N_clust2 = e(N_clust2)

        // Small option flag
        local is_small = (strpos("`suffix'", "small") > 0)
        gen double small = `is_small'

        export delimited using ///
            "`outdir'/`prefix'_diagnostics_`suffix'.csv", replace
        restore
    }

    // --- First-stage diagnostics (if IV model) ---
    capture confirm matrix e(first)
    if _rc == 0 {
        quietly {
            preserve
            clear

            matrix F = e(first)
            local frows = rowsof(F)
            local fcols = colsof(F)
            local fnames : colnames F
            local rnames : rownames F

            set obs `frows'
            gen str32 statistic = ""
            forvalues i = 1/`frows' {
                local rn : word `i' of `rnames'
                replace statistic = "`rn'" in `i'
            }

            forvalues j = 1/`fcols' {
                local cn : word `j' of `fnames'
                local cnclean = subinstr("`cn'", ".", "_", .)
                gen double `cnclean' = .
                forvalues i = 1/`frows' {
                    replace `cnclean' = F[`i', `j'] in `i'
                }
            }

            export delimited using ///
                "`outdir'/`prefix'_firststage_`suffix'.csv", replace
            restore
        }
    }
end


/*===========================================================================
  Data: cigar panel
  46 states x 30 years (1963-1992), Baltagi-Levin / Baltagi-Griffin-Xiong.
===========================================================================*/
display _newline(2) "=== Loading cigar panel data ==="

capture use "../validation/data/cigar.dta", clear
if _rc {
    display as error "cigar.dta cache not found at ../validation/data/."
    display as error "Regenerate it from plm::Cigar by running, from the repo root:"
    display as error `"  Rscript -e "setwd('pkg'); source('data-raw/cigar.R')""'
    exit 601
}
quietly tsset state year

// Derived columns in double precision, exactly as the R tests construct
// them from data(cigar) (F4 coherence discipline).
gen double lsales  = ln(sales)
gen double lrprice = ln(price/cpi)
gen double lrndi   = ln(ndi/cpi)
gen double lrpimin = ln(pimin/cpi)

// Deterministic aweight (the deterministic-weight convention used across the fixture generators)
gen double cwt = mod(state, 4) + 1


/*===========================================================================
  FIXTURE 1: IV overid, two-way cluster (small=FALSE)
  D5a: no upstream worked example of IV + two-way clustering exists; this
  is the IV x two-way parity cell, carrying the full diagnostics battery
  (the endogenous-price-demand specification follows the Baltagi-Levin /
  GFIC cigarette-demand lineage). This cell's firststage CSV is KEPT: it is
  consumed by the two-way first-stage-F parity test in
  test-vcov-cluster-twoway.R.
===========================================================================*/
display _newline(2) "=== Fixture: cigar_cl2 ==="
ivreg2 lsales lrndi (lrprice = lrpimin l.lrprice), cluster(state year) ///
    first endog(lrprice)
save_twoway_results, prefix(cigar) suffix(cl2) outdir(`outdir')


/*===========================================================================
  FIXTURE 2: IV overid, two-way cluster (small=TRUE)
  D5a; Baltagi-Levin lineage. Small-sample-corrected companion to cl2.
  firststage CSV erased below: first-stage-F fixture coverage lives in
  generate-firststage-fixtures.do,
  and the no-orphans rule says a CSV with no test consumer must
  not be checked in.
===========================================================================*/
display _newline(2) "=== Fixture: cigar_cl2_small ==="
ivreg2 lsales lrndi (lrprice = lrpimin l.lrprice), cluster(state year) ///
    first small endog(lrprice)
save_twoway_results, prefix(cigar) suffix(cl2_small) outdir(`outdir')
capture erase "`outdir'/cigar_firststage_cl2_small.csv"


/*===========================================================================
  FIXTURE 3: IV overid weighted, two-way cluster (small=FALSE)
  D5a; Baltagi-Levin lineage. This intersection (weighted x two-way IV) was
  known-dirty on the retired sim_twoway draw (Hansen J suppressed under a
  rank-deficient S); that gate was an artifact of the retired synthetic
  draw, not a property of weighted two-way clustering. Expected-clean here.
  firststage CSV erased: no test consumer.
===========================================================================*/
display _newline(2) "=== Fixture: cigar_cl2_wt ==="
ivreg2 lsales lrndi (lrprice = lrpimin l.lrprice) [aw=cwt], ///
    cluster(state year) first endog(lrprice)
save_twoway_results, prefix(cigar) suffix(cl2_wt) outdir(`outdir')
capture erase "`outdir'/cigar_firststage_cl2_wt.csv"


/*===========================================================================
  FIXTURE 4: IV overid dofminus(2) sdofminus(1), two-way cluster (small=FALSE)
  D5a; Baltagi-Levin lineage. firststage CSV erased: no test consumer.
===========================================================================*/
display _newline(2) "=== Fixture: cigar_cl2_dof ==="
ivreg2 lsales lrndi (lrprice = lrpimin l.lrprice), cluster(state year) ///
    first endog(lrprice) dofminus(2) sdofminus(1)
save_twoway_results, prefix(cigar) suffix(cl2_dof) outdir(`outdir')
capture erase "`outdir'/cigar_firststage_cl2_dof.csv"


/*===========================================================================
  Clean up and done
===========================================================================*/
display _newline(2) "=== All two-way cluster fixtures generated ==="
display "Output directory: `outdir'"
