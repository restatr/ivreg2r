/*===========================================================================
  generate-iweight-fixtures.do
  ----------------------------
  Generates CSV benchmark fixtures for importance weights (iweight; family
  M-12, re-based 2026-07-04 per planning/22-spec-matrix.md).

  Canonical base: mroz 2SLS baseline H31 (Stata ivreg2 help.txt line 1274).
  Stata documents the iweight option for ivreg2 but provides no worked
  example, so cells are D5a option-variation on the H31 base, using
  synthetic deterministic integer weights. The retired audit-20
  anti-pattern (passing the non-integer `wage` variable as an iweight) is
  not ported.

  Self-verifying property: an integral iweight is equivalent to a fweight
  of the same values. This is asserted at generation time (see the
  self-check block below) rather than merely assumed.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-iweight-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

/*===========================================================================
  Load mroz data
===========================================================================*/
capture use "../validation/data/mroz.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/wooldridge/mroz.dta, clear
    if _rc {
        display as error "Could not load mroz dataset (no local cache, no network)."
        exit 601
    }
}

// Deterministic synthetic weights derived from age (reproducible in R as
// (age %% 5) + 1). iwint is integer-valued in {1, ..., 5}. iwfrac adds 0.7
// so its sum over the 428-row estimation sample is genuinely non-integral
// (428 * 0.7 = 299.6) and the N = floor(sum(w)) convention binds -- a +1.5
// offset would sum integrally on an even sample size and leave floor() a
// no-op.
gen int iwint = mod(age, 5) + 1
gen double iwfrac = mod(age, 5) + 1.7


/*---------------------------------------------------------------------------
  Helper program: extract all results from ivreg2 and save to CSV
---------------------------------------------------------------------------*/
capture program drop save_ivreg2_results
program define save_ivreg2_results
    syntax, prefix(string) suffix(string) outdir(string)

    // Number of observations
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
        // Rename columns
        forvalues i = 1/`dim' {
            local nm : word `i' of `names'
            local cnm = subinstr("`nm'", ".", "_", .)
            rename V`i' vcov_`cnm'
        }
        export delimited using "`outdir'/`prefix'_vcov_`suffix'.csv", replace
        restore
    }

    // --- Diagnostics ---
    quietly {
        preserve
        clear
        set obs 1

        // Model stats
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

        // Initialize all conditional columns to missing
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

        // Endogeneity test
        gen double estat = .
        gen double estatp = .
        gen double estatdf = .
        capture replace estat = e(estat)
        capture replace estatp = e(estatp)
        capture replace estatdf = e(estatdf)

        // Stock-Wright S statistic
        gen double sstat = .
        gen double sstatp = .
        gen double sstatdf = .
        capture replace sstat = e(sstat)
        capture replace sstatp = e(sstatp)
        capture replace sstatdf = e(sstatdf)

        // Cluster info
        capture replace N_clust = e(N_clust)

        // Small option flag
        gen double small = ("`suffix'" == "iid_small")

        export delimited using ///
            "`outdir'/`prefix'_diagnostics_`suffix'.csv", replace
        restore
    }
end


/*===========================================================================
  FIXTURE SET 1: iweight OLS (integer weights)
  D5a option-variation, no upstream example; OLS companion to the H31 base.
  Model: lwage ~ exper expersq [iw=iwint]
===========================================================================*/
display _newline(2) "=== iweight OLS (integer weights) ==="

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq [iw=iwint]
save_ivreg2_results, prefix(mroz_iweight_ols) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq [iw=iwint], small
save_ivreg2_results, prefix(mroz_iweight_ols) suffix(iid_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 2: iweight 2SLS overidentified (integer weights)
  H31 model (help.txt:1274) with [iw=iwint], D5a.
  Model: lwage ~ exper expersq (educ = age kidslt6 kidsge6) [iw=iwint]
===========================================================================*/
display _newline(2) "=== iweight 2SLS overidentified (integer weights) ==="

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6) [iw=iwint]
save_ivreg2_results, prefix(mroz_iweight_overid) suffix(iid) outdir(`outdir')
// Capture e() from THIS fit for the self-check below, so the check is tied
// to the exact fixture fit rather than a re-typed copy that could drift.
matrix b_iw = e(b)
matrix V_iw = e(V)
local N_iw = e(N)

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6) [iw=iwint], small
save_ivreg2_results, prefix(mroz_iweight_overid) suffix(iid_small) outdir(`outdir')


/*===========================================================================
  Self-check: integral iweight is equivalent to fweight
  M-12's self-verifying property (planning/22-spec-matrix.md): an iweight
  whose values are all integers must give the same e(b), e(V), and e(N) as
  an fweight built from the same values. Compares the fweight fit against
  the b_iw/V_iw/N_iw captured from the FIXTURE SET 2 iid fit above.
===========================================================================*/
display _newline(2) "=== self-check: integral iweight == fweight ==="

ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6) [fw=iwint]
// assert/scalar expression parsers cannot take e(b)/e(V) (matrix-returning
// operators) inline; copy them to named matrices, then bind mreldif scalars.
matrix b_fw = e(b)
matrix V_fw = e(V)
scalar d_b = mreldif(b_iw, b_fw)
scalar d_V = mreldif(V_iw, V_fw)
assert d_b < 1e-12
assert d_V < 1e-12
assert `N_iw' == e(N)
display "self-check passed: integral iweight == fweight (b, V, N)"


/*===========================================================================
  FIXTURE SET 3: iweight fractional-weight N convention
  Pins Stata's N = floor(sum(w)) convention for non-integral iweight sums
  (the integer cells above cannot test this; previously pinned by the
  retired card_iweight_* [iw=wage] fixtures).
  Model: lwage ~ exper expersq [iw=iwfrac]
===========================================================================*/
display _newline(2) "=== iweight fractional-weight N convention ==="

ivreg2 lwage exper expersq [iw=iwfrac]
save_ivreg2_results, prefix(mroz_iweight_frac) suffix(iid) outdir(`outdir')


/*===========================================================================
  Done
===========================================================================*/
display ""
display "=========================================="
display "iweight fixtures generated successfully"
display "=========================================="
