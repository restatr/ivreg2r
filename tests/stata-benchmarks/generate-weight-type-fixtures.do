/*===========================================================================
  generate-weight-type-fixtures.do
  ---------------------------------
  Generates CSV benchmark fixtures for frequency weights (fweight),
  probability weights (pweight), and analytic weights (aweight).

  Canonical base: the Card (1995) returns-to-education overidentified spec,
  the same model used by M-10 (card_overid):
    lwage ~ exper expersq black south | educ | nearc4 nearc2
  `first endog(educ)` is required on every cell: `first` populates
  e(arf)/e(archi2)/e(sstat) and `endog(educ)` populates e(estat), both
  consumed by cross-family diagnostics tests. The e(first) matrix is NOT
  exported to CSV here -- M-25 owns first-stage-under-weights coverage.

  Stata documents the fweight/pweight/aweight options for ivreg2 but
  provides no worked weighted example, so cells are D6 option-variation on
  the M-10 base: pweight uses `weight`, Card's genuine NLS sampling weight;
  aweight is demonstrated via a self-verifying cell-means construction
  (collapse the micro data to cells, weight by cell count); fweight is
  demonstrated via a self-verifying duplicated-rows construction (weight by
  an integer replication count). Both self-verifying properties are
  asserted at generation time below, not merely assumed.

  Invariance retirement (verified 2026-07-05 by diffing the retired
  fixtures): Stata's [aw=w], robust and [pw=w] produce BYTE-IDENTICAL coef,
  vcov, and diagnostics CSVs, including the small and cluster(smsa)
  variants (card_aweight_overid_*_hc1{,_small} ==
  card_pweight_overid_*_hc0{,_small}; ..._cl{,_small} likewise). Therefore
  no aweight+robust or aweight+cluster cells are generated on the micro
  card data here: R pins those by fitting the aweight variant directly
  against the card_pw_* fixtures AND asserting the aweight fit equal to the
  pweight fit (see test-weight-types.R). The mandatory direct-equality
  assertion is the M-22 lesson: comparing both variants only against the
  fixture would loosen the invariant to Stata tolerance.

  Retired: the fweight/pweight/aweight cells built on `[?w=weight]` with
  binary smsa/smsa66 clusters (card_fweight_just_id/overid,
  card_pweight_just_id/overid, card_aweight_overid), and the
  card_just_id_weighted family in generate-fixtures.do. First-stage
  exports are dropped (M-25 owns first-stage-under-weights coverage), and
  just-identified cells are dropped (M-12 precedent: no weight-specific
  identification path; M-10 owns identification).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-weight-type-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

/*===========================================================================
  Load Card data BEFORE any program define (bcuse calls "clear all", which
  drops user-defined programs). Try the cached local copy first, then fall
  back to bcuse.
===========================================================================*/
capture use "../validation/data/card.dta", clear
if _rc != 0 {
    capture bcuse card, clear
    if _rc != 0 {
        display as error "Could not load Card dataset (no local cache, no network)."
        exit 601
    }
}

// Deterministic derived variables (reproduced identically in R):
//   region partitions the 3010 rows into 9 groups via the Census-region
//   dummies (M=9, non-degenerate, replacing the retired binary
//   smsa/smsa66 M=2 anti-pattern).
//   fwt is a deterministic integer replication count in {1, ..., 5} used
//   to build the fweight self-verifying demo.
gen byte region = reg661*1 + reg662*2 + reg663*3 + reg664*4 + reg665*5 + ///
                  reg666*6 + reg667*7 + reg668*8 + reg669*9
gen int fwt = mod(age, 5) + 1

tempfile cardbase
save `cardbase'


/*---------------------------------------------------------------------------
  Helper program: extract all results from ivreg2 and save to CSV
  (Adapted from generate-iweight-fixtures.do, the M-12 re-based template;
  no first-stage export block -- M-25 owns first-stage-under-weights
  coverage. The small-flag suffix list is extended to cover this family's
  four VCE families.)
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
        gen double small = ("`suffix'" == "iid_small" | ///
                            "`suffix'" == "hc1_small" | ///
                            "`suffix'" == "hc0_small" | ///
                            "`suffix'" == "cl_small")

        export delimited using ///
            "`outdir'/`prefix'_diagnostics_`suffix'.csv", replace
        restore
    }
end


/*===========================================================================
  FIXTURE SET 1: card_fw_* -- frequency weights (fweight), duplicated-rows
  construction (D6). Model: lwage ~ exper expersq black south | educ |
  nearc4 nearc2, [fw=fwt] where fwt = mod(age,5)+1. Six cells; the cluster
  cells use cluster(region) (M=9, non-degenerate, replacing the retired
  binary smsa/smsa66 M=2 anti-pattern).
===========================================================================*/
display _newline(2) "=== card_fw: frequency weights ==="

use `cardbase', clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=fwt], first endog(educ)
save_ivreg2_results, prefix(card_fw) suffix(iid) outdir(`outdir')
// Capture e() from THIS fit for the fweight-expand self-check below, so the
// check is tied to the exact fixture fit rather than a re-typed copy.
matrix b_fw_iid = e(b)
matrix V_fw_iid = e(V)
local N_fw_iid = e(N)

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=fwt], first endog(educ) small
save_ivreg2_results, prefix(card_fw) suffix(iid_small) outdir(`outdir')

// --- HC1 (robust), small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=fwt], first endog(educ) robust
save_ivreg2_results, prefix(card_fw) suffix(hc1) outdir(`outdir')
// Capture e() from THIS fit for the fweight-expand self-check below.
matrix b_fw_hc1 = e(b)
matrix V_fw_hc1 = e(V)
local N_fw_hc1 = e(N)

// --- HC1, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=fwt], first endog(educ) robust small
save_ivreg2_results, prefix(card_fw) suffix(hc1_small) outdir(`outdir')

// --- Cluster(region), small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=fwt], first endog(educ) cluster(region)
save_ivreg2_results, prefix(card_fw) suffix(cl) outdir(`outdir')

// --- Cluster(region), small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=fwt], first endog(educ) cluster(region) small
save_ivreg2_results, prefix(card_fw) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 2: card_pw_* -- probability weights (pweight), genuine NLS
  sampling weight (D6). Model: lwage ~ exper expersq black south | educ |
  nearc4 nearc2, [pw=weight] where `weight` is Card's genuine NLS sampling
  weight. pweight forces robust VCE, so there are no iid cells.
===========================================================================*/
display _newline(2) "=== card_pw: probability weights ==="

use `cardbase', clear

// --- HC0 (Stata's default for pweight), small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [pw=weight], first endog(educ)
save_ivreg2_results, prefix(card_pw) suffix(hc0) outdir(`outdir')

// --- HC0, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [pw=weight], first endog(educ) small
save_ivreg2_results, prefix(card_pw) suffix(hc0_small) outdir(`outdir')

// --- Cluster(region), small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [pw=weight], first endog(educ) cluster(region)
save_ivreg2_results, prefix(card_pw) suffix(cl) outdir(`outdir')

// --- Cluster(region), small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [pw=weight], first endog(educ) cluster(region) small
save_ivreg2_results, prefix(card_pw) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 3: card_aw_cells_* -- analytic weights (aweight), cell-means
  construction (D6). Micro data collapsed to
  (educ, exper, expersq, black, south, nearc4, nearc2) cells (1030 cells;
  all 3010 rows are complete cases on these variables), lwage averaged
  within cell, n = cell count. Model fit with [aw=n] on the collapsed data.
  No cluster cells here (see the invariance-retirement note in the file
  header).
===========================================================================*/
display _newline(2) "=== card_aw_cells: analytic weights (cell-means) ==="

use `cardbase', clear
// collapse stores (mean) results as float by default, which truncates the
// cell means to float32 (~1e-7 relative -- the F4/CSV-truncation lesson
// class; measured mreldif vs the micro fit: 1.1e-7 float vs 2.1e-12
// double). Take a double copy first so the means are full precision,
// matching the double-precision means R computes from the bundled data.
gen double lwage_d = lwage
collapse (mean) lwage_d (count) n=lwage_d, by(educ exper expersq black south nearc4 nearc2)
rename lwage_d lwage

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [aw=n], first endog(educ)
save_ivreg2_results, prefix(card_aw_cells) suffix(iid) outdir(`outdir')
// Capture e(b) from THIS fit for the aw-cells-recover-micro self-check below.
matrix b_cells = e(b)

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [aw=n], first endog(educ) small
save_ivreg2_results, prefix(card_aw_cells) suffix(iid_small) outdir(`outdir')

// --- HC1, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [aw=n], first endog(educ) robust
save_ivreg2_results, prefix(card_aw_cells) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [aw=n], first endog(educ) robust small
save_ivreg2_results, prefix(card_aw_cells) suffix(hc1_small) outdir(`outdir')


/*===========================================================================
  Self-check: fweight ([fw=fwt]) equivalent to duplicated rows (D6)
  Verifies e(b)/e(V)/e(N) from the ACTUAL card_fw_iid and card_fw_hc1
  fixture fits above equal the corresponding unweighted fit on the
  row-expanded data, for both the iid and robust VCE.
===========================================================================*/
display _newline(2) "=== self-check: fweight == expand(fwt) rows ==="

use `cardbase', clear
expand fwt

quietly ivreg2 lwage exper expersq black south (educ = nearc4 nearc2), first endog(educ)
matrix b_exp_iid = e(b)
matrix V_exp_iid = e(V)
scalar d_b_iid = mreldif(b_fw_iid, b_exp_iid)
scalar d_V_iid = mreldif(V_fw_iid, V_exp_iid)
assert d_b_iid < 1e-12
assert d_V_iid < 1e-12
assert `N_fw_iid' == e(N)

quietly ivreg2 lwage exper expersq black south (educ = nearc4 nearc2), first endog(educ) robust
matrix b_exp_hc1 = e(b)
matrix V_exp_hc1 = e(V)
scalar d_b_hc1 = mreldif(b_fw_hc1, b_exp_hc1)
scalar d_V_hc1 = mreldif(V_fw_hc1, V_exp_hc1)
assert d_b_hc1 < 1e-12
assert d_V_hc1 < 1e-12
assert `N_fw_hc1' == e(N)

display "self-check passed: fweight == expand(fwt) rows (b, V, N; iid & robust)"


/*===========================================================================
  Self-check: aweight cell-means recover the micro-data fit
  Coefficients of the collapsed [aw=n] fit must equal the unweighted
  micro-data fit (VCV is NOT compared -- it differs by construction; only
  b is preserved by the cell-means collapse). Compares b_cells captured
  from the card_aw_cells_iid fixture fit above against the unweighted
  micro fit computed here (not itself a fixture cell).
===========================================================================*/
display _newline(2) "=== self-check: aweight cell-means == micro-data fit (b only) ==="

use `cardbase', clear
quietly ivreg2 lwage exper expersq black south (educ = nearc4 nearc2), first endog(educ)
matrix b_micro = e(b)
scalar d_b_cells = mreldif(b_cells, b_micro)
assert d_b_cells < 1e-10

display "self-check passed: aweight cell-means == micro-data fit (b)"


/*===========================================================================
  Done
===========================================================================*/
display ""
display "=========================================="
display "Weight-type fixtures generated successfully"
display "=========================================="
