/*===========================================================================
  generate-fixtures.do
  --------------------
  Generates CSV benchmark fixtures for ivreg2r R package testing.
  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)
  Each fixture produces:
    - *_data.csv          (dataset, for simulated data only)
    - *_coef_*.csv        (coefficients and SEs per VCE x small combo)
    - *_vcov_*.csv        (full VCV matrix per combo)
    - *_diagnostics_*.csv (test statistics per combo)
    - *_firststage_*.csv  (first-stage diagnostics per combo)

  Execution-order dependency: the sim_multi_endo_data.csv and
  sim_cluster_data.csv exports written here are IMPORTED by the rf,
  redundancy, and orthog generators (see planning/18-fixture-audit.md's
  cross-file dependency table). Run this file before those; do not delete the
  sim data sections without checking the downstream consumers.

  M-13 re-base (2026-07-05): the sim_cluster ESTIMATION cells (FIXTURE 5's
  run_all_vce_combos call) were retired and replaced by new ab_cl cells fit
  on abdata (see planning/22-spec-matrix.md row M-13). sim_cluster_data.csv
  itself is RETAINED, since generate-rf-fixtures.do imports it (M-21) and
  FIXTURE 9 (sim_cluster_dofminus) re-imports it (M-29); only the direct
  estimation cells moved.

  M-10 re-base (2026-07-06): the card_just_id and card_overid baseline cells
  (formerly FIXTURE 1 and FIXTURE 2) were retired -- the card 2SLS baseline
  family is absorbed into the helpfile suite (planning/22-spec-matrix.md row
  M-10). 2SLS baseline Stata parity now lives in the hf suite (hf_mroz H31/
  H41, hf_gril H03, generate-helpfile-fixtures.do) plus the new mroz_justid
  D5a cells added there. The small variants were retired by invariance: only
  rmse/sigmasq (an exact N/(N-K) factor) and F_stat floating-point noise
  differed under small in the retired fixtures, with every other diagnostic
  byte-identical (verified 2026-07-06).

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

// Pre-load Card data and save before defining programs.
// bcuse internally calls "clear all" which drops user-defined programs,
// so we must call it before any program definitions.
capture bcuse card, clear
if _rc != 0 {
    display as error "Could not load Card dataset via bcuse."
    display as error "Install bcuse (ssc install bcuse) and rerun."
    exit 601
}
save "`outdir'/_card_temp.dta", replace

/*---------------------------------------------------------------------------
  Helper program: full-precision data export. The default `export delimited`
  writes numeric columns at 8-digit precision, which truncates float32 bit
  patterns (ticket F4 / the CSV float-truncation gotcha in CLAUDE.md).
  Defined here, immediately after the bcuse card load above and before its
  export below (bcuse's clear-all already happened, so program definitions
  are safe from this point on).
---------------------------------------------------------------------------*/
capture program drop export_full_precision
program define export_full_precision
    // Full-precision data export: the default `export delimited` writes
    // numeric columns at 8-digit precision, which truncates float32 bit
    // patterns (ticket F4 / the CSV float-truncation gotcha in CLAUDE.md).
    syntax using/
    quietly ds, has(type numeric)
    format `r(varlist)' %21.0g
    export delimited using "`using'", replace datafmt
end

export_full_precision using "`outdir'/card_data.csv"

/*---------------------------------------------------------------------------
  Pre-load abdata (H88 panel) and save before defining the remaining
  programs, for the new ab_cl cells (M-13 re-base). Imitates the card
  pre-load above; falls back to the bc.edu mirror when the bundled
  validation copy is absent (see the abdata block in
  generate-liml-fixtures.do for precedent).
---------------------------------------------------------------------------*/
capture use "../validation/data/abdata.dta", clear
if _rc {
    use http://fmwww.bc.edu/ec-p/data/macro/abdata.dta, clear
}
quietly tsset id year
save "`outdir'/_ab_temp.dta", replace

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
            // Clean name for variable naming (replace . with _)
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

        // Model stats (always returned by ivreg2)
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

        // Initialize all conditional columns to missing (.) so every CSV
        // has the same schema regardless of model specification.
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

        // Overidentification (conditionally returned)
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

        // Endogeneity / C-statistic (orthog() option → e(cstat))
        capture replace cstat = e(cstat)
        capture replace cstatp = e(cstatp)
        capture replace cstatdf = e(cstatdf)

        // Endogeneity test (endog() option → e(estat))
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

        // Cluster info (only when clustered)
        capture replace N_clust = e(N_clust)

        // Small option flag
        gen double small = ("`suffix'" == "iid_small" | ///
                            "`suffix'" == "hc1_small" | ///
                            "`suffix'" == "cl_small")

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

            // e(first) matrix: rows are statistics, columns are endogenous vars
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


/*---------------------------------------------------------------------------
  Helper: run one ivreg2 specification, save results
  Uses globals $ivreg2_depvar $ivreg2_exog $ivreg2_endo $ivreg2_iv
  set by the caller before invoking run_all_vce_combos.
---------------------------------------------------------------------------*/
capture program drop run_all_vce_combos
program define run_all_vce_combos
    syntax, prefix(string) outdir(string) [clustervar(string) ///
            endogopt(string) modelopts(string)]

    // Build endog() option string
    local endogstr ""
    if "`endogopt'" != "" {
        local endogstr "endog(`endogopt')"
    }

    // Build the ivreg2 command from globals.
    // Globals are used because Stata's syntax parser cannot handle
    // nested parentheses inside string options.
    local cmd "ivreg2 $ivreg2_depvar $ivreg2_exog"
    if "$ivreg2_endo" != "" {
        local cmd "`cmd' ($ivreg2_endo = $ivreg2_iv)"
    }

    // --- IID, small=FALSE ---
    `cmd', first `modelopts' `endogstr'
    save_ivreg2_results, prefix(`prefix') suffix(iid) outdir(`outdir')

    // --- IID, small=TRUE ---
    `cmd', first `modelopts' small `endogstr'
    save_ivreg2_results, prefix(`prefix') suffix(iid_small) outdir(`outdir')

    // --- HC1, small=FALSE ---
    `cmd', first `modelopts' robust `endogstr'
    save_ivreg2_results, prefix(`prefix') suffix(hc1) outdir(`outdir')

    // --- HC1, small=TRUE ---
    `cmd', first `modelopts' robust small `endogstr'
    save_ivreg2_results, prefix(`prefix') suffix(hc1_small) outdir(`outdir')

    // --- Cluster (if cluster variable specified) ---
    if "`clustervar'" != "" {
        // Cluster, small=FALSE
        `cmd', first `modelopts' cluster(`clustervar') `endogstr'
        save_ivreg2_results, prefix(`prefix') suffix(cl) outdir(`outdir')

        // Cluster, small=TRUE
        `cmd', first `modelopts' cluster(`clustervar') small `endogstr'
        save_ivreg2_results, prefix(`prefix') suffix(cl_small) outdir(`outdir')
    }
end


/*===========================================================================
  FIXTURE 3: sim_multi_endo
  Simulated data with 2 endogenous variables and 4 excluded IVs
  Purpose: Tests SW/AP, Shea partial R2, multi-endo eigenvalue paths
===========================================================================*/
display _newline(2) "=== FIXTURE 3: sim_multi_endo ==="

clear
set seed 20260214
set obs 1000

// Exogenous variables
gen x1 = rnormal()
gen x2 = rnormal()

// Excluded instruments (4 IVs for 2 endogenous)
gen z1 = rnormal()
gen z2 = rnormal()
gen z3 = rnormal()
gen z4 = rnormal()

// Structural errors (correlated with endogenous vars)
gen u = rnormal()
gen v1 = rnormal()
gen v2 = rnormal()

// Endogenous variables (correlated with u via v1, v2)
gen endo1 = 0.5*z1 + 0.3*z2 + 0.2*x1 + 0.8*v1 + 0.3*u
gen endo2 = 0.4*z2 + 0.3*z3 + 0.2*z4 + 0.1*x2 + 0.7*v2 + 0.2*u

// Outcome
gen y = 1.0 + 0.5*x1 - 0.3*x2 + 1.2*endo1 - 0.8*endo2 + u

// Save simulated data. This data.csv is re-imported downstream
// (generate-rf-fixtures.do), so full precision matters.
export_full_precision using "`outdir'/sim_multi_endo_data.csv"

global ivreg2_depvar "y"
global ivreg2_exog "x1 x2"
global ivreg2_endo "endo1 endo2"
global ivreg2_iv "z1 z2 z3 z4"

run_all_vce_combos, ///
    prefix(sim_multi_endo) ///
    outdir(`outdir') ///
    endogopt(endo1 endo2)

// e(first) exports are orphans since the M-25 re-base (M-25 owns first-stage fixtures) -- erase them so a full re-run cannot recreate untracked orphans.
foreach s in iid iid_small hc1 hc1_small {
    erase "`outdir'/sim_multi_endo_firststage_`s'.csv"
}


/*===========================================================================
  FIXTURE 4: sim_no_constant
  Simulated data, model without intercept
  Purpose: Tests uncentered R-squared, no-constant code path
===========================================================================*/
display _newline(2) "=== FIXTURE 4: sim_no_constant ==="

clear
set seed 20260215
set obs 500

gen z1 = rnormal()
gen z2 = rnormal()
gen x1 = rnormal()
gen v = rnormal()
gen u = rnormal()

gen endo1 = 0.6*z1 + 0.4*z2 + 0.3*x1 + 0.7*v + 0.2*u
gen y = 2.0*x1 + 1.5*endo1 + u

// Save simulated data.
export_full_precision using "`outdir'/sim_no_constant_data.csv"

global ivreg2_depvar "y"
global ivreg2_exog "x1"
global ivreg2_endo "endo1"
global ivreg2_iv "z1 z2"

// noconstant model
run_all_vce_combos, ///
    prefix(sim_no_constant) ///
    outdir(`outdir') ///
    endogopt(endo1) ///
    modelopts(noconstant)

// e(first) exports are orphans since the M-25 re-base (M-25 owns first-stage fixtures) -- erase them so a full re-run cannot recreate untracked orphans.
foreach s in iid iid_small hc1 hc1_small {
    erase "`outdir'/sim_no_constant_firststage_`s'.csv"
}


/*===========================================================================
  FIXTURE 5: sim_cluster
  Simulated data with group structure (50 clusters)
  Purpose: Tests cluster VCV, cluster-robust diagnostics
===========================================================================*/
display _newline(2) "=== FIXTURE 5: sim_cluster ==="

clear
set seed 20260216
set obs 1000

// 50 clusters, 20 obs each
gen cluster_id = ceil(_n / 20)

// Cluster-level effects
gen alpha = rnormal() if mod(_n - 1, 20) == 0
bysort cluster_id: replace alpha = alpha[1]

// Individual-level variables
gen x1 = rnormal() + 0.3*alpha
gen z1 = rnormal() + 0.2*alpha
gen z2 = rnormal() + 0.1*alpha
gen v = rnormal()
gen u = rnormal() + 0.5*alpha

// Endogenous variable
gen endo1 = 0.5*z1 + 0.3*z2 + 0.2*x1 + 0.6*v + 0.3*u

// Outcome
gen y = 1.0 + 0.8*x1 + 1.5*endo1 + u

// Save simulated data. This data.csv is re-imported downstream
// (generate-rf-fixtures.do and FIXTURE 9 below), so full precision matters.
export_full_precision using "`outdir'/sim_cluster_data.csv"

// M-13 re-base (2026-07-05): the estimation cells formerly generated here
// via run_all_vce_combos (iid/iid_small/hc1/hc1_small/cl/cl_small) were
// retired -- IV 2SLS + one-way cluster Stata parity now lives on the
// abdata-based ab_cl cells below (FIXTURE 5b), per planning/22-spec-matrix.md
// row M-13. The DGP and sim_cluster_data.csv export above are KEPT: they
// feed generate-rf-fixtures.do (M-21) and FIXTURE 9 (M-29).


/*===========================================================================
  FIXTURE 5b: ab_cl (M-13 re-base, 2026-07-05)
  Dataset: abdata (H88 panel, help.txt:1541)
  Model: n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys) -- the H88 model minus
  gmm2s (M-15/M-16/M-17 precedent), i.e. plain 2SLS.
  cluster(id) anchor: H91 (help.txt:1558).
  Purpose: D5a option-variation (plain 2SLS x one-way cluster x K1=3
  multi-endogenous x endog-subset). This is the IV 2SLS + one-way cluster
  Stata-parity cell: it does not exist anywhere else post-sim-deletion
  (ab_liml/ab_cue/tsop_ab_gmm2s are other estimators). OLS + one-way
  cluster Stata parity is owned by the hf suite (H91/H103); the
  griliches cluster(year) M=7 known-dirty lesson (J + Stock-Wright
  suppressed) is owned by hf H27/H28 -- neither is duplicated here.
  `first endog(w)` is required: `first` populates e(sstat)/e(arf)/
  e(archi2) and endog(w) populates e(estat), feeding the five
  cross-family diagnostics consumers (M-11 card_fw precedent). The
  e(first) matrix itself is NOT shipped -- first-stage fixtures are
  owned by M-25 -- so the saver-written firststage CSVs are erased
  immediately below (hf H22 erase precedent).
===========================================================================*/
display _newline(2) "=== FIXTURE 5b: ab_cl ==="

use "`outdir'/_ab_temp.dta", clear
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), cluster(id) first endog(w)
save_ivreg2_results, prefix(ab_cl) suffix(cl) outdir(`outdir')
erase "`outdir'/ab_cl_firststage_cl.csv"

use "`outdir'/_ab_temp.dta", clear
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), cluster(id) small first endog(w)
save_ivreg2_results, prefix(ab_cl) suffix(cl_small) outdir(`outdir')
erase "`outdir'/ab_cl_firststage_cl_small.csv"


// FIXTURE 6 (card_just_id_weighted) retired 2026-07-05 by the M-11 re-base;
// weighted cells now live in generate-weight-type-fixtures.do.


/*===========================================================================
  FIXTURE 7: card_just_id_dofminus
  Dataset: Card (1995) — returns to education
  Model: lwage ~ exper expersq black south | educ | nearc4
  Purpose: dofminus(1) sdofminus(1) adjustments (just-identified)
===========================================================================*/
display _newline(2) "=== FIXTURE 7: card_just_id_dofminus ==="

// Reload Card data
use "`outdir'/_card_temp.dta", clear

global ivreg2_depvar "lwage"
global ivreg2_exog "exper expersq black south"
global ivreg2_endo "educ"
global ivreg2_iv "nearc4"

run_all_vce_combos, ///
    prefix(card_just_id_dofminus) ///
    outdir(`outdir') ///
    clustervar(smsa66) ///
    endogopt(educ) ///
    modelopts(dofminus(1) sdofminus(1))


/*===========================================================================
  FIXTURE 8: card_overid_dofminus
  Dataset: Card (1995)
  Model: lwage ~ exper expersq black south | educ | nearc2 nearc4
  Purpose: dofminus(1) sdofminus(1) adjustments (overidentified)
===========================================================================*/
display _newline(2) "=== FIXTURE 8: card_overid_dofminus ==="

// Reload Card data
use "`outdir'/_card_temp.dta", clear
global ivreg2_depvar "lwage"
global ivreg2_exog "exper expersq black south"
global ivreg2_endo "educ"
global ivreg2_iv "nearc2 nearc4"

run_all_vce_combos, ///
    prefix(card_overid_dofminus) ///
    outdir(`outdir') ///
    clustervar(smsa66) ///
    endogopt(educ) ///
    modelopts(dofminus(1) sdofminus(1))


/*===========================================================================
  FIXTURE 9: sim_cluster_dofminus
  Simulated data with group structure (50 clusters)
  Purpose: dofminus(1) sdofminus(1) with cluster-robust
===========================================================================*/
display _newline(2) "=== FIXTURE 9: sim_cluster_dofminus ==="

// Reload simulated cluster data
import delimited using "`outdir'/sim_cluster_data.csv", clear

global ivreg2_depvar "y"
global ivreg2_exog "x1"
global ivreg2_endo "endo1"
global ivreg2_iv "z1 z2"

// Run all combos including cluster
run_all_vce_combos, ///
    prefix(sim_cluster_dofminus) ///
    outdir(`outdir') ///
    clustervar(cluster_id) ///
    endogopt(endo1) ///
    modelopts(dofminus(1) sdofminus(1))


/*===========================================================================
  Clean up globals and done
===========================================================================*/
macro drop ivreg2_depvar ivreg2_exog ivreg2_endo ivreg2_iv
capture erase "`outdir'/_card_temp.dta"
capture erase "`outdir'/_ab_temp.dta"

display _newline(2) "=== All fixtures generated ==="
display "Output directory: `outdir'"
display "Run the R test suite to compare against these fixtures."
