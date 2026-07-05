/*===========================================================================
  generate-cue-fixtures.do
  ------------------------
  Generates CSV benchmark fixtures for CUE estimation testing in ivreg2r.
  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-cue-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

// Pre-load Card data and save before defining programs.
// bcuse internally calls "clear all" which drops user-defined programs.
capture bcuse card, clear
if _rc != 0 {
    display as error "Could not load Card dataset via bcuse."
    display as error "Install bcuse (ssc install bcuse) and rerun."
    exit 601
}
save "`outdir'/_card_cue_temp.dta", replace

/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV
  (Same structure as other fixture generators)
---------------------------------------------------------------------------*/
capture program drop save_ivreg2_results
program define save_ivreg2_results
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
        local vr = rowsof(V)
        local vc = colsof(V)
        set obs `vr'
        forvalues j = 1/`vc' {
            gen double v`j' = .
            forvalues i = 1/`vr' {
                replace v`j' = V[`i', `j'] in `i'
            }
        }
        export delimited using "`outdir'/`prefix'_vcov_`suffix'.csv", replace
        restore
    }

    // --- Diagnostics ---
    quietly {
        preserve
        clear
        set obs 1

        // Overidentification
        gen double overid_stat = e(sargan)
        if overid_stat == . {
            replace overid_stat = e(j)
        }
        gen double overid_p = e(sarganp)
        if overid_p == . {
            replace overid_p = e(jp)
        }
        gen int overid_df = e(sargandf)
        if overid_df == . {
            replace overid_df = e(jdf)
        }

        // Underidentification (Anderson / KP)
        gen double underid_stat = e(idstat)
        gen double underid_p = e(idp)
        gen int underid_df = e(iddf)

        // Weak identification
        gen double weak_id_cd_f = e(cdf)
        gen double weak_id_kp_f = e(widstat)

        // First-stage F
        gen double first_stage_f = e(rkf)

        // Anderson-Rubin
        gen double ar_f = e(arf)
        gen double ar_f_p = e(arfp)
        gen int ar_f_df1 = e(ardf)
        gen int ar_f_df2 = e(ardf_r)
        gen double ar_chi2 = e(archi2)
        gen double ar_chi2_p = e(archi2p)

        // Stock-Wright S
        gen double sw_stat = e(sstat)
        gen double sw_p = e(sstatp)
        gen int sw_df = e(sstatdf)

        // Endogeneity
        gen double endog_stat = e(estat)
        gen double endog_p = e(estatp)
        gen int endog_df = e(estatdf)

        // Model F
        gen double model_f = e(F)
        gen double model_f_p = e(Fp)
        gen int model_f_df1 = e(Fdf1)
        gen int model_f_df2 = e(Fdf2)

        // Summary stats
        gen double sigma = e(rmse)
        gen double rss = e(rss)
        gen double r2 = e(r2)
        gen double r2_a = e(r2_a)
        gen double r2u = e(r2u)
        gen double r2c = e(r2c)

        // Counts
        gen int N = e(N)
        gen int K = `K'
        gen int L = `L'

        export delimited using "`outdir'/`prefix'_diagnostics_`suffix'.csv", replace
        restore
    }
end


/*===========================================================================
  CUE fixtures using Card data (overidentified model)
  Model: lwage ~ exper expersq black south (educ = nearc2 nearc4)
===========================================================================*/

// --- CUE IID (no robust) ---
use "`outdir'/_card_cue_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cue
save_ivreg2_results, prefix(card_overid_cue) suffix(iid) outdir(`outdir')

// --- CUE IID small ---
use "`outdir'/_card_cue_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cue small
save_ivreg2_results, prefix(card_overid_cue) suffix(iid_small) outdir(`outdir')

// --- CUE robust (HC0) ---
use "`outdir'/_card_cue_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cue robust
save_ivreg2_results, prefix(card_overid_cue) suffix(robust) outdir(`outdir')

// --- CUE robust small (HC1) ---
use "`outdir'/_card_cue_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cue robust small
save_ivreg2_results, prefix(card_overid_cue) suffix(robust_small) outdir(`outdir')

// --- CUE cluster(age) ---
use "`outdir'/_card_cue_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cue cluster(age)
save_ivreg2_results, prefix(card_overid_cue) suffix(cluster) outdir(`outdir')

// --- CUE cluster(age) small ---
use "`outdir'/_card_cue_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cue cluster(age) small
save_ivreg2_results, prefix(card_overid_cue) suffix(cluster_small) outdir(`outdir')


/*===========================================================================
  Just-identified CUE (should equal 2SLS)
  Model: lwage ~ exper expersq black south (educ = nearc4)
===========================================================================*/

// --- Just-identified CUE robust ---
use "`outdir'/_card_cue_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc4), cue robust
save_ivreg2_results, prefix(card_justid_cue) suffix(robust) outdir(`outdir')


/*===========================================================================
  Weighted CUE
===========================================================================*/

// --- CUE aweight robust ---
use "`outdir'/_card_cue_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4) [aw=age], cue robust
save_ivreg2_results, prefix(card_overid_cue_weighted) suffix(aw_robust) outdir(`outdir')


/*===========================================================================
  CUE with dofminus
===========================================================================*/

// --- CUE robust dofminus(2) ---
use "`outdir'/_card_cue_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cue robust dofminus(2)
save_ivreg2_results, prefix(card_overid_cue_dofminus) suffix(robust) outdir(`outdir')


/*===========================================================================
  CUE with endog test
===========================================================================*/

// --- CUE robust with endog test ---
use "`outdir'/_card_cue_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cue robust endog(educ)
save_ivreg2_results, prefix(card_overid_cue_endog) suffix(robust) outdir(`outdir')


/*===========================================================================
  CUE with orthog test
===========================================================================*/

// --- CUE robust with orthog test ---
use "`outdir'/_card_cue_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cue robust orthog(nearc2)
save_ivreg2_results, prefix(card_overid_cue_orthog) suffix(robust) outdir(`outdir')


/*===========================================================================
  HAC CUE (using time-series data)
===========================================================================*/

// Generate simulated time-series IV data (same DGP as GMM fixtures)
clear
set seed 12345
set obs 200
gen t = _n
tsset t

// Generate instruments
gen z1 = rnormal()
gen z2 = rnormal()

// AR(1) errors
gen u = rnormal()
replace u = 0.5 * u[_n-1] + rnormal() if t > 1

// Endogenous regressor
gen x = 0.3*z1 + 0.2*z2 + 0.5*u + rnormal()

// Outcome
gen y = 1.5 + 2.0*x + u

// Exogenous regressor
gen w = rnormal()

save "`outdir'/_ts_cue_temp.dta", replace

// ts_gmm_data.csv was historically exported by generate-gmm-fixtures.do from this same seed-12345 DGP; the M-16 re-base removed that generator's synthetic block, so the export moves here where test-cue.R's consumers live until the M-17 CUE re-base retires it. (Deliberately plain export without %21.0g to remain byte-compatible with the existing committed file -- regenerating must not churn it.)
export delimited using "`outdir'/ts_gmm_data.csv", replace

// --- HAC CUE Bartlett bw=3 ---
use "`outdir'/_ts_cue_temp.dta", clear
ivreg2 y w (x = z1 z2), cue bw(3) kernel(bartlett) robust
save_ivreg2_results, prefix(ts_cue) suffix(hac_bartlett_bw3) outdir(`outdir')


/*===========================================================================
  CUE b0 — evaluate at 2SLS estimates
  First run 2SLS to get estimates, then pass them as b0()
===========================================================================*/

use "`outdir'/_card_cue_temp.dta", clear
// First get 2SLS estimates
quietly ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust
matrix b_2sls = e(b)
// Now evaluate CUE at those estimates using b0() (b0 implies cue)
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust b0(b_2sls)
save_ivreg2_results, prefix(card_overid_b0) suffix(robust) outdir(`outdir')

// Also save the b0 vector used
quietly {
    preserve
    clear
    local ncols = colsof(b_2sls)
    local names : colnames b_2sls
    set obs `ncols'
    gen str32 term = ""
    gen double b0_value = .
    forvalues i = 1/`ncols' {
        local nm : word `i' of `names'
        replace term = "`nm'" in `i'
        replace b0_value = b_2sls[1, `i'] in `i'
    }
    export delimited using "`outdir'/card_overid_b0_vector.csv", replace
    restore
}


/*===========================================================================
  Clean up temp files
===========================================================================*/
capture erase "`outdir'/_card_cue_temp.dta"
capture erase "`outdir'/_ts_cue_temp.dta"

display "CUE fixture generation complete."
