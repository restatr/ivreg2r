/*===========================================================================
  generate-center-fixtures.do
  --------------------------
  Generates CSV benchmark fixtures for center option testing in ivreg2r.
  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-center-fixtures.do
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
save "`outdir'/_card_center_temp.dta", replace

/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV
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
  Center fixtures using Card data (overidentified model)
  Model: lwage ~ exper expersq black south (educ = nearc2 nearc4)
===========================================================================*/

// --- 1. HC0 + center ---
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust center
save_ivreg2_results, prefix(card_overid) suffix(center_hc0) outdir(`outdir')

// --- 2. HC1 + small + center ---
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust small center
save_ivreg2_results, prefix(card_overid) suffix(center_hc1_small) outdir(`outdir')

// --- 3. Cluster + center ---
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cluster(smsa) center
save_ivreg2_results, prefix(card_overid) suffix(center_cl) outdir(`outdir')

// --- 4. Cluster + small + center ---
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cluster(smsa) small center
save_ivreg2_results, prefix(card_overid) suffix(center_cl_small) outdir(`outdir')

// --- 5. Just-identified HC0 + center ---
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc4), robust center
save_ivreg2_results, prefix(card_justid) suffix(center_hc0) outdir(`outdir')

// --- 6. GMM2S robust + center ---
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), gmm2s robust center
save_ivreg2_results, prefix(card_overid) suffix(center_gmm2s_hc0) outdir(`outdir')

// --- 7. GMM2S cluster(age) + center ---
// (cluster(smsa) causes rank-deficient omega; use age instead)
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), gmm2s cluster(age) center
save_ivreg2_results, prefix(card_overid) suffix(center_gmm2s_cl) outdir(`outdir')

// --- 8. CUE robust + center ---
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cue robust center
save_ivreg2_results, prefix(card_overid) suffix(center_cue_hc0) outdir(`outdir')

// --- 9. CUE cluster(age) + center ---
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cue cluster(age) center
save_ivreg2_results, prefix(card_overid) suffix(center_cue_cl) outdir(`outdir')

// --- 10. Endogeneity test + center ---
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust center endog(educ)
save_ivreg2_results, prefix(card_overid) suffix(center_endog) outdir(`outdir')

// --- 11. Orthogonality test + center ---
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust center orthog(nearc2)
save_ivreg2_results, prefix(card_overid) suffix(center_orthog) outdir(`outdir')

// --- 12. dofminus + center ---
use "`outdir'/_card_center_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust center dofminus(1)
save_ivreg2_results, prefix(card_overid) suffix(center_dofminus) outdir(`outdir')


/*===========================================================================
  HAC + center fixture (time-series Card data, Panel sim not needed since
  we use sorted card obs as a pseudo-time-series)
===========================================================================*/

// --- 13. HAC Bartlett + center ---
use "`outdir'/_card_center_temp.dta", clear
gen t = _n
tsset t
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), bw(3) kernel(bartlett) center
save_ivreg2_results, prefix(card_overid) suffix(center_hac_bartlett) outdir(`outdir')


/*===========================================================================
  Cleanup
===========================================================================*/
capture erase "`outdir'/_card_center_temp.dta"

display as result "Center fixtures generated successfully."
