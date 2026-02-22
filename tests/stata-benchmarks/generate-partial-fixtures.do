/*===========================================================================
  generate-partial-fixtures.do
  ----------------------------
  Generates CSV benchmark fixtures for partial option testing in ivreg2r.
  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-partial-fixtures.do
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
save "`outdir'/_card_partial_temp.dta", replace

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

        // Underidentification
        gen double underid_stat = e(idstat)
        gen double underid_p = e(idp)
        gen int underid_df = e(iddf)

        // Weak identification
        gen double weak_id_cd_f = e(cdf)
        gen double weak_id_kp_f = e(rkf)

        // Anderson-Rubin
        gen double ar_f = e(arf)
        gen double ar_fp = e(arfp)
        gen double ar_chi2 = e(archi2)
        gen double ar_chi2p = e(archi2p)
        gen int ar_df = e(ardf)
        gen int ar_df_r = e(ardf_r)

        // Endogeneity
        gen double endog_stat = e(estat)
        gen double endog_p = e(estatp)
        gen int endog_df = e(estatdf)

        // Model F
        gen double model_f = e(F)
        gen double model_f_p = e(Fp)
        gen int model_f_df1 = e(Fdf1)
        gen int model_f_df2 = e(Fdf2)

        // Fit stats
        gen double sigma = e(rmse)
        gen double rss = e(rss)
        gen double r2 = e(r2)
        gen double r2_a = e(r2_a)
        gen double r2u = e(r2u)
        gen double r2c = e(r2c)
        gen double mss = e(mss)
        gen int N = e(N)
        gen int K = e(rankxx)
        gen int L = e(rankzz)
        gen int df_r = e(df_r)
        gen int sdofminus = e(sdofminus)
        gen int partial_ct = e(partial_ct)
        gen int partialcons = e(partialcons)

        export delimited using "`outdir'/`prefix'_diagnostics_`suffix'.csv", replace
        restore
    }
end

/*===========================================================================
  Load data
===========================================================================*/
use "`outdir'/_card_partial_temp.dta", clear

/*===========================================================================
  1. Basic partial: partial(black south smsa)
     IV regression with 3 vars partialled out + constant
===========================================================================*/
display "--- 1a. Basic partial, IID ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(black south smsa)
save_ivreg2_results, prefix(card_partial_basic) suffix(iid) outdir(`outdir')

display "--- 1b. Basic partial, IID small ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(black south smsa) small
save_ivreg2_results, prefix(card_partial_basic) suffix(iid_small) outdir(`outdir')

/*===========================================================================
  2. partial(_cons) — demean only
===========================================================================*/
display "--- 2. partial(_cons), IID ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(_cons)
save_ivreg2_results, prefix(card_partial_cons) suffix(iid) outdir(`outdir')

/*===========================================================================
  3. partial(_all) — partial all exogenous
===========================================================================*/
display "--- 3. partial(_all), IID ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(_all)
save_ivreg2_results, prefix(card_partial_all) suffix(iid) outdir(`outdir')

/*===========================================================================
  4. Weighted partial
===========================================================================*/
display "--- 4. Weighted partial, IID ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4) [aw=weight], partial(black south smsa)
save_ivreg2_results, prefix(card_partial_weighted) suffix(iid) outdir(`outdir')

/*===========================================================================
  5. nopartialsmall
===========================================================================*/
display "--- 5a. nopartialsmall, IID ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(black south smsa) nopartialsmall
save_ivreg2_results, prefix(card_partial_nosmall) suffix(iid) outdir(`outdir')

display "--- 5b. nopartialsmall, IID small ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(black south smsa) nopartialsmall small
save_ivreg2_results, prefix(card_partial_nosmall) suffix(iid_small) outdir(`outdir')

/*===========================================================================
  6. OLS with partial
===========================================================================*/
display "--- 6. OLS partial, IID ---"
ivreg2 lwage exper expersq black south smsa, partial(black south)
save_ivreg2_results, prefix(card_partial_ols) suffix(iid) outdir(`outdir')

/*===========================================================================
  7. LIML + partial
===========================================================================*/
display "--- 7. LIML + partial, IID ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(black south smsa) liml
save_ivreg2_results, prefix(card_partial_liml) suffix(iid) outdir(`outdir')

/*===========================================================================
  8. GMM2S + partial
===========================================================================*/
display "--- 8. GMM2S + partial, robust ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(black south smsa) gmm2s robust
save_ivreg2_results, prefix(card_partial_gmm2s) suffix(robust) outdir(`outdir')

/*===========================================================================
  9. Robust + partial
===========================================================================*/
display "--- 9. Robust + partial ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(black south smsa) robust
save_ivreg2_results, prefix(card_partial_robust) suffix(hc0) outdir(`outdir')

display "--- 9b. Robust small + partial ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(black south smsa) robust small
save_ivreg2_results, prefix(card_partial_robust) suffix(hc1_small) outdir(`outdir')

/*===========================================================================
  10. Cluster + partial
===========================================================================*/
display "--- 10. Cluster + partial ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(black south smsa) cluster(smsa66)
save_ivreg2_results, prefix(card_partial_cluster) suffix(cl) outdir(`outdir')

display "--- 10b. Cluster small + partial ---"
ivreg2 lwage exper expersq black south smsa (educ = nearc2 nearc4), partial(black south smsa) cluster(smsa66) small
save_ivreg2_results, prefix(card_partial_cluster) suffix(cl_small) outdir(`outdir')

/*===========================================================================
  Cleanup
===========================================================================*/
capture erase "`outdir'/_card_partial_temp.dta"

display "=== All partial fixtures generated ==="
