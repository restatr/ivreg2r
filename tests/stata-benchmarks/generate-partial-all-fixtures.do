/*===========================================================================
  generate-partial-all-fixtures.do
  --------------------------------
  Generates CSV benchmark fixtures for the F5 bug fix: identification
  tests with an empty exogenous block (partial(_all) with K1 > 0, and
  noconstant models with no included exogenous regressors).
  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Data: mroz_data.csv — the same CSV the R tests read (753 rows incl.
  missing lwage; ivreg2 drops those rows -> N = 428).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-partial-all-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV
  (same schema as generate-partial-fixtures.do)
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

/*---------------------------------------------------------------------------
  Helper program: save redundancy diagnostics to CSV
  (same schema as generate-redundancy-fixtures.do)
---------------------------------------------------------------------------*/
capture program drop save_redundancy_results
program define save_redundancy_results
    syntax, prefix(string) suffix(string) outdir(string)

    quietly {
        preserve
        clear
        set obs 1

        gen double redstat = .
        gen double redp = .
        gen double reddf = .
        gen str80 redlist = ""

        capture replace redstat = e(redstat)
        capture replace redp = e(redp)
        capture replace reddf = e(reddf)
        capture replace redlist = "`e(redlist)'"

        export delimited using ///
            "`outdir'/`prefix'_redundancy_`suffix'.csv", replace
        restore
    }
end

/*===========================================================================
  Load data
===========================================================================*/
import delimited "`outdir'/mroz_data.csv", clear case(preserve)

/*===========================================================================
  1. partial(_all) with K1 > 0 — F5 ticket repro, IID
===========================================================================*/
display "--- 1. partial(_all), IID ---"
ivreg2 lwage exper expersq (educ = age kidslt6 kidsge6), partial(_all)
save_ivreg2_results, prefix(mroz_partial_all) suffix(iid) outdir(`outdir')

/*===========================================================================
  2. partial(_all) with K1 > 0 — robust (KP rk LM / rk Wald F)
===========================================================================*/
display "--- 2. partial(_all), robust ---"
ivreg2 lwage exper expersq (educ = age kidslt6 kidsge6), partial(_all) robust
save_ivreg2_results, prefix(mroz_partial_all) suffix(robust) outdir(`outdir')

/*===========================================================================
  3. noconstant with no included exogenous regressors (K2 = 0, no
     partialling) — same empty-exog ranktest path arrived at directly
===========================================================================*/
display "--- 3. noconstant, K2 = 0, IID ---"
ivreg2 lwage (educ = age kidslt6 kidsge6), noconstant
save_ivreg2_results, prefix(mroz_nocons) suffix(iid) outdir(`outdir')

/*===========================================================================
  4. Weighted partial(_all) — the guard's weighted entry path
===========================================================================*/
gen double wtvar = hours + 1

display "--- 4a. weighted partial(_all), IID ---"
ivreg2 lwage exper expersq (educ = age kidslt6 kidsge6) [aw=wtvar], partial(_all)
save_ivreg2_results, prefix(mroz_partial_all_w) suffix(iid) outdir(`outdir')

display "--- 4b. weighted partial(_all), robust ---"
ivreg2 lwage exper expersq (educ = age kidslt6 kidsge6) [aw=wtvar], partial(_all) robust
save_ivreg2_results, prefix(mroz_partial_all_w) suffix(robust) outdir(`outdir')

/*===========================================================================
  5. Redundancy test with an empty partial set: noconstant, K2 = 0, ALL
     excluded instruments tested (ranktest gets an empty partial varlist)
===========================================================================*/
display "--- 5a. noconstant, redundant(all instruments), IID ---"
ivreg2 lwage (educ = age kidslt6 kidsge6), noconstant redundant(age kidslt6 kidsge6)
save_redundancy_results, prefix(mroz_nocons) suffix(iid) outdir(`outdir')

display "--- 5b. noconstant, redundant(all instruments), robust ---"
ivreg2 lwage (educ = age kidslt6 kidsge6), noconstant redundant(age kidslt6 kidsge6) robust
save_redundancy_results, prefix(mroz_nocons) suffix(robust) outdir(`outdir')

display "Done: partial-all fixtures written to `outdir'"
