/*===========================================================================
  generate-wmatrix-fixtures.do
  ----------------------------
  Generates CSV benchmark fixtures for wmatrix/smatrix GMM testing in ivreg2r.
  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-wmatrix-fixtures.do
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
save "`outdir'/_card_wmatrix_temp.dta", replace

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


/*---------------------------------------------------------------------------
  Helper program: save a matrix to CSV
---------------------------------------------------------------------------*/
capture program drop save_matrix
program define save_matrix
    syntax, mat(string) outfile(string)
    matrix M = `mat'
    local mr = rowsof(M)
    local mc = colsof(M)
    quietly {
        preserve
        clear
        set obs `mr'
        forvalues j = 1/`mc' {
            gen double v`j' = .
            forvalues i = 1/`mr' {
                replace v`j' = M[`i', `j'] in `i'
            }
        }
        export delimited using "`outfile'", replace
        restore
    }
end


/*===========================================================================
  Card overidentified model:
  lwage ~ exper expersq black south (educ = nearc2 nearc4)
  L = 7 instruments (nearc2 nearc4 exper expersq black south _cons)
===========================================================================*/

/*---------------------------------------------------------------------------
  Fixture 1: wmatrix = I(L), robust
  Identity weighting matrix — inefficient GMM
---------------------------------------------------------------------------*/
use "`outdir'/_card_wmatrix_temp.dta", clear
matrix I_L = I(7)
// Column/row names must match instrument list in Stata's variable order
matrix colnames I_L = exper expersq black south _cons nearc2 nearc4
matrix rownames I_L = exper expersq black south _cons nearc2 nearc4
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust wmatrix(I_L)
save_ivreg2_results, prefix(card_overid_wmat_ident) suffix(robust) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 2: wmatrix = I(L), robust small
---------------------------------------------------------------------------*/
use "`outdir'/_card_wmatrix_temp.dta", clear
matrix I_L = I(7)
matrix colnames I_L = exper expersq black south _cons nearc2 nearc4
matrix rownames I_L = exper expersq black south _cons nearc2 nearc4
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust small wmatrix(I_L)
save_ivreg2_results, prefix(card_overid_wmat_ident) suffix(robust_small) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 3: smatrix = S from robust 2-step
  First run robust GMM2S to get S, then re-estimate with that S
---------------------------------------------------------------------------*/
use "`outdir'/_card_wmatrix_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust gmm2s
matrix S_robust = e(S)
save_matrix, mat(S_robust) outfile("`outdir'/card_overid_smat_S_robust.csv")
// Now use smatrix with that S
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust smatrix(S_robust)
save_ivreg2_results, prefix(card_overid_smat) suffix(robust) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 4: smatrix = S from robust + gmm2s
  Should give same as fixture 3 (smatrix implies efficient GMM)
---------------------------------------------------------------------------*/
use "`outdir'/_card_wmatrix_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust gmm2s
matrix S_robust = e(S)
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust gmm2s smatrix(S_robust)
save_ivreg2_results, prefix(card_overid_smat_gmm2s) suffix(robust) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 5: wmatrix = I(L) + gmm2s, robust
  Identity first step → efficient second step
---------------------------------------------------------------------------*/
use "`outdir'/_card_wmatrix_temp.dta", clear
matrix I_L = I(7)
matrix colnames I_L = exper expersq black south _cons nearc2 nearc4
matrix rownames I_L = exper expersq black south _cons nearc2 nearc4
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust gmm2s wmatrix(I_L)
save_ivreg2_results, prefix(card_overid_wmat_gmm2s) suffix(robust) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 6: wmatrix = I(L) + smatrix, robust
  W for beta, user S for VCV/J
---------------------------------------------------------------------------*/
use "`outdir'/_card_wmatrix_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust gmm2s
matrix S_robust = e(S)
save_matrix, mat(S_robust) outfile("`outdir'/card_overid_wmat_smat_S_robust.csv")
matrix I_L = I(7)
matrix colnames I_L = exper expersq black south _cons nearc2 nearc4
matrix rownames I_L = exper expersq black south _cons nearc2 nearc4
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust wmatrix(I_L) smatrix(S_robust)
save_ivreg2_results, prefix(card_overid_wmat_smat) suffix(robust) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 7: wmatrix from IID (Z'Z/N), robust
  IID-based W but robust VCV
---------------------------------------------------------------------------*/
use "`outdir'/_card_wmatrix_temp.dta", clear
// Get Z'Z/N as W from the IID 2SLS run
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4)
matrix W_iid = e(S)
save_matrix, mat(W_iid) outfile("`outdir'/card_overid_wmat_from_iid_W.csv")
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), robust wmatrix(W_iid)
save_ivreg2_results, prefix(card_overid_wmat_from_iid) suffix(robust) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 8: smatrix, robust — just-identified model
  Model: lwage ~ exper expersq black south (educ = nearc4)
---------------------------------------------------------------------------*/
use "`outdir'/_card_wmatrix_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc4), robust gmm2s
matrix S_robust_ji = e(S)
save_matrix, mat(S_robust_ji) outfile("`outdir'/card_justid_smat_S_robust.csv")
ivreg2 lwage exper expersq black south (educ = nearc4), robust smatrix(S_robust_ji)
save_ivreg2_results, prefix(card_justid_smat) suffix(robust) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 9: wmatrix = I(L), cluster(age)
---------------------------------------------------------------------------*/
use "`outdir'/_card_wmatrix_temp.dta", clear
matrix I_L = I(7)
matrix colnames I_L = exper expersq black south _cons nearc2 nearc4
matrix rownames I_L = exper expersq black south _cons nearc2 nearc4
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cluster(age) wmatrix(I_L)
save_ivreg2_results, prefix(card_overid_wmat_ident) suffix(cluster) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 10: smatrix from cluster, cluster(age)
---------------------------------------------------------------------------*/
use "`outdir'/_card_wmatrix_temp.dta", clear
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cluster(age) gmm2s
matrix S_cluster = e(S)
save_matrix, mat(S_cluster) outfile("`outdir'/card_overid_smat_S_cluster.csv")
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), cluster(age) smatrix(S_cluster)
save_ivreg2_results, prefix(card_overid_smat) suffix(cluster) outdir(`outdir')


/*===========================================================================
  Clean up temp files
===========================================================================*/
capture erase "`outdir'/_card_wmatrix_temp.dta"

display "wmatrix/smatrix fixture generation complete."
