/*===========================================================================
  generate-wmatrix-autobiw-fixtures.do
  ------------------------------------
  Generates CSV benchmark fixtures for wmatrix + bw(auto) HAC testing.
  Exercises the code path where auto-bandwidth uses W-step residuals.

  Uses the same time-series DGP as generate-hac-fixtures.do (seed 12345).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-wmatrix-autobiw-fixtures.do
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

        // Counts
        gen int N = e(N)
        gen int K = `K'
        gen int L = `L'

        export delimited using "`outdir'/`prefix'_diagnostics_`suffix'.csv", replace
        restore
    }
end

/*---------------------------------------------------------------------------
  Helper: save auto-bw + standard results
---------------------------------------------------------------------------*/
capture program drop save_autobiw_results
program define save_autobiw_results
    syntax, prefix(string) suffix(string) outdir(string)

    // Save standard results
    save_ivreg2_results, prefix(`prefix') suffix(`suffix') outdir(`outdir')

    // Save bandwidth value
    quietly {
        preserve
        clear
        set obs 1
        gen double bw = e(bw)
        export delimited using "`outdir'/`prefix'_bw_`suffix'.csv", replace
        restore
    }
end


/*---------------------------------------------------------------------------
  Generate time-series IV data (identical to generate-hac-fixtures.do)
---------------------------------------------------------------------------*/
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

save "`outdir'/_ts_wmat_autobiw_temp.dta", replace


/*===========================================================================
  Overidentified model: y w (x = z1 z2)
  L = 4 instruments (w _cons z1 z2)
===========================================================================*/

/*---------------------------------------------------------------------------
  Fixture 1: wmatrix = I(L), bw(auto), Bartlett, robust
  Inefficient GMM with auto-bandwidth from W-step residuals
---------------------------------------------------------------------------*/
use "`outdir'/_ts_wmat_autobiw_temp.dta", clear
matrix I_L = I(4)
matrix colnames I_L = w _cons z1 z2
matrix rownames I_L = w _cons z1 z2
ivreg2 y w (x = z1 z2), robust wmatrix(I_L) bw(auto) kernel(bartlett)
save_autobiw_results, prefix(ts_wmat_autobiw) suffix(ident_bartlett) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 2: wmatrix = I(L) + gmm2s, bw(auto), Bartlett, robust
  W as first step, then efficient second step, auto-bw from W-step residuals
---------------------------------------------------------------------------*/
use "`outdir'/_ts_wmat_autobiw_temp.dta", clear
matrix I_L = I(4)
matrix colnames I_L = w _cons z1 z2
matrix rownames I_L = w _cons z1 z2
ivreg2 y w (x = z1 z2), robust gmm2s wmatrix(I_L) bw(auto) kernel(bartlett)
save_autobiw_results, prefix(ts_wmat_autobiw) suffix(gmm2s_bartlett) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 3: wmatrix = I(L), bw(auto), Parzen, robust
---------------------------------------------------------------------------*/
use "`outdir'/_ts_wmat_autobiw_temp.dta", clear
matrix I_L = I(4)
matrix colnames I_L = w _cons z1 z2
matrix rownames I_L = w _cons z1 z2
ivreg2 y w (x = z1 z2), robust wmatrix(I_L) bw(auto) kernel(parzen)
save_autobiw_results, prefix(ts_wmat_autobiw) suffix(ident_parzen) outdir(`outdir')

/*---------------------------------------------------------------------------
  Fixture 4: wmatrix = I(L) + gmm2s, bw(auto), Parzen, robust
---------------------------------------------------------------------------*/
use "`outdir'/_ts_wmat_autobiw_temp.dta", clear
matrix I_L = I(4)
matrix colnames I_L = w _cons z1 z2
matrix rownames I_L = w _cons z1 z2
ivreg2 y w (x = z1 z2), robust gmm2s wmatrix(I_L) bw(auto) kernel(parzen)
save_autobiw_results, prefix(ts_wmat_autobiw) suffix(gmm2s_parzen) outdir(`outdir')


/*===========================================================================
  Clean up
===========================================================================*/
capture erase "`outdir'/_ts_wmat_autobiw_temp.dta"

display "wmatrix + auto-bw fixture generation complete."
