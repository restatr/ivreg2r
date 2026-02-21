/*===========================================================================
  generate-hac-fixtures.do
  ------------------------
  Generates CSV benchmark fixtures for HAC/AC VCE testing in ivreg2r.
  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-hac-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV
  (Identical to the helper in generate-fixtures.do)
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
  Generate simulated time-series IV data
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

// Save dataset
export delimited using "`outdir'/ts_hac_data.csv", replace
save "`outdir'/_ts_hac_temp.dta", replace


/*---------------------------------------------------------------------------
  HAC fixtures: overidentified model (y x w (x = z1 z2))
---------------------------------------------------------------------------*/

// --- HAC Bartlett bw=3 ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(3) kernel(bartlett) robust
save_ivreg2_results, prefix(ts_hac) suffix(bartlett_bw3) outdir(`outdir')

// --- HAC Bartlett bw=5 ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(5) kernel(bartlett) robust
save_ivreg2_results, prefix(ts_hac) suffix(bartlett_bw5) outdir(`outdir')

// --- HAC Parzen bw=4 ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(4) kernel(parzen) robust
save_ivreg2_results, prefix(ts_hac) suffix(parzen_bw4) outdir(`outdir')

// --- HAC Truncated bw=3 ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(3) kernel(truncated) robust
save_ivreg2_results, prefix(ts_hac) suffix(truncated_bw3) outdir(`outdir')

// --- HAC Tukey-Hanning bw=3 ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(3) kernel(thann) robust
save_ivreg2_results, prefix(ts_hac) suffix(thanning_bw3) outdir(`outdir')

// --- HAC Tukey-Hamming bw=3 ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(3) kernel(thamm) robust
save_ivreg2_results, prefix(ts_hac) suffix(thamming_bw3) outdir(`outdir')

// --- HAC Daniell bw=4 ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(4) kernel(dan) robust
save_ivreg2_results, prefix(ts_hac) suffix(daniell_bw4) outdir(`outdir')

// --- HAC Tent bw=4 ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(4) kernel(tent) robust
save_ivreg2_results, prefix(ts_hac) suffix(tent_bw4) outdir(`outdir')

// --- HAC QS bw=3 (just-identified: y (x = z1)) ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1), bw(3) kernel(qs) robust
save_ivreg2_results, prefix(ts_hac) suffix(qs_bw3_justid) outdir(`outdir')

// --- HAC QS bw=3 (overidentified: y w (x = z1 z2)) ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(3) kernel(qs) robust
save_ivreg2_results, prefix(ts_hac) suffix(qs_bw3) outdir(`outdir')


/*---------------------------------------------------------------------------
  AC fixtures: overidentified model
---------------------------------------------------------------------------*/

// --- AC Bartlett bw=3 ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(3) kernel(bartlett)
save_ivreg2_results, prefix(ts_ac) suffix(bartlett_bw3) outdir(`outdir')

// --- AC Bartlett bw=5 (just-identified) ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1), bw(5) kernel(bartlett)
save_ivreg2_results, prefix(ts_ac) suffix(bartlett_bw5_justid) outdir(`outdir')

// --- AC Parzen bw=4 ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(4) kernel(parzen)
save_ivreg2_results, prefix(ts_ac) suffix(parzen_bw4) outdir(`outdir')


/*---------------------------------------------------------------------------
  Gappy time-series HAC fixture
---------------------------------------------------------------------------*/
use "`outdir'/_ts_hac_temp.dta", clear
// Drop some observations to create gaps
drop if t == 50 | t == 100 | t == 150
tsset t
ivreg2 y w (x = z1 z2), bw(3) kernel(bartlett) robust
save_ivreg2_results, prefix(ts_gap_hac) suffix(bartlett_bw3) outdir(`outdir')

// Save gappy dataset
export delimited using "`outdir'/ts_gap_hac_data.csv", replace


/*---------------------------------------------------------------------------
  Helper: save auto-bw results (extends save_ivreg2_results with e(bw))
---------------------------------------------------------------------------*/
capture program drop save_autobiw_results
program define save_autobiw_results
    syntax, prefix(string) suffix(string) outdir(string)

    // Save standard results via existing helper
    save_ivreg2_results, prefix(`prefix') suffix(`suffix') outdir(`outdir')

    // Save bandwidth value separately
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
  Auto-bandwidth fixtures: HAC overidentified (y w (x = z1 z2))
---------------------------------------------------------------------------*/

// --- HAC auto Bartlett ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(auto) kernel(bartlett) robust
save_autobiw_results, prefix(ts_hac) suffix(auto_bartlett) outdir(`outdir')

// --- HAC auto Parzen ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(auto) kernel(parzen) robust
save_autobiw_results, prefix(ts_hac) suffix(auto_parzen) outdir(`outdir')

// --- HAC auto QS ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(auto) kernel(qs) robust
save_autobiw_results, prefix(ts_hac) suffix(auto_qs) outdir(`outdir')

// --- AC auto Bartlett (iid + kernel) ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1 z2), bw(auto) kernel(bartlett)
save_autobiw_results, prefix(ts_ac) suffix(auto_bartlett) outdir(`outdir')

// --- HAC auto Bartlett just-identified ---
use "`outdir'/_ts_hac_temp.dta", clear
ivreg2 y w (x = z1), bw(auto) kernel(bartlett) robust
save_autobiw_results, prefix(ts_hac) suffix(auto_bartlett_justid) outdir(`outdir')


/*---------------------------------------------------------------------------
  Auto-bandwidth fixture: gappy data
---------------------------------------------------------------------------*/
use "`outdir'/_ts_hac_temp.dta", clear
drop if t == 50 | t == 100 | t == 150
tsset t
ivreg2 y w (x = z1 z2), bw(auto) kernel(bartlett) robust
save_autobiw_results, prefix(ts_gap_hac) suffix(auto_bartlett) outdir(`outdir')


/*---------------------------------------------------------------------------
  Clean up temp files
---------------------------------------------------------------------------*/
capture erase "`outdir'/_ts_hac_temp.dta"

display "HAC/AC fixture generation complete."
