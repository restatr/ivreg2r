/*===========================================================================
  generate-ck-fixtures.do
  -----------------------
  Generates CSV benchmark fixtures for Kiefer, Driscoll-Kraay, and
  Thompson (cluster+kernel) VCE testing in ivreg2r.
  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Uses wagepan data (bcuse wagepan).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-ck-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

/*---------------------------------------------------------------------------
  Load wagepan data (bcuse calls clear all — must come before program define)
---------------------------------------------------------------------------*/
bcuse wagepan, clear
save "`outdir'/_wp_temp.dta", replace

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
  Export wagepan data as CSV (for R tests)
---------------------------------------------------------------------------*/
use "`outdir'/_wp_temp.dta", clear
export delimited using "`outdir'/wp_ck_data.csv", replace


/*---------------------------------------------------------------------------
  Kiefer fixtures
  Kiefer = AC with kernel=truncated, bw=T (time span)
  Panel: nr (individual), year (time)
  Model overid: lwage exper expersq married union (hours = educ black)
  Model justid: lwage exper expersq married union (hours = educ)
---------------------------------------------------------------------------*/

// --- Kiefer overid ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ black), kiefer
save_ivreg2_results, prefix(wp_kiefer) suffix(overid) outdir(`outdir')

// --- Kiefer overid small ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ black), kiefer small
save_ivreg2_results, prefix(wp_kiefer) suffix(overid_small) outdir(`outdir')

// --- Kiefer justid ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ), kiefer
save_ivreg2_results, prefix(wp_kiefer) suffix(justid) outdir(`outdir')


/*---------------------------------------------------------------------------
  Driscoll-Kraay fixtures
  dkraay(bw) clusters on tvar + kernel smoothing across time lags
---------------------------------------------------------------------------*/

// --- DK Bartlett bw=3 overid ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ black), dkraay(3)
save_ivreg2_results, prefix(wp_dk) suffix(bar3_overid) outdir(`outdir')

// --- DK Bartlett bw=3 overid small ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ black), dkraay(3) small
save_ivreg2_results, prefix(wp_dk) suffix(bar3_overid_small) outdir(`outdir')

// --- DK Parzen bw=4 overid ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ black), dkraay(4) kernel(parzen)
save_ivreg2_results, prefix(wp_dk) suffix(par4_overid) outdir(`outdir')

// --- DK Truncated bw=2 overid ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ black), dkraay(2) kernel(truncated)
save_ivreg2_results, prefix(wp_dk) suffix(tru2_overid) outdir(`outdir')

// --- DK Bartlett bw=3 justid ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ), dkraay(3)
save_ivreg2_results, prefix(wp_dk) suffix(bar3_justid) outdir(`outdir')


/*---------------------------------------------------------------------------
  Thompson (two-way cluster + kernel) fixtures
  cluster(ivar tvar) kernel(bar) bw(3)
---------------------------------------------------------------------------*/

// --- Thompson Bartlett bw=3 overid ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ black), cluster(nr year) kernel(bar) bw(3)
save_ivreg2_results, prefix(wp_ck) suffix(bar3_overid) outdir(`outdir')

// --- Thompson Bartlett bw=3 overid small ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ black), cluster(nr year) kernel(bar) bw(3) small
save_ivreg2_results, prefix(wp_ck) suffix(bar3_overid_small) outdir(`outdir')

// --- Thompson Parzen bw=4 overid ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ black), cluster(nr year) kernel(parzen) bw(4)
save_ivreg2_results, prefix(wp_ck) suffix(par4_overid) outdir(`outdir')

// --- Thompson Bartlett bw=3 justid ---
use "`outdir'/_wp_temp.dta", clear
xtset nr year
ivreg2 lwage exper expersq married union (hours = educ), cluster(nr year) kernel(bar) bw(3)
save_ivreg2_results, prefix(wp_ck) suffix(bar3_justid) outdir(`outdir')


/*---------------------------------------------------------------------------
  Clean up temp files
---------------------------------------------------------------------------*/
capture erase "`outdir'/_wp_temp.dta"

display "Cluster+kernel fixture generation complete."
