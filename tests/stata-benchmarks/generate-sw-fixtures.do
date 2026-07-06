/*===========================================================================
  generate-sw-fixtures.do
  -----------------------
  Generates CSV benchmark fixtures for Stock-Watson (2008) panel-robust VCE
  testing in ivreg2r.

  Uses wagepan data (bcuse wagepan).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-sw-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

// Load wagepan: cache-first from the source of record that also feeds the bundled data(wagepan), with a bcuse fallback (bcuse calls clear all -> must come before program define). Guarded per the generate-fixtures.do card pattern: bc.edu rate limiting can return rc 0 with empty memory.
capture use "../validation/data/wagepan.dta", clear
if _rc {
    capture bcuse wagepan, clear
}
quietly describe
if r(N) == 0 | r(k) == 0 {
    display as error "Could not load wagepan (no local cache, bcuse failed)."
    exit 601
}
save "`outdir'/_wagepan_sw_temp.dta", replace

/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV
---------------------------------------------------------------------------*/
capture program drop save_sw_results
program define save_sw_results
    syntax, prefix(string) outdir(string)

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
        export delimited using "`outdir'/`prefix'_coef.csv", replace
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
        export delimited using "`outdir'/`prefix'_vcov.csv", replace
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

        // Underidentification (KP)
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

        export delimited using "`outdir'/`prefix'_diagnostics.csv", replace
        restore
    }
end


/*===========================================================================
  SW Fixtures using wagepan data
  Just-id model: lwage exper (hours = educ)
  Overid model:  lwage exper (hours = educ married)
  Panel: nr (individual), year (time)
===========================================================================*/

// --- 1. SW just-identified ---
use "`outdir'/_wagepan_sw_temp.dta", clear
tsset nr year
ivreg2 lwage exper (hours = educ), sw
save_sw_results, prefix(wp_sw_justid) outdir(`outdir')

// --- 2. SW just-identified + small ---
use "`outdir'/_wagepan_sw_temp.dta", clear
tsset nr year
ivreg2 lwage exper (hours = educ), sw small
save_sw_results, prefix(wp_sw_justid_small) outdir(`outdir')

// --- 3. SW overidentified ---
use "`outdir'/_wagepan_sw_temp.dta", clear
tsset nr year
ivreg2 lwage exper (hours = educ married), sw
save_sw_results, prefix(wp_sw_overid) outdir(`outdir')

// --- 4. SW overidentified + small ---
use "`outdir'/_wagepan_sw_temp.dta", clear
tsset nr year
ivreg2 lwage exper (hours = educ married), sw small
save_sw_results, prefix(wp_sw_overid_small) outdir(`outdir')

// --- 5. SW + aweight ---
use "`outdir'/_wagepan_sw_temp.dta", clear
tsset nr year
gen aw = hours + 10
ivreg2 lwage exper (hours = educ married) [aw=aw], sw
save_sw_results, prefix(wp_sw_aweight) outdir(`outdir')

// --- 6. SW + LIML ---
use "`outdir'/_wagepan_sw_temp.dta", clear
tsset nr year
ivreg2 lwage exper (hours = educ married), sw liml
save_sw_results, prefix(wp_sw_liml) outdir(`outdir')

// --- 7. SW + dofminus(1) ---
use "`outdir'/_wagepan_sw_temp.dta", clear
tsset nr year
ivreg2 lwage exper (hours = educ married), sw dofminus(1)
save_sw_results, prefix(wp_sw_dofminus) outdir(`outdir')

// --- 8. SW + endogeneity test ---
use "`outdir'/_wagepan_sw_temp.dta", clear
tsset nr year
ivreg2 lwage exper (hours = educ married), sw endog(hours)
save_sw_results, prefix(wp_sw_endog) outdir(`outdir')

// --- 9. SW + center ---
use "`outdir'/_wagepan_sw_temp.dta", clear
tsset nr year
ivreg2 lwage exper (hours = educ married), sw center
save_sw_results, prefix(wp_sw_center) outdir(`outdir')

// --- 10. SW + gmm2s ---
use "`outdir'/_wagepan_sw_temp.dta", clear
tsset nr year
ivreg2 lwage exper (hours = educ married), sw gmm2s
save_sw_results, prefix(wp_sw_gmm2s) outdir(`outdir')

// --- 11. SW + cue ---
use "`outdir'/_wagepan_sw_temp.dta", clear
tsset nr year
ivreg2 lwage exper (hours = educ married), sw cue
save_sw_results, prefix(wp_sw_cue) outdir(`outdir')

// Clean up temp file
capture erase "`outdir'/_wagepan_sw_temp.dta"

display ""
display "=== SW fixture generation complete ==="
display ""
