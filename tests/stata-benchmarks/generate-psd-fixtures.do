/*===========================================================================
  generate-psd-fixtures.do
  ------------------------
  Generates CSV benchmark fixtures for the psd0/psda correction-site fix
  (ticket F2). The binding configuration is Driscoll-Kraay + truncated
  kernel on the wagepan panel, whose moment covariance S is indefinite, so
  the correction actually fires and the fixtures pin down BOTH the
  correction site (S, not the final VCV) and the correction space (the
  L x L instrument-space S, not the K x K meat — the model is
  overidentified, so psd(A'SA) != A'psd(S)A is observable).

  Each block also logs the minimum eigenvalue of e(S) ("MINEIG <label>:")
  so binding status is documented in the .log.

  Data: wagepan loads cache-first from ../validation/data/wagepan.dta, the file of record from which pkg/data-raw/wagepan.R builds the bundled data(wagepan) bit-identically, so Stata and the R tests still compute on identical inputs, now at full precision (wp_ck_data.csv was a truncated 8-digit export and was retired 2026-07-06; the wp_psd fixtures were verified byte-identical under the full-precision load at that retirement). The load is hoisted above the program define blocks so its bcuse fallback (which calls clear all) is safe. phillips loads cache-first from ../validation/data/phillips.dta (its constructed lag columns ship inside the .dta), re-pointed off the retired phillips_data.csv snapshot at the M-27 re-base 2026-07-06.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-psd-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

local outdir "tests/stata-benchmarks/fixtures"

/*---------------------------------------------------------------------------
  Load wagepan: cache-first from the source of record that also feeds the bundled data(wagepan), with a bcuse fallback — hoisted above the program define blocks so the bcuse fallback (which calls clear all) is safe. Guarded per the generate-fixtures.do card pattern: bc.edu rate limiting can return rc 0 with empty memory.
---------------------------------------------------------------------------*/
capture use "../validation/data/wagepan.dta", clear
if _rc {
    capture bcuse wagepan, clear
}
quietly describe
if r(N) == 0 | r(k) == 0 {
    display as error "Could not load wagepan (no local cache, bcuse failed)."
    exit 601
}
save "`outdir'/_wp_psd_temp.dta", replace

/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV
  (same schema as generate-ck-fixtures.do)
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

    // --- e(S) matrix ---
    matrix S = e(S)
    quietly {
        preserve
        clear
        local sr = rowsof(S)
        set obs `sr'
        forvalues j = 1/`sr' {
            gen double v`j' = .
            forvalues i = 1/`sr' {
                replace v`j' = S[`i', `j'] in `i'
            }
        }
        export delimited using "`outdir'/`prefix'_eS_`suffix'.csv", replace
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
  Helper: save e(S) column names to a one-column CSV (same as
  generate-stored-matrices-fixtures.do). The R tests map "_cons" ->
  "(Intercept)" and align by name.
---------------------------------------------------------------------------*/
capture program drop save_s_names
program define save_s_names
    syntax, outfile(string)
    matrix Smat = e(S)
    local cn : colnames Smat
    tempname fh
    file open `fh' using "`outfile'", write replace
    file write `fh' "name" _n
    foreach c of local cn {
        file write `fh' "`c'" _n
    }
    file close `fh'
end

/*---------------------------------------------------------------------------
  Helper program: log the minimum eigenvalue of e(S)
  (symeigen returns eigenvalues in descending order -> min is the last)
---------------------------------------------------------------------------*/
capture program drop report_mineig
program define report_mineig
    syntax, label(string)
    matrix Smat = e(S)
    matrix symeigen EVEC eval = Smat
    local mc = colsof(eval)
    local mineig = eval[1, `mc']
    display as text "MINEIG `label': `mineig'"
end

/*===========================================================================
  Driscoll-Kraay + truncated kernel (the binding repro)
  Overid model: lwage exper expersq married union (hours = educ black)
  Panel: nr (individual), year (time)
===========================================================================*/
// --- Baseline (no psd): document binding via min eigenvalue of e(S) ---
use "`outdir'/_wp_psd_temp.dta", clear
tsset nr year
ivreg2 lwage exper expersq married union (hours = educ black), dkraay(2) kernel(truncated)
report_mineig, label(dk_tru2_nopsd)

// --- DK truncated psd0 ---
use "`outdir'/_wp_psd_temp.dta", clear
tsset nr year
ivreg2 lwage exper expersq married union (hours = educ black), dkraay(2) kernel(truncated) psd0
report_mineig, label(dk_tru2_psd0)
save_ivreg2_results, prefix(wp_psd) suffix(dk_tru2_psd0) outdir(`outdir')
save_s_names, outfile("`outdir'/wp_psd_inames.csv")

// --- DK truncated psda ---
use "`outdir'/_wp_psd_temp.dta", clear
tsset nr year
ivreg2 lwage exper expersq married union (hours = educ black), dkraay(2) kernel(truncated) psda
report_mineig, label(dk_tru2_psda)
save_ivreg2_results, prefix(wp_psd) suffix(dk_tru2_psda) outdir(`outdir')

// --- DK truncated psd0 small ---
use "`outdir'/_wp_psd_temp.dta", clear
tsset nr year
ivreg2 lwage exper expersq married union (hours = educ black), dkraay(2) kernel(truncated) psd0 small
save_ivreg2_results, prefix(wp_psd) suffix(dk_tru2_psd0_small) outdir(`outdir')

/*===========================================================================
  Stock-Watson + psda (regression fixture for the SW correction-space fix;
  psda is the Stock-Watson 2008 Remark 8 correction)
===========================================================================*/
use "`outdir'/_wp_psd_temp.dta", clear
tsset nr year
ivreg2 lwage exper expersq married union (hours = educ black), sw psda
report_mineig, label(sw_psda)
save_ivreg2_results, prefix(wp_psd) suffix(sw_psda) outdir(`outdir')

/*===========================================================================
  Two-way cluster + psd0 (CGM subtraction is a classic indefiniteness
  source; saved regardless of binding as regression coverage)
===========================================================================*/
use "`outdir'/_wp_psd_temp.dta", clear
tsset nr year
ivreg2 lwage exper expersq married union (hours = educ black), cluster(nr year) psd0
report_mineig, label(twoway_psd0)
save_ivreg2_results, prefix(wp_psd) suffix(twoway_psd0) outdir(`outdir')

/*===========================================================================
  HAC truncated kernel + psd0 probe (overidentified time series, phillips)
===========================================================================*/
capture use "../validation/data/phillips.dta", clear
if _rc {
    use http://fmwww.bc.edu/ec-p/data/wooldridge/phillips.dta, clear
}
tsset year

ivreg2 cinf (unem = unem_1 unem_2), robust kernel(tru) bw(2) psd0
report_mineig, label(phil_hac_tru2_psd0)
save_ivreg2_results, prefix(phil_psd) suffix(hac_tru2_psd0) outdir(`outdir')
save_s_names, outfile("`outdir'/phil_psd_inames.csv")

// Clean up temp file
capture erase "`outdir'/_wp_psd_temp.dta"

display "generate-psd-fixtures.do complete"
