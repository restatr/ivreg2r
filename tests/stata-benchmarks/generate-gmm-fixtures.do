/*===========================================================================
  generate-gmm-fixtures.do
  ------------------------
  Generates CSV benchmark fixtures for GMM2S estimation testing (M-16 re-base): D5a option-variations layered on the griliches76 H06, abdata H88, and phillips GMM+HAC canonical bases. The canonical cells themselves live elsewhere: hf_gril H06/H25/H26 in generate-helpfile-fixtures.do, and tsop_ab_gmm2s_cl / tsop_phil in generate-ts-operator-fixtures.do.

  gmm2s+iid and just-identified gmm2s cells retired: GMM2S with an IID weighting matrix equals 2SLS algebraically and just-identified GMM2S equals 2SLS for any weighting matrix; both identities are asserted fixture-free in test-gmm2s.R, Stata's own agreement was verified at ~3e-11 in the retired fixtures, and 2SLS Stata parity is owned by the M-10 family.

  gmm2s robust cell retired to hf_gril H06 (compare_hf asserts coef/SE/VCV and every diagnostics column); gmm2s+orthog cells retired to hf H25/H26 (M-22 ruling); gmm2s cluster cell retired to tsop_ab H88; gmm2s+HAC cells retired to tsop_phil H84/H85 plus the M-19 autobandwidth gmm2s cells.

  pweight cell retired by invariance: Stata `[aw=w] ... gmm2s robust` and `[pw=w] ... gmm2s` produced byte-identical coef/vcov/diagnostics CSVs in the retired card fixtures (verified 2026-07-05); R fits the pweight variant against the aweight fixture and asserts direct equality with the aweight fit (M-11 rule).

  The synthetic seed-12345 time-series DGP moved to generate-cue-fixtures.do, which still consumes it for the ts_cue cells (M-17 will retire both).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-gmm-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"


/*===========================================================================
  Load all three datasets FIRST, before the program is defined (data loads can interact with programs, so mirror the orthog generator's ordering: every use/save block ahead of "program define").
===========================================================================*/

// --- griliches76: local cache, then BC.edu web fallback (per the H06 base model) ---
capture use "../validation/data/griliches76.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/hayashi/griliches76.dta, clear
    if _rc {
        display as error "Could not load griliches76 dataset (no local cache, no network)."
        exit 601
    }
}
quietly xi i.year
gen double awt = mod(age, 5) + 1
save "`outdir'/_gmm_gril_temp.dta", replace

// --- abdata: local cache, then BC.edu web fallback (per the H88 base model) ---
capture use "../validation/data/abdata.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/macro/abdata.dta, clear
    if _rc {
        display as error "Could not load abdata dataset (no local cache, no network)."
        exit 601
    }
}
quietly tsset id year
gen double abwt = mod(year, 3) + 1
save "`outdir'/_gmm_ab_temp.dta", replace

// --- phillips: local cache, then BC.edu web fallback (per the H83-H85 base model) ---
capture use "../validation/data/phillips.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/macro/phillips.dta, clear
    if _rc {
        display as error "Could not load phillips dataset (no local cache, no network)."
        exit 601
    }
}
quietly tsset year
save "`outdir'/_gmm_phil_temp.dta", replace


/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV
  (Same structure as other fixture generators; reused verbatim -- the R tests parse these column names, so do not change this program.)
---------------------------------------------------------------------------*/
capture program drop save_ivreg2_results
program define save_ivreg2_results
    syntax, prefix(string) suffix(string) outdir(string) [diagonly]

    local N = e(N)
    local K = e(rankxx)
    local L = e(rankzz)

    // --- Coefficients and SEs ---
    matrix b = e(b)
    matrix V = e(V)
    local names : colnames b
    local ncols = colsof(b)

    if "`diagonly'" == "" {
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

    // --- Full VCV matrix (term + vcov_<name> headers, the shared expect_vcov_fixture / audit-reader format; replaced the legacy v1..vK term-less export at the M-16 review) ---
    matrix V = e(V)
    quietly {
        preserve
        clear
        local dim = rowsof(V)
        svmat double V
        gen str32 term = ""
        forvalues i = 1/`dim' {
            local nm : word `i' of `names'
            replace term = "`nm'" in `i'
        }
        order term
        forvalues i = 1/`dim' {
            local nm : word `i' of `names'
            local cnm = subinstr("`nm'", ".", "_", .)
            rename V`i' vcov_`cnm'
        }
        export delimited using "`outdir'/`prefix'_vcov_`suffix'.csv", replace
        restore
    }
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
  BLOCK 1: griliches76
  Base model = help-file H06:
    ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), gmm2s robust
===========================================================================*/

// --- D5a: gmm2s robust small -- GMM2S applies its small-sample corrections on a method-specific path, so small x gmm2s needs its own fixture ---
use "`outdir'/_gmm_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + gmm2s robust small ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), gmm2s robust small
save_ivreg2_results, prefix(gril_gmm2s) suffix(robust_small) outdir(`outdir')

// --- D5a: deterministic aweight; also serves the retired pweight cell via the byte-identity noted above ---
use "`outdir'/_gmm_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + [aw=awt] gmm2s robust ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt) [aw=awt], gmm2s robust
save_ivreg2_results, prefix(gril_gmm2s_aw) suffix(robust) outdir(`outdir')

// --- D5a: dofminus(2) ---
use "`outdir'/_gmm_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + gmm2s robust dofminus(2) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), gmm2s robust dofminus(2)
save_ivreg2_results, prefix(gril_gmm2s_dof) suffix(robust) outdir(`outdir')

// --- D5a: endog(iq) -- iq is the endogenous regressor. Diagnostics-only export: endog() does not change estimation, so the coef/vcov CSVs were byte-identical to hf_gril H06's (verified at the M-16 review); test-gmm2s.R pins that invariance fixture-free instead ---
use "`outdir'/_gmm_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + gmm2s robust endog(iq) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), gmm2s robust endog(iq)
save_ivreg2_results, prefix(gril_gmm2s_endog) suffix(robust) outdir(`outdir') diagonly


/*===========================================================================
  BLOCK 2: abdata
  Base model = help-file H88:
    ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), gmm2s cluster(id)
===========================================================================*/

// --- D5a: gmm2s cluster(id) small -- the non-small cluster cell is tsop_ab_gmm2s_cl ---
use "`outdir'/_gmm_ab_temp.dta", clear
display _newline(2) "=== abdata H88 + gmm2s cluster(id) small ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), gmm2s cluster(id) small
save_ivreg2_results, prefix(ab_gmm2s) suffix(cl_small) outdir(`outdir')

// --- D5a: weighted x cluster x gmm2s intersection, replacing the retired card [aw=weight] cluster(age) cell ---
use "`outdir'/_gmm_ab_temp.dta", clear
display _newline(2) "=== abdata H88 + [aw=abwt] gmm2s cluster(id) ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys) [aw=abwt], gmm2s cluster(id)
save_ivreg2_results, prefix(ab_gmm2s_aw) suffix(cl) outdir(`outdir')


/*===========================================================================
  BLOCK 3: phillips
  Base model family = help-file H83-H85: ivreg2 cinf (unem = l(1/3).unem), ...
===========================================================================*/

// --- D5a variation on H83/H84: kernel defaults to Bartlett, NO robust, so this exercises the AC (homoskedastic-kernel) weighting-matrix path -- regression guard for the GMM2S-AC weighting-matrix bug fixed in commit 0a405e8 ---
use "`outdir'/_gmm_phil_temp.dta", clear
display _newline(2) "=== phillips + gmm2s bw(3) (AC path, no robust) ==="
ivreg2 cinf (unem = l(1/3).unem), gmm2s bw(3)
save_ivreg2_results, prefix(phil_gmm2s) suffix(ac_bw3) outdir(`outdir')


/*===========================================================================
  Clean up temp files and done
===========================================================================*/
capture erase "`outdir'/_gmm_gril_temp.dta"
capture erase "`outdir'/_gmm_ab_temp.dta"
capture erase "`outdir'/_gmm_phil_temp.dta"
display _newline(2) "=== All GMM fixtures generated ==="
display "Output directory: `outdir'"
