/*===========================================================================
  generate-center-fixtures.do
  ---------------------------
  Generates CSV benchmark fixtures for the center-option family (M-18 re-base, planning/22-spec-matrix.md): the old 13 card cells are replaced by 10 cells layered as D5a option-variations on the M-16 canonical bases -- griliches76 H06, abdata H88, and the phillips H83-H85 family. Stata documents the center option with no worked example, so every cell here is a D5a option-variation on a help-file base model.

  Both CUE+center cells are DELETED, not ported: the CUE+cluster+center cell is the known-dirty ubuntu CI basin-pathology cell per the spec-matrix delete table, and the CUE robust+center cell is the same card-CUE known-dirty class -- CUE x center coverage moves to M-17's canonical bases.

  The just-identified center cell is DELETED: center enters only through the shared meat helpers, with no identification-specific branch (M-12 precedent), so a just-identified center cell adds no distinct code path.

  endog() does not change estimation -- the retired card endog cell's coef/vcov CSVs were byte-identical to the plain hc0 cell's (verified 2026-07-05) -- and Stata does NOT forward center to the recursive endog call (ivreg2.ado ~1576-1601 omit it; documented intentional-match row in planning/07c line 198), so cell 5 exports diagnostics only and its e(estat) pins the uncentered-inner-call behavior both implementations share.

  The retired card orthog+center cell was vacuous -- all three of its CSVs were byte-identical to the hc0 cell's because the old export never recorded e(cstat) -- so cell 6 pins the orthog C-statistic under center for the FIRST time; Stata's recursive orthog call receives the centered full-model S via smatrix() (ivreg2.ado ~1522-1536), so cstat IS center-affected.

  Generation-time small-invariance assertions (cells 5 and 6): this file re-runs cell 5 and cell 6 with `small` added and asserts e(estat)/e(estatp) and e(cstat)/e(cstatp) reproduce to reldif < 1e-12 with exact df equality. This is the evidence the R tests' double-small-fit harness relies on. M-22 already verified cstat small-invariance byte-identical in its retired fixtures; the estat analog has no retired byte-identity pair in this family, so this assert is its evidence.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-center-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"


/*===========================================================================
  Load all three datasets FIRST, before the program is defined (data loads can interact with programs, so mirror the M-16 generator's ordering: every use/save block ahead of "program define").
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
save "`outdir'/_center_gril_temp.dta", replace

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
save "`outdir'/_center_ab_temp.dta", replace

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
save "`outdir'/_center_phil_temp.dta", replace


/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV
  (Copied verbatim from generate-gmm-fixtures.do -- the R tests parse these column names, including the term + vcov_<name> vcov headers, so do not change this program.)
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


/*---------------------------------------------------------------------------
  Helper program: save orthog-specific diagnostics to CSV.
  (Copied from generate-orthog-fixtures.do -- exports cstat/cstatp/cstatdf plus j/jdf/jp so read_diagnostics() consumers work identically to the M-22 fixtures.)
---------------------------------------------------------------------------*/
capture program drop save_orthog_results
program define save_orthog_results
    syntax, prefix(string) suffix(string) outdir(string)

    quietly {
        preserve
        clear
        set obs 1

        gen double cstat = .
        gen double cstatp = .
        gen double cstatdf = .
        gen str80 clist = ""
        gen double j = .
        gen double jp = .
        gen double jdf = .

        // Orthogonality C-statistic
        capture replace cstat = e(cstat)
        capture replace cstatp = e(cstatp)
        capture replace cstatdf = e(cstatdf)
        capture replace clist = "`e(clist)'"

        // Overidentification (for reference / sanity checks)
        capture replace j = e(j)
        capture replace jp = e(jp)
        capture replace jdf = e(jdf)
        // IID path uses sargan instead of j
        capture {
            if missing(j) {
                replace j = e(sargan)
                replace jp = e(sarganp)
                replace jdf = e(sargandf)
            }
        }

        export delimited using ///
            "`outdir'/`prefix'_orthog_`suffix'.csv", replace
        restore
    }
end


/*===========================================================================
  BLOCK 1: griliches76
  Base model = help-file H06 (help.txt:1154) minus gmm2s except where noted:
    ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), ...
  All center cells here are D5a option-variations (Stata documents center with no worked example).
===========================================================================*/

// --- 1. D5a: robust center -- replaces the retired card hc0 cell; 2SLS + centered HC meat ---
use "`outdir'/_center_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + robust center ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), robust center
save_ivreg2_results, prefix(gril_center) suffix(robust) outdir(`outdir')

// --- 2. D5a: robust small center -- small x center on the plain-2SLS path (small correction composes with the centered S) ---
use "`outdir'/_center_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + robust small center ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), robust small center
save_ivreg2_results, prefix(gril_center) suffix(robust_small) outdir(`outdir')

// --- 3. D5a: gmm2s robust center -- center changes the GMM2S weighting matrix and hence the coefficients; this is the gmm2s+center replacement direction named in the spec-matrix M-18 row ---
use "`outdir'/_center_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + gmm2s robust center ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), gmm2s robust center
save_ivreg2_results, prefix(gril_gmm2s_center) suffix(robust) outdir(`outdir')

// --- 4. D5a: robust center dofminus(2) -- dofminus(2) matches the M-16 dof-cell precedent ---
use "`outdir'/_center_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + robust center dofminus(2) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), robust center dofminus(2)
save_ivreg2_results, prefix(gril_center_dof) suffix(robust) outdir(`outdir')

// --- 5. D5a: robust center endog(iq) -- diagnostics-only export. endog() does not change estimation (retired card endog cell's coef/vcov CSVs were byte-identical to the plain hc0 cell's, verified 2026-07-05), and Stata does NOT forward center to the recursive endog call (ivreg2.ado ~1576-1601 omit it; intentional-match row in planning/07c line 198), so e(estat) here pins the uncentered-inner-call behavior both implementations share ---
use "`outdir'/_center_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + robust center endog(iq) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), robust center endog(iq)
scalar estat_base = e(estat)
scalar estatp_base = e(estatp)
scalar estatdf_base = e(estatdf)
save_ivreg2_results, prefix(gril_center_endog) suffix(robust) outdir(`outdir') diagonly

// --- 6. D5a: robust center orthog(age mrt) -- orthog C-stat export. The retired card orthog cell was vacuous (all three CSVs byte-identical to the hc0 cell's because the old export never recorded e(cstat)), so this cell pins the orthog C-stat under center for the FIRST time; Stata's recursive orthog call receives the centered full-model S via smatrix() (ivreg2.ado ~1522-1536), so cstat IS center-affected ---
use "`outdir'/_center_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + robust center orthog(age mrt) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), robust center orthog(age mrt)
scalar cstat_base = e(cstat)
scalar cstatp_base = e(cstatp)
scalar cstatdf_base = e(cstatdf)
save_orthog_results, prefix(gril_center) suffix(robust) outdir(`outdir')


/*===========================================================================
  Generation-time small-invariance assertions (cells 5 and 6).
  This is the evidence the R tests' double-small-fit harness relies on: it certifies that adding `small` leaves the endogeneity C-statistic (cell 5) and the orthogonality C-statistic (cell 6) unchanged, so a single fixture serves both small = FALSE and small = TRUE R fits. On assert failure the batch run errors and the main loop treats it as a gate mismatch.
===========================================================================*/

// --- Cell 5 small-invariance: re-run robust center endog(iq) WITH small and certify e(estat)/e(estatp) reldif < 1e-12, e(estatdf) exact ---
use "`outdir'/_center_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + robust center endog(iq) small (invariance check) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), robust center endog(iq) small
scalar estat_small = e(estat)
scalar estatp_small = e(estatp)
scalar estatdf_small = e(estatdf)
// certifies the endogeneity statistic is invariant to small
assert reldif(estat_base, estat_small) < 1e-12
// certifies the endogeneity p-value is invariant to small
assert reldif(estatp_base, estatp_small) < 1e-12
// certifies the endogeneity df is exactly equal under small
assert estatdf_base == estatdf_small

// --- Cell 6 small-invariance: re-run robust center orthog(age mrt) WITH small and certify e(cstat)/e(cstatp) reldif < 1e-12, e(cstatdf) exact ---
use "`outdir'/_center_gril_temp.dta", clear
display _newline(2) "=== griliches H06 + robust center orthog(age mrt) small (invariance check) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), robust center orthog(age mrt) small
scalar cstat_small = e(cstat)
scalar cstatp_small = e(cstatp)
scalar cstatdf_small = e(cstatdf)
// certifies the orthogonality C-statistic is invariant to small
assert reldif(cstat_base, cstat_small) < 1e-12
// certifies the orthogonality C-stat p-value is invariant to small
assert reldif(cstatp_base, cstatp_small) < 1e-12
// certifies the orthogonality C-stat df is exactly equal under small
assert cstatdf_base == cstatdf_small


/*===========================================================================
  BLOCK 2: abdata
  Base model = help-file H88 (help.txt:1541), tsset id year already applied:
    ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), ...
===========================================================================*/

// --- 7. D5a: cluster(id) center -- replaces the retired card cluster(smsa) M=2 anti-pattern cell; id gives M=140 clusters ---
use "`outdir'/_center_ab_temp.dta", clear
display _newline(2) "=== abdata H88 + cluster(id) center ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), cluster(id) center
save_ivreg2_results, prefix(ab_center) suffix(cl) outdir(`outdir')

// --- 8. D5a: cluster(id) small center ---
use "`outdir'/_center_ab_temp.dta", clear
display _newline(2) "=== abdata H88 + cluster(id) small center ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), cluster(id) small center
save_ivreg2_results, prefix(ab_center) suffix(cl_small) outdir(`outdir')

// --- 9. D5a: gmm2s cluster(id) center -- the designated replacement for the deleted CUE+cluster+center headline cell ("drop ... in favor of gmm2s+center", spec-matrix M-18 row + delete table) ---
use "`outdir'/_center_ab_temp.dta", clear
display _newline(2) "=== abdata H88 + gmm2s cluster(id) center ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), gmm2s cluster(id) center
save_ivreg2_results, prefix(ab_gmm2s_center) suffix(cl) outdir(`outdir')


/*===========================================================================
  BLOCK 3: phillips
  Base model family = help-file H83-H85: ivreg2 cinf (unem = l(1/3).unem), ...
  tsset year already applied.
===========================================================================*/

// --- 10. D5a: robust bw(3) kernel(bartlett) center -- replaces the retired card pseudo-time-series HAC cell (gen t=_n; tsset t on cross-sectional card, an anti-pattern) with a genuine time series ---
use "`outdir'/_center_phil_temp.dta", clear
display _newline(2) "=== phillips + robust bw(3) kernel(bartlett) center ==="
ivreg2 cinf (unem = l(1/3).unem), robust bw(3) kernel(bartlett) center
save_ivreg2_results, prefix(phil_center) suffix(hac_bw3) outdir(`outdir')


/*===========================================================================
  Clean up temp files and done
===========================================================================*/
capture erase "`outdir'/_center_gril_temp.dta"
capture erase "`outdir'/_center_ab_temp.dta"
capture erase "`outdir'/_center_phil_temp.dta"
display _newline(2) "=== All center fixtures generated ==="
display "Output directory: `outdir'"
