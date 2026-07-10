/*===========================================================================
  generate-cue-fixtures.do
  ------------------------
  Generates CSV benchmark fixtures for the CUE family: 10 cells layered as D5a option-variations on the canonical bases -- the stockwatson BSS 2007 pp. 476/480 CUE arc and the abdata H88 model minus gmm2s. Stata's help file has no worked CUE examples beyond H22/H74, so the BSS Stata Journal paper (Baum, Schaffer & Stillman 2007, 7(4):465-506) supplies the canonical worked CUE arc that every new stockwatson cell varies.

  Two upstream fixtures are REUSED, not regenerated here: klein H74 CUE-IID parity (CUE = LIML on a just-identified model) is owned by the tsop_klein_cue fixtures in generate-ts-operator-fixtures.do, and mroz H68 b0 parity (CUE evaluated at a fixed b0 vector) is owned by the hf_mroz H68 fixtures in generate-helpfile-fixtures.do. Neither is touched by this file.

  Old cells DELETED, not ported (omitted entirely, not commented out): all card CUE cells (known-dirty optimizer-basin class; the cluster(age) cells were likewise deleted); the card just-identified CUE cell (retired to the fixture-free just-id CUE = 2SLS identity); the [aw=age] weighted cell (non-deterministic weight anti-pattern, replaced below by the deterministic swwt weight); the card b0 cells and card_overid_b0_vector.csv (b0 Stata parity is owned by the hf_mroz H68 fixture); the synthetic seed-12345 DGP with its ts_cue cell and the ts_gmm_data.csv data export.

  The retired iid_small cell is retired COMPOSITIONALLY: the small correction applies in the shared VCV assembly regardless of omega type, so a dedicated CUE-IID small cell adds no distinct code path; the sw_cue robust_small cell (cell 3) pins the small correction for CUE.

  stockwatson construction-coherence rationale: this file constructs the derived columns EXACTLY as pkg/data-raw/stockwatson.R builds them, in double precision, so the Stata data is numerically coherent with the bundled data(stockwatson) the R tests fit on. Both compute in double from the same float32 macrodat source, following the F4 precision-lesson discipline (float32 truncation of a CSV round-trip previously caused a false CUE non-convergence). No data CSV is exported by this file; the R side loads data(stockwatson) directly.

  Generation-time small-invariance assertions (cells 7 and 8): this file re-runs cell 7 (endog) and cell 8 (orthog) with `small` added and asserts e(estat)/e(estatp) and e(cstat)/e(cstatp) reproduce to reldif < 1e-12 with exact df equality. This is the evidence the R tests' double-small-fit harness relies on: a single fixture serves both small = FALSE and small = TRUE R fits. On assert failure the batch run errors and the main loop treats it as a gate finding.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-cue-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"


/*===========================================================================
  Load both datasets FIRST, before the program is defined (data loads can interact with programs, so mirror the M-16/M-18 generator ordering: every use/save block ahead of "program define").
===========================================================================*/

// --- stockwatson: local cache, then BC.edu web fallback (per the BSS p. 480 base model) ---
capture use "../validation/data/macrodat.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/stockwatson/macrodat.dta, clear
    if _rc {
        display as error "Could not load stockwatson (macrodat) dataset (no local cache, no network)."
        exit 601
    }
}
// Construct the derived columns EXACTLY as pkg/data-raw/stockwatson.R does, in double precision, so the Stata data is numerically coherent with the bundled data(stockwatson) the R tests fit on (F4 precision-lesson discipline: both sides compute in double from the same float32 macrodat source).
sort date
tsset date
gen double inf = 100 * ln(CPI / L4.CPI)
gen double ggdp = 100 * ln(GDP / L4.GDP)
gen double dinf = inf - L.inf
gen double ggdp_2 = L2.ggdp
gen double TBILL_1 = L.TBILL
gen double ER_1 = L.ER
gen double TBON_1 = L.TBON
// swwt = mod(date,4)+1 is the deterministic analytic weight (M-12 precedent): date is the Stata quarterly index, so this is identical to R's date %% 4 + 1, and it replaces the retired non-deterministic [aw=age] anti-pattern cell.
gen double swwt = mod(date, 4) + 1
save "`outdir'/_cue_sw_temp.dta", replace

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
save "`outdir'/_cue_ab_temp.dta", replace


/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV
  (Copied verbatim from generate-center-fixtures.do -- the R tests parse these column names, including the term + vcov_<name> vcov headers, so do not change this program.)
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
  (Copied verbatim from generate-center-fixtures.do -- exports cstat/cstatp/cstatdf plus j/jdf/jp so read_diagnostics() consumers work identically to the M-22 fixtures.)
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
  BLOCK 1: stockwatson
  Base model = BSS 2007 p. 480 (Section 4.3), the canonical worked CUE arc:
    ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue ...
  All cells here are D5a option-variations on that base (Stata's help file documents CUE with no worked example beyond H22/H74).
===========================================================================*/

// --- 1. D5a: cue robust kernel(bartlett) bw(5) -- THE canonical cell (BSS 2007 Section 4.3, p. 480; same command as the stockwatson stored-matrices cell in generate-stored-matrices-fixtures.do:141). Replaces the retired synthetic seed-12345 ts_cue HAC cell with a genuine BSS-documented HAC-CUE run ---
use "`outdir'/_cue_sw_temp.dta", clear
display _newline(2) "=== stockwatson BSS p.480 + cue robust kernel(bartlett) bw(5) ==="
ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue robust kernel(bartlett) bw(5)
save_ivreg2_results, prefix(sw_cue) suffix(hac_bw5) outdir(`outdir')

// --- 2. D5a: cue robust -- HC0 CUE on the BSS base; replaces the retired card CUE robust cell ---
use "`outdir'/_cue_sw_temp.dta", clear
display _newline(2) "=== stockwatson BSS p.480 + cue robust ==="
ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue robust
save_ivreg2_results, prefix(sw_cue) suffix(robust) outdir(`outdir')

// --- 3. D5a: cue robust small -- the small correction applies in the shared VCV assembly regardless of omega type, so this cell pins the small correction for CUE and compositionally retires the old iid_small cell ---
use "`outdir'/_cue_sw_temp.dta", clear
display _newline(2) "=== stockwatson BSS p.480 + cue robust small ==="
ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue robust small
save_ivreg2_results, prefix(sw_cue) suffix(robust_small) outdir(`outdir')

// --- 4. D5a: cue robust center -- the CUE x center variation deferred from M-18 (M-18 deleted both card CUE x center cells without porting; CUE x center coverage was explicitly moved to M-17's canonical base) ---
use "`outdir'/_cue_sw_temp.dta", clear
display _newline(2) "=== stockwatson BSS p.480 + cue robust center ==="
ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue robust center
save_ivreg2_results, prefix(sw_cue_center) suffix(robust) outdir(`outdir')

// --- 5. D5a: cue robust dofminus(2) -- dofminus(2) matches the M-16/M-18 dof-cell precedent ---
use "`outdir'/_cue_sw_temp.dta", clear
display _newline(2) "=== stockwatson BSS p.480 + cue robust dofminus(2) ==="
ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue robust dofminus(2)
save_ivreg2_results, prefix(sw_cue_dof) suffix(robust) outdir(`outdir')

// --- 6. D5a/D6: [aw=swwt] cue robust -- deterministic analytic weight (swwt = mod(date,4)+1), replacing the retired non-deterministic [aw=age] anti-pattern cell ---
use "`outdir'/_cue_sw_temp.dta", clear
display _newline(2) "=== stockwatson BSS p.480 + [aw=swwt] cue robust ==="
ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1) [aw=swwt], cue robust
save_ivreg2_results, prefix(sw_cue_aw) suffix(robust) outdir(`outdir')

// --- 7. D5a: cue robust endog(UR) -- diagnostics-only export. endog() does not change estimation (M-16/M-18 byte-identity precedent), so e(estat) here pins the endogeneity C-statistic both implementations share; capture the estat scalars for the small-invariance assert below ---
use "`outdir'/_cue_sw_temp.dta", clear
display _newline(2) "=== stockwatson BSS p.480 + cue robust endog(UR) ==="
ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue robust endog(UR)
scalar estat_base = e(estat)
scalar estatp_base = e(estatp)
scalar estatdf_base = e(estatdf)
save_ivreg2_results, prefix(sw_cue_endog) suffix(robust) outdir(`outdir') diagonly

// --- 8. D5a: cue robust orthog(TBON_1) -- orthog C-stat export. The retired card CUE orthog cell was read by the tests but asserted nothing from it (the old export never recorded e(cstat)), so this cell pins the CUE orthog C-statistic (computed against the CUE-defining omega, the F5 fix path) for the FIRST time; capture the cstat scalars for the small-invariance assert below ---
use "`outdir'/_cue_sw_temp.dta", clear
display _newline(2) "=== stockwatson BSS p.480 + cue robust orthog(TBON_1) ==="
ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue robust orthog(TBON_1)
scalar cstat_base = e(cstat)
scalar cstatp_base = e(cstatp)
scalar cstatdf_base = e(cstatdf)
save_orthog_results, prefix(sw_cue) suffix(robust) outdir(`outdir')


/*===========================================================================
  Generation-time small-invariance assertions (cells 7 and 8) — see the file header for the full rationale.
===========================================================================*/

// --- Cell 7 small-invariance: re-run cue robust endog(UR) WITH small and certify e(estat)/e(estatp) reldif < 1e-12, e(estatdf) exact ---
use "`outdir'/_cue_sw_temp.dta", clear
display _newline(2) "=== stockwatson BSS p.480 + cue robust endog(UR) small (invariance check) ==="
ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue robust endog(UR) small
scalar estat_small = e(estat)
scalar estatp_small = e(estatp)
scalar estatdf_small = e(estatdf)
// certifies the endogeneity statistic is invariant to small
assert reldif(estat_base, estat_small) < 1e-12
// certifies the endogeneity p-value is invariant to small
assert reldif(estatp_base, estatp_small) < 1e-12
// certifies the endogeneity df is exactly equal under small
assert estatdf_base == estatdf_small

// --- Cell 8 small-invariance: re-run cue robust orthog(TBON_1) WITH small and certify e(cstat)/e(cstatp) reldif < 1e-12, e(cstatdf) exact ---
use "`outdir'/_cue_sw_temp.dta", clear
display _newline(2) "=== stockwatson BSS p.480 + cue robust orthog(TBON_1) small (invariance check) ==="
ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue robust orthog(TBON_1) small
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
  Base model = help-file H88 (help.txt:1541) minus gmm2s, tsset id year already applied:
    ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), cue ...
  These cells carry the CUE x cluster and CUE x (K1=3 multi-endogenous) intersections, replacing the retired card cluster(age) cells.
===========================================================================*/

// --- 9. D5a: cue cluster(id) -- CUE x cluster and CUE x multi-endogenous intersections; id gives M=140 clusters ---
use "`outdir'/_cue_ab_temp.dta", clear
display _newline(2) "=== abdata H88 + cue cluster(id) ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), cue cluster(id)
save_ivreg2_results, prefix(ab_cue) suffix(cl) outdir(`outdir')

// --- 10. D5a: cue cluster(id) small ---
use "`outdir'/_cue_ab_temp.dta", clear
display _newline(2) "=== abdata H88 + cue cluster(id) small ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), cue cluster(id) small
save_ivreg2_results, prefix(ab_cue) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  Clean up temp files and done
===========================================================================*/
capture erase "`outdir'/_cue_sw_temp.dta"
capture erase "`outdir'/_cue_ab_temp.dta"
display _newline(2) "=== All CUE fixtures generated ==="
display "Output directory: `outdir'"
