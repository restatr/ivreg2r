/*===========================================================================
  generate-partial-fixtures.do
  -----------------------------
  Generates CSV benchmark fixtures for the partial()/FWL fixture family: 10 cells on the griliches76 partial arc, help.txt:1250-1262. The base model is H28 (help.txt:1253) minus its cluster(year) option: `ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), partial(_I*)`. Stata documents partial() with worked examples only in the cluster context (H27-H29), so every cell here is a D5a option-variation on that base.

  NOT duplicated here: cluster(year)+partial(_I*) Stata parity is owned by the helpfile suite (hf_gril H28, generate-helpfile-fixtures.do, with the M = 7 clusters vs L = 8 instruments rank-deficiency warnings asserted there as the documented lesson); the H29 gmm2s+cluster+partial command is a DOCUMENTED HELP-FILE BUG owned by the hf suite as an assert-only cell (Stata rejects it with r(506)). The gmm2s cell below uses robust instead of cluster -- with robust the partialled GMM2S model is feasible; H29's infeasibility is specific to #clusters < #moments.

  Deleted, not ported: card_partial_all (partial(_all) Stata parity is owned by generate-partial-all-fixtures.do, on mroz); card_partial_cluster cl/cl_small (cluster(smsa66) is the binary M=2 cluster anti-pattern; cluster x partial lives at hf H28); the remaining card cells are replaced 1:1 by the griliches cells below.

  Invariance retirement (byte-identity verified in the retired card fixtures, 2026-07-06): nopartialsmall under plain IID leaves e(b) and e(V) byte-identical to the basic cell -- it only zeroes the sdofminus increment, which moves df-dependent diagnostics (e(cdf), e(F), e(Fdf2), e(r2_a)) -- so the nosmall cell exports DIAGNOSTICS ONLY and the R test asserts the nosmall fit's coef/vcov directly equal to the basic fit. The nosmall x small cell is retired compositionally: under IID small the basic and nosmall fits share coefficients, rss, and (X'X)^-1, so their VCVs differ by exactly the ratio df_r(nosmall)/df_r(basic); R pins this as a fixture-free correction-factor identity (M-14 cl2_ols_small precedent), break-it-once verified.

  The aweight cell uses the deterministic weight awt = mod(age,5)+1 (M-12/M-23 precedent; the retired card cell's [aw=weight] survey weight is replaced). R mirrors it via the shared griliches_awt helper object.

  No data CSV is exported by this generator (R consumes bundled data(griliches); the %21.0g+datafmt data-export rule is therefore moot here).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-partial-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"


/*===========================================================================
  Load griliches76 FIRST, before the program is defined (data loads can interact with programs, so mirror the M-16/M-18 generators' ordering: every use/save block ahead of "program define").
===========================================================================*/

// --- griliches76: local cache, then BC.edu web fallback (per the H28 base model) ---
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
save "`outdir'/_partial_gril_temp.dta", replace


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

        // Weak identification. Deliberate deviation from the center-generator copy: weak_id_kp_f is taken from e(rkf), not e(widstat), because this family has iid cells -- under iid, e(widstat) falls back to the Cragg-Donald F (which R stores in weak_id, already pinned by weak_id_cd_f), while e(rkf) is posted only for robust-class VCE, exactly matching when R populates weak_id_robust. The first_stage_f = e(rkf) column of the center copy is dropped here as it would duplicate weak_id_kp_f.
        gen double weak_id_cd_f = e(cdf)
        gen double weak_id_kp_f = e(rkf)

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

        // Partial bookkeeping (this family's raison d'etre).
        // These four columns are the M-24 additions to the shared schema; test-partial.R's bookkeeping assertions consume them. Stata posts e(df_r) only under the small option, so df_r is blank in the non-small cells and the R assertion is conditional on it.
        gen int df_r = e(df_r)
        gen int sdofminus = e(sdofminus)
        gen int partial_ct = e(partial_ct)
        gen int partialcons = e(partialcons)

        export delimited using "`outdir'/`prefix'_diagnostics_`suffix'.csv", replace
        restore
    }
end


/*===========================================================================
  griliches76 partial() arc
  Base model = help-file H28 (help.txt:1253) minus cluster(year):
    ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), partial(_I*)
  All cells here are D5a option-variations (Stata documents partial() with
  worked examples only in the cluster context, H27-H29).
===========================================================================*/

// --- 1. D5a base: iid ---
use "`outdir'/_partial_gril_temp.dta", clear
display _newline(2) "=== griliches H28 (minus cluster) + partial(_I*), iid ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), partial(_I*)
save_ivreg2_results, prefix(gril_partial) suffix(iid) outdir(`outdir')

// --- 2. D5a base: iid small ---
use "`outdir'/_partial_gril_temp.dta", clear
display _newline(2) "=== griliches H28 (minus cluster) + partial(_I*), iid small ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), partial(_I*) small
save_ivreg2_results, prefix(gril_partial) suffix(iid_small) outdir(`outdir')

// --- 3. D5a: robust ---
use "`outdir'/_partial_gril_temp.dta", clear
display _newline(2) "=== griliches H28 (minus cluster) + partial(_I*), robust ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), partial(_I*) robust
save_ivreg2_results, prefix(gril_partial) suffix(robust) outdir(`outdir')

// --- 4. D5a: robust small ---
use "`outdir'/_partial_gril_temp.dta", clear
display _newline(2) "=== griliches H28 (minus cluster) + partial(_I*), robust small ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), partial(_I*) robust small
save_ivreg2_results, prefix(gril_partial) suffix(robust_small) outdir(`outdir')

// --- 5. D5a: aweight, deterministic weight awt = mod(age,5)+1 (M-12/M-23 precedent) ---
use "`outdir'/_partial_gril_temp.dta", clear
display _newline(2) "=== griliches H28 (minus cluster) + partial(_I*), aweight ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age) [aw=awt], partial(_I*)
save_ivreg2_results, prefix(gril_partial_aw) suffix(iid) outdir(`outdir')

// --- 6. D5a: LIML ---
use "`outdir'/_partial_gril_temp.dta", clear
display _newline(2) "=== griliches H28 (minus cluster) + partial(_I*), LIML ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), partial(_I*) liml
save_ivreg2_results, prefix(gril_partial_liml) suffix(iid) outdir(`outdir')

// --- 7. D5a: GMM2S robust -- feasible with robust (contrast H29's cluster infeasibility, r(506), owned by the hf suite) ---
use "`outdir'/_partial_gril_temp.dta", clear
display _newline(2) "=== griliches H28 (minus cluster) + partial(_I*), GMM2S robust ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), partial(_I*) gmm2s robust
save_ivreg2_results, prefix(gril_partial_gmm2s) suffix(robust) outdir(`outdir')

// --- 8. D5a: partial(_cons) -- demean-only path, partial_ct = 1, partialcons = 1 ---
use "`outdir'/_partial_gril_temp.dta", clear
display _newline(2) "=== griliches H28 (minus cluster) + partial(_cons) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), partial(_cons)
save_ivreg2_results, prefix(gril_partial_cons) suffix(iid) outdir(`outdir')

// --- 9. D5a: OLS -- no endogenous regressor, pins the OLS x partial fit path ---
use "`outdir'/_partial_gril_temp.dta", clear
display _newline(2) "=== griliches H28 (minus cluster), OLS + partial(_I*) ==="
ivreg2 lw s expr tenure rns smsa _I*, partial(_I*)
save_ivreg2_results, prefix(gril_partial_ols) suffix(iid) outdir(`outdir')

// --- 10. D5a: nopartialsmall, DIAGONLY -- coef/vcov byte-identical to cell 1 (invariance retirement, see header) ---
use "`outdir'/_partial_gril_temp.dta", clear
display _newline(2) "=== griliches H28 (minus cluster) + partial(_I*) nopartialsmall ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), partial(_I*) nopartialsmall
save_ivreg2_results, prefix(gril_partial_nosmall) suffix(iid) outdir(`outdir') diagonly


/*===========================================================================
  Clean up temp file and done
===========================================================================*/
capture erase "`outdir'/_partial_gril_temp.dta"
display _newline(2) "=== All partial fixtures generated ==="
display "Output directory: `outdir'"
