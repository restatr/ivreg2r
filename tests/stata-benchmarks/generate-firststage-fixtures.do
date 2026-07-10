/*===========================================================================
  generate-firststage-fixtures.do
  --------------------------------
  Generates CSV benchmark fixtures for first-stage diagnostics and first-stage model objects on the canonical help-file bases: mroz (H31/H41, help.txt:1274/1325-1359), griliches76 (H03/H04, help.txt:1138/1141), and abdata (H88 minus gmm2s, help.txt:1541).

  `small` variants are deliberately NOT generated for any cell in this file: Stata's first-stage diagnostics (e(first)) were verified byte-identical under `small` across iid/hc1/cluster in the retired card/sim fixtures, so the R tests fit both small = FALSE and small = TRUE against each fixture here and assert the whole first_stage object is equal either way.
  No liml/kclass/gmm2s cells are generated: the saved first-stage regression is the same OLS fit regardless of the main-equation estimator -- fs_overid_liml_* was verified byte-identical to fs_overid_iid_* in the retired card fixtures, so the R tests assert this estimator-invariance as a fixture-free identity rather than pinning a duplicate fixture.

  Two kinds of fixtures are produced:
    - DIAGNOSTICS (13 cells, via save_firststage): the e(first) summary matrix (rmse, sheapr2, pr2, F, df, df_r, pvalue, Sanderson-Windmeijer and Angrist-Pischke conditional-F/chi2/partial-R2 rows), one row per statistic and one column per endogenous regressor.
    - OBJECT (5 cells, via save_fs_fixture): the saved first-stage regression equation itself (coefficients, VCV, RMSE, df_r, N, N_clust) for a single named endogenous variable, restored from e(firsteqs) via `estimates restore`.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-firststage-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"


/*===========================================================================
  Load all three datasets FIRST, before the program is defined (bcuse-style loaders can issue "clear all" internally; keeping every use/save block ahead of "program define" avoids any ordering hazard).
===========================================================================*/

// --- mroz: local cache, then BC.edu web fallback (per the H31 base model, help.txt:1274) ---
capture use "../validation/data/mroz.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/wooldridge/mroz.dta, clear
    if _rc {
        display as error "Could not load mroz dataset (no local cache, no network)."
        exit 601
    }
}
save "`outdir'/_fs_mroz_temp.dta", replace

// --- griliches76: local cache, then BC.edu web fallback (per the H03/H04 base model, help.txt:1138/1141) ---
capture use "../validation/data/griliches76.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/hayashi/griliches76.dta, clear
    if _rc {
        display as error "Could not load griliches76 dataset (no local cache, no network)."
        exit 601
    }
}
quietly xi i.year
// Deterministic synthetic analytic weight (M-12/M-23 precedent): mod(age,5)+1 gives a reproducible, non-degenerate aweight without inventing an economically motivated weight.
gen awt = mod(age,5)+1
save "`outdir'/_fs_gril_temp.dta", replace

// --- abdata: local cache, then BC.edu web fallback (per the H88 base model minus gmm2s, help.txt:1541) ---
capture use "../validation/data/abdata.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/macro/abdata.dta, clear
    if _rc {
        display as error "Could not load abdata dataset (no local cache, no network)."
        exit 601
    }
}
// Deterministic synthetic analytic weight (M-12/M-23 precedent): mod(year,3)+1 gives a reproducible, non-degenerate aweight without inventing an economically motivated weight.
gen abwt = mod(year,3)+1
save "`outdir'/_fs_ab_temp.dta", replace


/*---------------------------------------------------------------------------
  Helper program: save first-stage summary diagnostics (e(first)) to CSV
  Extracted from the "First-stage diagnostics" block of save_ivreg2_results
  in generate-fixtures.do, parameterized with outdir() instead of the
  hardcoded "fixtures/" literal used there.
---------------------------------------------------------------------------*/
capture program drop save_firststage
program define save_firststage
    syntax, prefix(string) suffix(string) outdir(string)

    quietly {
        preserve
        clear

        // e(first) matrix: rows are statistics, columns are endogenous vars
        matrix F = e(first)
        local frows = rowsof(F)
        local fcols = colsof(F)
        local fnames : colnames F
        local rnames : rownames F

        set obs `frows'
        gen str32 statistic = ""
        forvalues i = 1/`frows' {
            local rn : word `i' of `rnames'
            replace statistic = "`rn'" in `i'
        }

        forvalues j = 1/`fcols' {
            local cn : word `j' of `fnames'
            local cnclean = subinstr("`cn'", ".", "_", .)
            gen double `cnclean' = .
            forvalues i = 1/`frows' {
                replace `cnclean' = F[`i', `j'] in `i'
            }
        }

        export delimited using ///
            "`outdir'/`prefix'_firststage_`suffix'.csv", replace
        restore
    }
end


/*---------------------------------------------------------------------------
  Helper program: save a stored first-stage equation (coefs + VCV + scalars)
  Reused from the prior version of this file, with outdir parameterized
  as an outdir() option instead of the hardcoded "fixtures/" literal.
---------------------------------------------------------------------------*/
capture program drop save_fs_fixture
program define save_fs_fixture
    syntax, prefix(string) eqname(string) outdir(string) [NOCOEF]

    estimates restore `eqname'

    * --- Coefficients (skipped under nocoef: the saved first-stage OLS coefficients are VCE-invariant, so robust cells would duplicate their iid counterpart byte-for-byte -- verified at re-base; the R tests compare the robust fit's coefficients against the iid coef fixture instead) ---
    matrix b = e(b)
    local K = colsof(b)
    local names : colnames b
    if "`nocoef'" == "" {
        file open fh using "`outdir'/`prefix'_coef.csv", write replace
        file write fh "term,estimate" _n
        forval i = 1/`K' {
            local nm : word `i' of `names'
            file write fh "`nm'," %20.12f (b[1,`i']) _n
        }
        file close fh
    }

    * --- VCV matrix ---
    matrix V = e(V)
    file open fh using "`outdir'/`prefix'_vcov.csv", write replace
    * Header row
    file write fh "term"
    forval j = 1/`K' {
        local nm : word `j' of `names'
        file write fh ",vcov_`nm'"
    }
    file write fh _n
    * Data rows
    forval i = 1/`K' {
        local rn : word `i' of `names'
        file write fh "`rn'"
        forval j = 1/`K' {
            file write fh "," %20.12f (V[`i',`j'])
        }
        file write fh _n
    }
    file close fh

    * --- Scalars ---
    file open fh using "`outdir'/`prefix'_scalars.csv", write replace
    file write fh "quantity,value" _n
    file write fh "rmse," %20.12f (e(rmse)) _n
    file write fh "df_r," %10.0f (e(df_r)) _n
    file write fh "N," %10.0f (e(N)) _n
    if !missing(e(N_clust)) {
        file write fh "N_clust," %10.0f (e(N_clust)) _n
    }
    file close fh

end


/*===========================================================================
  BLOCK 1: mroz
  Base model = H31 (help.txt:1274): ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6)
===========================================================================*/
use "`outdir'/_fs_mroz_temp.dta", clear

// --- Cell 1 (DIAGNOSTICS) + OBJECT fs_mroz_iid: H31 base, IID (help.txt:1274) ---
display _newline(2) "=== mroz first-stage, IID (H31) ==="
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6), ffirst savefirst
// save_firststage must run BEFORE save_fs_fixture: the "estimates restore"
// inside save_fs_fixture changes the active estimation results in memory.
save_firststage, prefix(mroz_fs) suffix(iid) outdir(`outdir')
save_fs_fixture, prefix(fs_mroz_iid) eqname(_ivreg2_educ) outdir(`outdir')

// --- Cell 2 (DIAGNOSTICS) + OBJECT fs_mroz_robust: H31 base, robust (H41/H46 arc, help.txt:1325) ---
display _newline(2) "=== mroz first-stage, robust (H41/H46) ==="
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6), robust ffirst savefirst
save_firststage, prefix(mroz_fs) suffix(hc1) outdir(`outdir')
save_fs_fixture, prefix(fs_mroz_robust) eqname(_ivreg2_educ) outdir(`outdir') nocoef

// --- Cell 3: D5a WITH DISCLOSURE -- H31 restricted to a single excluded instrument to exercise the just-identified (L1=1) first-stage path; instrument choice is arbitrary, not economically motivated ---
display _newline(2) "=== mroz first-stage, just-identified, IID (D5a) ==="
ivreg2 lwage exper expersq (educ=age), ffirst
save_firststage, prefix(mroz_fs_justid) suffix(iid) outdir(`outdir')

// --- Cell 4 (DIAGNOSTICS) + OBJECT fs_mroz_justid_robust: just-identified, robust (D5a) -- the object files restore the saved-equation coverage of the retired just-identified fs_just_robust config on the canonical base (review finding, 2026-07-04) ---
display _newline(2) "=== mroz first-stage, just-identified, robust (D5a) ==="
ivreg2 lwage exper expersq (educ=age), robust ffirst savefirst
save_firststage, prefix(mroz_fs_justid) suffix(hc1) outdir(`outdir')
save_fs_fixture, prefix(fs_mroz_justid_robust) eqname(_ivreg2_educ) outdir(`outdir')

// --- Cell 5: D5a -- the M-04 K2=0 noconstant spec -- no-intercept first stage (replaces retired sim_no_constant cells) ---
display _newline(2) "=== mroz first-stage, noconstant, IID (D5a) ==="
ivreg2 lwage (educ=age kidslt6 kidsge6), noconstant ffirst
save_firststage, prefix(mroz_fs_nocons) suffix(iid) outdir(`outdir')

// --- Cell 6: noconstant, robust (D5a) ---
display _newline(2) "=== mroz first-stage, noconstant, robust (D5a) ==="
ivreg2 lwage (educ=age kidslt6 kidsge6), noconstant robust ffirst
save_firststage, prefix(mroz_fs_nocons) suffix(hc1) outdir(`outdir')


/*===========================================================================
  BLOCK 2: griliches76
  Base model = H03 (help.txt:1138): ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt)
===========================================================================*/
use "`outdir'/_fs_gril_temp.dta", clear

// --- Cell 7: H03 base, IID, ffirst -- H04 (help.txt:1141) is this command plus small+ffirst; small-invariance means this fixture also verifies the H04 small fit ---
display _newline(2) "=== griliches first-stage, IID (H03/H04) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), ffirst
save_firststage, prefix(gril_fs) suffix(iid) outdir(`outdir')

// --- Cell 7b: D5a WITH DISCLOSURE -- cluster(tenure80) on the H03 base restores the K1=1 x cluster first-stage intersection lost with the retired sim_cluster cells (review finding, 2026-07-04); tenure80 (23 groups) is an arbitrary-but-real grouping variable unrelated to the model, not an economically motivated cluster structure ---
display _newline(2) "=== griliches first-stage, cluster(tenure80) (D5a) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), cluster(tenure80) ffirst
save_firststage, prefix(gril_fs) suffix(cl) outdir(`outdir')

// --- Cell 8 (DIAGNOSTICS) + OBJECT fs_gril_aw_iid: D5a -- deterministic synthetic aweight (replaces retired card weighted cells) ---
display _newline(2) "=== griliches first-stage, aweight, IID (D5a) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt) [aw=awt], ffirst savefirst
save_firststage, prefix(gril_fs_aw) suffix(iid) outdir(`outdir')
save_fs_fixture, prefix(fs_gril_aw_iid) eqname(_ivreg2_iq) outdir(`outdir')

// --- Cell 9 (DIAGNOSTICS) + OBJECT fs_gril_aw_robust: aweight, robust (D5a) ---
display _newline(2) "=== griliches first-stage, aweight, robust (D5a) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt) [aw=awt], robust ffirst savefirst
save_firststage, prefix(gril_fs_aw) suffix(hc1) outdir(`outdir')
save_fs_fixture, prefix(fs_gril_aw_robust) eqname(_ivreg2_iq) outdir(`outdir') nocoef


/*===========================================================================
  BLOCK 3: abdata
  Base model = H88 minus gmm2s (help.txt:1541): ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys)
  K1=3 so SW/AP conditional F genuinely differ from plain F and Shea != partial R2.
===========================================================================*/
use "`outdir'/_fs_ab_temp.dta", clear
quietly tsset id year

// --- Cell 10: H88 minus gmm2s, IID (help.txt:1541) ---
display _newline(2) "=== abdata first-stage, IID (H88 minus gmm2s) ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), ffirst
save_firststage, prefix(ab_fs) suffix(iid) outdir(`outdir')

// --- Cell 11: robust (help.txt:1541) ---
display _newline(2) "=== abdata first-stage, robust (H88 minus gmm2s) ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), robust ffirst
save_firststage, prefix(ab_fs) suffix(hc1) outdir(`outdir')

// --- Cell 12 (DIAGNOSTICS) + OBJECT fs_ab_cl: cluster(id), K1=3 multi-endogenous cluster-robust first stage (help.txt:1541) ---
display _newline(2) "=== abdata first-stage, cluster(id) (H88 minus gmm2s) ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), cluster(id) ffirst savefirst
save_firststage, prefix(ab_fs) suffix(cl) outdir(`outdir')
save_fs_fixture, prefix(fs_ab_cl) eqname(_ivreg2_w) outdir(`outdir')

// --- Cell 13: D5a -- weighted x cluster x first-stage intersection on a clean base (replaces the retired card M=2 smsa66 weighted-cluster cells) ---
display _newline(2) "=== abdata first-stage, aweight + cluster(id) (D5a) ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys) [aw=abwt], cluster(id) ffirst
save_firststage, prefix(ab_fs_aw) suffix(cl) outdir(`outdir')


/*===========================================================================
  Clean up temp files and done
===========================================================================*/
capture erase "`outdir'/_fs_mroz_temp.dta"
capture erase "`outdir'/_fs_gril_temp.dta"
capture erase "`outdir'/_fs_ab_temp.dta"
display _newline(2) "=== All first-stage fixtures generated ==="
display "Output directory: `outdir'"
