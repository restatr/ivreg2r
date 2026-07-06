/*===========================================================================
  generate-redundancy-fixtures.do
  --------------------------------
  Generates CSV benchmark fixtures for instrument redundancy testing (M-23, planning/22-spec-matrix.md), re-based off the Card/sim_multi_endo/ts_hac ad hoc data onto the canonical help-file / literature bases:

    - griliches76 (BSS 2007 pp. 490-493): the weak-instrument schooling-IQ wage model, redundant(mrt) and redundant(age mrt), used as the canonical anchor — BSS report redundant(mrt) chi2(1) = 0.002, p = 0.9665 on this exact spec.
    - abdata (H88, help.txt:1541): the Arellano-Bond employment model with first- and second-difference instruments, minus gmm2s, with redundant() added as a D5a option-variation to exercise multi-endogenous (K1=3) and cluster-robust redundancy.
    - phillips (M-19 base, help.txt:1489-1520): the HAC-VCE Phillips-curve model, with redundant() added as a D5a option-variation to exercise the HAC-redundancy code path.

  `small` variants are deliberately NOT generated for any cell in this file: Stata's e(redstat)/e(redp)/e(reddf) are invariant to the `small` option (verified byte-identical in the retired card fixtures hc0 vs hc1_small), so the R tests fit both small = FALSE and small = TRUE against each fixture here and assert both reproduce the same redundancy statistic.
  gmm2s cells are deliberately NOT generated: e(redstat) under gmm2s robust is byte-identical to 2SLS robust (verified in the retired card fixtures), so the R tests assert estimator-invariance as a fixture-free identity rather than pinning a duplicate fixture.

  Redundancy test: H0 = specified excluded instruments have zero explanatory power for endogenous regressors, conditional on maintained instruments.  Stata stores e(redstat), e(redp), e(reddf), e(redlist).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-redundancy-fixtures.do
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

// --- griliches76: local cache, then BC.edu web fallback (per the BSS 2007 pp. 490-493 base model) ---
capture use "../validation/data/griliches76.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/hayashi/griliches76.dta, clear
    if _rc {
        display as error "Could not load griliches76 dataset (no local cache, no network)."
        exit 601
    }
}
quietly xi i.year
// Deterministic synthetic analytic weight (M-12 precedent): mod(age,5)+1 gives a reproducible, non-degenerate aweight without inventing an economically motivated weight.
gen awt = mod(age,5)+1
save "`outdir'/_red_gril_temp.dta", replace

// --- abdata: local cache, then BC.edu web fallback (per the H88 base model, help.txt:1541) ---
capture use "../validation/data/abdata.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/macro/abdata.dta, clear
    if _rc {
        display as error "Could not load abdata dataset (no local cache, no network)."
        exit 601
    }
}
save "`outdir'/_red_ab_temp.dta", replace

// --- phillips: local cache, then BC.edu web fallback (per the M-19 base model, help.txt:1489-1520) ---
capture use "../validation/data/phillips.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/macro/phillips.dta, clear
    if _rc {
        display as error "Could not load phillips dataset (no local cache, no network)."
        exit 601
    }
}
save "`outdir'/_red_phil_temp.dta", replace


/*---------------------------------------------------------------------------
  Helper program: save redundancy-specific diagnostics to CSV
---------------------------------------------------------------------------*/
capture program drop save_redundancy_results
program define save_redundancy_results
    syntax, prefix(string) suffix(string) outdir(string)

    quietly {
        preserve
        clear
        set obs 1

        gen double redstat = .
        gen double redp = .
        gen double reddf = .
        gen str80 redlist = ""

        // Redundancy test
        capture replace redstat = e(redstat)
        capture replace redp = e(redp)
        capture replace reddf = e(reddf)
        capture replace redlist = "`e(redlist)'"

        export delimited using ///
            "`outdir'/`prefix'_redundancy_`suffix'.csv", replace
        restore
    }
end


/*===========================================================================
  BLOCK 1: griliches76
  Base model = the BSS 2007 pp. 490-493 weak-instrument specification:
    ivreg2 lw s expr tenure rns smsa _I* (iq=age mrt)
  (excluded instruments ONLY age and mrt -- the weak-ID demo spec, NOT the H03 4-instrument model)
===========================================================================*/
use "`outdir'/_red_gril_temp.dta", clear

// --- Canonical anchor: BSS 2007 p. 493 report redundant(mrt) chi2(1) = 0.002, p = 0.9665 on this spec -- the generated fixture should reproduce those values ---
display _newline(2) "=== griliches redundant(mrt), robust (BSS 2007 p. 493 anchor) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=age mrt), robust redundant(mrt)
save_redundancy_results, prefix(gril_red) suffix(mrt_robust) outdir(`outdir')

// --- D5a: BSS command minus robust -- IID canonical-correlations path ---
display _newline(2) "=== griliches redundant(mrt), IID ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=age mrt), redundant(mrt)
save_redundancy_results, prefix(gril_red) suffix(mrt_iid) outdir(`outdir')

// --- D5a: all excluded IVs tested (df=2); consumed by the R-side redundancy==underid K1=1 identity test ---
display _newline(2) "=== griliches redundant(age mrt), IID ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=age mrt), redundant(age mrt)
save_redundancy_results, prefix(gril_red) suffix(both_iid) outdir(`outdir')

// --- D5a: dofminus threading (redstat verified NOT dofminus-invariant in the retired card fixtures) ---
display _newline(2) "=== griliches redundant(mrt), IID, dofminus(1) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=age mrt), redundant(mrt) dofminus(1)
save_redundancy_results, prefix(gril_red_dof) suffix(mrt_iid) outdir(`outdir')

display _newline(2) "=== griliches redundant(mrt), robust, dofminus(1) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=age mrt), robust redundant(mrt) dofminus(1)
save_redundancy_results, prefix(gril_red_dof) suffix(mrt_robust) outdir(`outdir')

// --- D5a: deterministic synthetic aweight on the canonical base (M-12 precedent) -- weighted-redundancy path ---
display _newline(2) "=== griliches redundant(mrt), IID, aweight ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=age mrt) [aw=awt], redundant(mrt)
save_redundancy_results, prefix(gril_red_aw) suffix(mrt_iid) outdir(`outdir')


/*===========================================================================
  BLOCK 2: abdata
  Base model = help-file H88 (help.txt:1541), minus gmm2s:
    ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys)
===========================================================================*/
use "`outdir'/_red_ab_temp.dta", clear
quietly tsset id year

// --- D5a: K1=3 multi-endogenous redundancy, IID path (df = K1*L_tested = 3); replaces the retired sim_multi_endo cells ---
display _newline(2) "=== abdata redundant(d2.ys), IID ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), redundant(d2.ys)
save_redundancy_results, prefix(ab_red) suffix(d2ys_iid) outdir(`outdir')

// --- D5a: cluster-robust redundancy on a real canonical base; replaces the retired card cluster(smsa) M=2 anti-pattern cell ---
display _newline(2) "=== abdata redundant(d2.ys), cluster(id) ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), cluster(id) redundant(d2.ys)
save_redundancy_results, prefix(ab_red) suffix(d2ys_cl) outdir(`outdir')

// --- D5a: K1=3 with plain heteroskedasticity-robust (non-cluster) VCE; restores the retired sim_multi_endo z1_hc1 intersection (multi-endo x HC) on a canonical base ---
display _newline(2) "=== abdata redundant(d2.ys), robust ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), robust redundant(d2.ys)
save_redundancy_results, prefix(ab_red) suffix(d2ys_robust) outdir(`outdir')

// --- D5a: K1=3 with L_tested=2 (df = 6); restores the retired sim_multi_endo z1z2 intersection (multi-endo x multi-instrument) on a canonical base ---
display _newline(2) "=== abdata redundant(d2.k d2.ys), IID ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), redundant(d2.k d2.ys)
save_redundancy_results, prefix(ab_red) suffix(d2kd2ys_iid) outdir(`outdir')

// --- dofminus x cluster: no separate fixture is generated because Stata ignores dofminus in the cluster-VCE redundancy stat — `cluster(id) dofminus(1) redundant(d2.ys)` was verified byte-identical to the d2ys_cl cell at re-base (2026-07-04); the R tests pin this behavior by fitting dofminus = 1L against the d2ys_cl fixture ---


/*===========================================================================
  BLOCK 3: phillips
  Base model = the M-19 phillips IV model, help.txt:1489-1520
===========================================================================*/
use "`outdir'/_red_phil_temp.dta", clear
quietly tsset year

// --- D5a: HAC-redundancy option-variation on the phillips base (replaces the retired synthetic ts_hac cell) ---
display _newline(2) "=== phillips redundant(l3.unem), HAC robust bw(3) Bartlett ==="
ivreg2 cinf (unem = l(1/3).unem), robust bw(3) kernel(bartlett) redundant(l3.unem)
save_redundancy_results, prefix(phil_red) suffix(l3_bartlett_bw3) outdir(`outdir')


/*===========================================================================
  Clean up temp files and done
===========================================================================*/
capture erase "`outdir'/_red_gril_temp.dta"
capture erase "`outdir'/_red_ab_temp.dta"
capture erase "`outdir'/_red_phil_temp.dta"
display _newline(2) "=== All redundancy fixtures generated ==="
display "Output directory: `outdir'"
