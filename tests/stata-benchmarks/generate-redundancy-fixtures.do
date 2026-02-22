/*===========================================================================
  generate-redundancy-fixtures.do
  --------------------------------
  Generates CSV benchmark fixtures for instrument redundancy testing
  (Ticket P1).

  Redundancy test: H0 = specified excluded instruments have zero
  explanatory power for endogenous regressors, conditional on maintained
  instruments.  Stata stores e(redstat), e(redp), e(reddf), e(redlist).

  Card model: ivreg2 lwage exper expersq (educ = nearc2 nearc4)
  Multi-endo model: ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4)
  HAC model: ivreg2 y w (x = z1 z2), kernel() bw()

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-redundancy-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"


/*===========================================================================
  Load Card data FIRST — bcuse calls "clear all" which drops programs
===========================================================================*/
capture bcuse card, clear
if _rc != 0 {
    capture use "`outdir'/_card_temp.dta", clear
    if _rc != 0 {
        display as error "Could not load Card dataset."
        exit 601
    }
}
save "`outdir'/_card_redundancy_temp.dta", replace


/*---------------------------------------------------------------------------
  Helper program: save redundancy-specific diagnostics to CSV
  Defined AFTER bcuse to avoid being dropped by "clear all".
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
  Reload Card data
===========================================================================*/
use "`outdir'/_card_redundancy_temp.dta", clear


/*===========================================================================
  FIXTURE SET 1: Card, redundant(nearc2) — IID, HC0, HC1 small, cluster
  Model: ivreg2 lwage exper expersq (educ = nearc2 nearc4)
===========================================================================*/
display _newline(2) "=== redundant(nearc2), IID ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), redundant(nearc2)
save_redundancy_results, prefix(card_overid) suffix(nearc2_iid) outdir(`outdir')

display _newline(2) "=== redundant(nearc2), HC0 ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), redundant(nearc2) robust
save_redundancy_results, prefix(card_overid) suffix(nearc2_hc0) outdir(`outdir')

display _newline(2) "=== redundant(nearc2), HC1 small ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), redundant(nearc2) robust small
save_redundancy_results, prefix(card_overid) suffix(nearc2_hc1_small) outdir(`outdir')

display _newline(2) "=== redundant(nearc2), cluster ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), redundant(nearc2) cluster(smsa)
save_redundancy_results, prefix(card_overid) suffix(nearc2_cluster) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 2: Card, redundant(nearc4) — IID (different instrument)
===========================================================================*/
display _newline(2) "=== redundant(nearc4), IID ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), redundant(nearc4)
save_redundancy_results, prefix(card_overid) suffix(nearc4_iid) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 3: Card, redundant(nearc2 nearc4) — IID (test all excluded)
===========================================================================*/
display _newline(2) "=== redundant(nearc2 nearc4), IID ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), redundant(nearc2 nearc4)
save_redundancy_results, prefix(card_overid) suffix(both_iid) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 4: Card, redundant(nearc2) + dofminus — IID and cluster
===========================================================================*/
display _newline(2) "=== redundant(nearc2), IID, dofminus(1) ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), redundant(nearc2) dofminus(1)
save_redundancy_results, prefix(card_overid_dof) suffix(nearc2_iid) outdir(`outdir')

display _newline(2) "=== redundant(nearc2), cluster, dofminus(1) ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), redundant(nearc2) cluster(smsa) dofminus(1)
save_redundancy_results, prefix(card_overid_dof) suffix(nearc2_cluster) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 5: Card, redundant(nearc2) + aweight — IID
===========================================================================*/
display _newline(2) "=== redundant(nearc2), IID, aweight ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4) [aw=weight], redundant(nearc2)
save_redundancy_results, prefix(card_overid_aw) suffix(nearc2_iid) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 6: Card, redundant(nearc2) + gmm2s — robust
===========================================================================*/
display _newline(2) "=== redundant(nearc2), gmm2s robust ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), redundant(nearc2) gmm2s robust
save_redundancy_results, prefix(card_overid_gmm) suffix(nearc2_robust) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 7: sim_multi_endo, redundant — IID and HC1
  Model: ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4)
  K1=2 path
===========================================================================*/
import delimited using "`outdir'/sim_multi_endo_data.csv", clear

display _newline(2) "=== sim_multi_endo redundant(z1), IID ==="
ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4), redundant(z1)
save_redundancy_results, prefix(sim_multi_endo) suffix(z1_iid) outdir(`outdir')

display _newline(2) "=== sim_multi_endo redundant(z1), HC1 ==="
ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4), redundant(z1) robust
save_redundancy_results, prefix(sim_multi_endo) suffix(z1_hc1) outdir(`outdir')

display _newline(2) "=== sim_multi_endo redundant(z1 z2), IID ==="
ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4), redundant(z1 z2)
save_redundancy_results, prefix(sim_multi_endo) suffix(z1z2_iid) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 8: ts_hac, redundant — HAC path
  Model: ivreg2 y w (x = z1 z2), kernel(bartlett) bw(3) robust
===========================================================================*/
import delimited using "`outdir'/ts_hac_data.csv", clear
tsset t

display _newline(2) "=== ts_hac redundant(z1), Bartlett bw3 ==="
ivreg2 y w (x = z1 z2), redundant(z1) bw(3) kernel(bartlett) robust
save_redundancy_results, prefix(ts_hac) suffix(z1_bartlett_bw3) outdir(`outdir')


/*===========================================================================
  Clean up temp file and done
===========================================================================*/
capture erase "`outdir'/_card_redundancy_temp.dta"
display _newline(2) "=== All redundancy fixtures generated ==="
display "Output directory: `outdir'"
