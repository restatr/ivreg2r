/*===========================================================================
  generate-orthog-fixtures.do
  ---------------------------
  Generates CSV benchmark fixtures for orthogonality (instrument-subset)
  C-statistic testing (Ticket J1).

  Card model: ivreg2 lwage exper expersq (educ = nearc2 nearc4)
  Overidentified: L-K = 1 overid df.

  Test configurations:
    - orthog(nearc2): test 1 of 2 excluded instruments
    - orthog(nearc2 nearc4): test both (restricted model underidentified)

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-orthog-fixtures.do
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
save "`outdir'/_card_orthog_temp.dta", replace


/*---------------------------------------------------------------------------
  Helper program: save orthog-specific diagnostics to CSV
  Defined AFTER bcuse to avoid being dropped by "clear all".
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
  Reload Card data
===========================================================================*/
use "`outdir'/_card_orthog_temp.dta", clear


/*===========================================================================
  FIXTURE SET 1: Card, orthog(nearc2) — IID and HC1
  Model: ivreg2 lwage exper expersq (educ = nearc2 nearc4)
===========================================================================*/
display _newline(2) "=== orthog(nearc2), IID ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), orthog(nearc2)
save_orthog_results, prefix(card_orthog1) suffix(iid) outdir(`outdir')

display _newline(2) "=== orthog(nearc2), IID small ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), orthog(nearc2) small
save_orthog_results, prefix(card_orthog1) suffix(iid_small) outdir(`outdir')

display _newline(2) "=== orthog(nearc2), HC1 ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), orthog(nearc2) robust
save_orthog_results, prefix(card_orthog1) suffix(hc1) outdir(`outdir')

display _newline(2) "=== orthog(nearc2), HC1 small ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), orthog(nearc2) robust small
save_orthog_results, prefix(card_orthog1) suffix(hc1_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 2: Card, orthog(nearc2 nearc4) — test both
  Restricted model underidentified → Stata sets stat=0
===========================================================================*/
display _newline(2) "=== orthog(nearc2 nearc4), IID ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), orthog(nearc2 nearc4)
save_orthog_results, prefix(card_orthog2) suffix(iid) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 3: Card, orthog(nearc2) with dofminus(1) — IID and HC1
===========================================================================*/
display _newline(2) "=== orthog(nearc2), IID, dofminus(1) ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), orthog(nearc2) dofminus(1)
save_orthog_results, prefix(card_orthog1_dof) suffix(iid) outdir(`outdir')

display _newline(2) "=== orthog(nearc2), IID small, dofminus(1) ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), orthog(nearc2) dofminus(1) small
save_orthog_results, prefix(card_orthog1_dof) suffix(iid_small) outdir(`outdir')

display _newline(2) "=== orthog(nearc2), HC1, dofminus(1) ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), orthog(nearc2) dofminus(1) robust
save_orthog_results, prefix(card_orthog1_dof) suffix(hc1) outdir(`outdir')

display _newline(2) "=== orthog(nearc2), HC1 small, dofminus(1) ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), orthog(nearc2) dofminus(1) robust small
save_orthog_results, prefix(card_orthog1_dof) suffix(hc1_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 4: sim_cluster data, orthog(z1) — cluster-robust
  Model: ivreg2 y x1 (endo1 = z1 z2), cluster(cluster_id)
  Uses pre-generated sim_cluster_data.csv (50 clusters, 20 obs each)
===========================================================================*/
import delimited using "`outdir'/sim_cluster_data.csv", clear

display _newline(2) "=== sim_cluster orthog(z1), cluster ==="
ivreg2 y x1 (endo1 = z1 z2), orthog(z1) cluster(cluster_id)
save_orthog_results, prefix(sim_cluster_orthog1) suffix(cl) outdir(`outdir')

display _newline(2) "=== sim_cluster orthog(z1), cluster small ==="
ivreg2 y x1 (endo1 = z1 z2), orthog(z1) cluster(cluster_id) small
save_orthog_results, prefix(sim_cluster_orthog1) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 5: sim_twoway data, orthog(z1) — two-way cluster
  Model: ivreg2 y x1 (endo1 = z1 z2), cluster(firm_id year_id)
  Uses pre-generated sim_twoway_data.csv (25 firms, 10 years)
===========================================================================*/
import delimited using "`outdir'/sim_twoway_data.csv", clear

display _newline(2) "=== sim_twoway orthog(z1), cl2 ==="
ivreg2 y x1 (endo1 = z1 z2), orthog(z1) cluster(firm_id year_id)
save_orthog_results, prefix(sim_twoway_orthog1) suffix(cl2) outdir(`outdir')

display _newline(2) "=== sim_twoway orthog(z1), cl2 small ==="
ivreg2 y x1 (endo1 = z1 z2), orthog(z1) cluster(firm_id year_id) small
save_orthog_results, prefix(sim_twoway_orthog1) suffix(cl2_small) outdir(`outdir')


/*===========================================================================
  Clean up temp file and done
===========================================================================*/
capture erase "`outdir'/_card_orthog_temp.dta"
display _newline(2) "=== All orthog fixtures generated ==="
display "Output directory: `outdir'"
