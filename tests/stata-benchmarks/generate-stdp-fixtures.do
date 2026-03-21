/*===========================================================================
  generate-stdp-fixtures.do
  -------------------------
  Generates CSV benchmark fixtures for predict(..., se.fit = TRUE) testing
  (Ticket Q1).

  For each ivreg2 specification, saves fitted values and prediction SEs
  (Stata's "predict, stdp") for (a) all observations and (b) first 10 obs
  as "newdata" subset.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-stdp-fixtures.do
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
save "`outdir'/_card_stdp_temp.dta", replace


/*---------------------------------------------------------------------------
  Helper program: save fitted values and prediction SEs to CSV.
  Saves two files:
    - *_stdp_<suffix>.csv       (all obs: obs, fit, se_fit)
    - *_stdp_<suffix>_sub.csv   (first 10 obs only)
---------------------------------------------------------------------------*/
capture program drop save_stdp_results
program define save_stdp_results
    syntax, prefix(string) suffix(string) outdir(string)

    // Predict fitted values and stdp for estimation sample
    tempvar yhat se
    predict double `yhat', xb
    predict double `se', stdp

    quietly {
        preserve
        gen long obs = _n
        keep if e(sample)
        keep obs `yhat' `se'
        rename `yhat' fit
        rename `se' se_fit
        export delimited using "`outdir'/`prefix'_stdp_`suffix'.csv", replace

        // Subset: first 10 observations
        keep if _n <= 10
        export delimited using "`outdir'/`prefix'_stdp_`suffix'_sub.csv", replace
        restore
    }
end


/*===========================================================================
  Reload Card data
===========================================================================*/
use "`outdir'/_card_stdp_temp.dta", clear


/*===========================================================================
  FIXTURE SET 1: Card overid, IID
  Model: ivreg2 lwage exper expersq (educ = nearc2 nearc4)
===========================================================================*/
display _newline(2) "=== stdp: Card overid, IID ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4)
save_stdp_results, prefix(card_overid) suffix(iid) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 2: Card overid, robust (HC0)
===========================================================================*/
display _newline(2) "=== stdp: Card overid, robust ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), robust
save_stdp_results, prefix(card_overid) suffix(hc0) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 3: Card overid, robust small (HC1)
===========================================================================*/
display _newline(2) "=== stdp: Card overid, robust small ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), robust small
save_stdp_results, prefix(card_overid) suffix(hc1_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 4: Card overid, cluster(smsa)
===========================================================================*/
display _newline(2) "=== stdp: Card overid, cluster ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), cluster(smsa)
save_stdp_results, prefix(card_overid) suffix(cluster) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 5: Card overid, cluster(smsa) small
===========================================================================*/
display _newline(2) "=== stdp: Card overid, cluster small ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), cluster(smsa) small
save_stdp_results, prefix(card_overid) suffix(cluster_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 6: Card overid, aweight
===========================================================================*/
display _newline(2) "=== stdp: Card overid, aweight ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4) [aw=weight]
save_stdp_results, prefix(card_overid) suffix(aw) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 7: Card overid, LIML
===========================================================================*/
display _newline(2) "=== stdp: Card overid, LIML ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), liml
save_stdp_results, prefix(card_overid) suffix(liml) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 8: Card overid, Fuller(1)
===========================================================================*/
display _newline(2) "=== stdp: Card overid, Fuller(1) ==="
ivreg2 lwage exper expersq (educ = nearc2 nearc4), fuller(1)
save_stdp_results, prefix(card_overid) suffix(fuller1) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 9: Card just-identified, IID
  Model: ivreg2 lwage exper expersq (educ = nearc4)
===========================================================================*/
display _newline(2) "=== stdp: Card justid, IID ==="
ivreg2 lwage exper expersq (educ = nearc4)
save_stdp_results, prefix(card_justid) suffix(iid) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 10: Card just-identified, robust
===========================================================================*/
display _newline(2) "=== stdp: Card justid, robust ==="
ivreg2 lwage exper expersq (educ = nearc4), robust
save_stdp_results, prefix(card_justid) suffix(hc0) outdir(`outdir')


/*===========================================================================
  Clean up temp file and done
===========================================================================*/
capture erase "`outdir'/_card_stdp_temp.dta"
display _newline(2) "=== All stdp fixtures generated ==="
display "Output directory: `outdir'"
