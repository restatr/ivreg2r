/*===========================================================================
  generate-stdp-fixtures.do
  -------------------------
  Generates CSV benchmark fixtures for predict(..., se.fit = TRUE) parity
  testing.

  Canonical base: mroz 2SLS baseline H31 (Stata ivreg2 help.txt line 1274).
  `predict, stdp` is Stata post-estimation with no worked example in the
  help file, so cells are option-variation on the canonical base.

  For each ivreg2 specification, saves fitted values and prediction SEs
  (Stata's "predict, stdp") for (a) all observations and (b) first 10 obs
  as a "newdata" subset.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-stdp-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

/*===========================================================================
  Load mroz data
===========================================================================*/
capture use "../validation/data/mroz.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/wooldridge/mroz.dta, clear
    if _rc {
        display as error "Could not load mroz dataset (no local cache, no network)."
        exit 601
    }
}


/*---------------------------------------------------------------------------
  Helper program: save fitted values and prediction SEs to CSV.
  Saves two files:
    - *_stdp_<suffix>.csv       (all obs: obs, fit, se_fit)
    - *_stdp_<suffix>_sub.csv   (first 10 obs only)
  The obs column (original row index) is not read by the tests, which
  compare positionally; it is kept for row-alignment debugging.
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
        format fit se_fit %21.0g
        export delimited using "`outdir'/`prefix'_stdp_`suffix'.csv", ///
            replace datafmt

        // Subset: first 10 observations (format set above persists)
        keep if _n <= 10
        export delimited using "`outdir'/`prefix'_stdp_`suffix'_sub.csv", ///
            replace datafmt
        restore
    }
end


/*===========================================================================
  FIXTURE 1: mroz IID
  Canonical base H31 (help.txt:1274); predict, stdp per D5a (no worked
  example in the help file).
  Model: ivreg2 lwage exper expersq (educ = age kidslt6 kidsge6)
===========================================================================*/
display _newline(2) "=== stdp: mroz iid ==="
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6)
save_stdp_results, prefix(mroz) suffix(iid) outdir(`outdir')


/*===========================================================================
  FIXTURE 2: mroz robust
  H41/H46 robust spec (help.txt:1325/1350) on the H31 base.
===========================================================================*/
display _newline(2) "=== stdp: mroz robust ==="
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6), robust
save_stdp_results, prefix(mroz) suffix(robust) outdir(`outdir')


/*===========================================================================
  FIXTURE 3: mroz robust small
  D5a option-variation (small) on the H41 robust spec; verifies the
  small-sample e(V) convention propagates into stdp.
===========================================================================*/
display _newline(2) "=== stdp: mroz robust small ==="
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6), robust small
save_stdp_results, prefix(mroz) suffix(robust_small) outdir(`outdir')


/*===========================================================================
  Done
===========================================================================*/
display _newline(2) "=== All stdp fixtures generated ==="
display "Output directory: `outdir'"
