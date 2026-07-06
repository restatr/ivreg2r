/*===========================================================================
  generate-vignette-data.do
  -------------------------
  Exports Mroz and wagepan datasets for use in ivreg2r vignettes.
  Run in Stata with bcuse installed (ssc install bcuse).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-vignette-data.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

// --- Mroz (1987) female labor supply data ---
capture bcuse mroz, clear
if _rc != 0 {
    display as error "Could not load mroz dataset via bcuse."
    display as error "Install bcuse (ssc install bcuse) and rerun."
    exit 601
}
// Full-precision export: the default `export delimited` writes numeric
// columns at 8-digit precision, which truncates float32 bit patterns
// (see ticket F4 / the CSV float-truncation gotcha in CLAUDE.md).
quietly ds, has(type numeric)
format `r(varlist)' %21.0g
export delimited using "`outdir'/mroz_data.csv", replace datafmt

// --- Wooldridge wagepan panel data ---
capture bcuse wagepan, clear
if _rc != 0 {
    display as error "Could not load wagepan dataset via bcuse."
    display as error "Install bcuse (ssc install bcuse) and rerun."
    exit 601
}
// Full-precision export: the default `export delimited` writes numeric
// columns at 8-digit precision, which truncates float32 bit patterns
// (see ticket F4 / the CSV float-truncation gotcha in CLAUDE.md).
// WARNING (2026-07-04): do NOT regenerate wagepan_data.csv until the
// CUE convergence criterion is investigated at M-17 -- on the exact
// float32 data the wp_sw_cue fit (test-vcov-sw.R) warns of CUE
// non-convergence (estimates still match the fixture). See the
// data-CSV retirement item in planning/22-spec-matrix.md.
quietly ds, has(type numeric)
format `r(varlist)' %21.0g
export delimited using "`outdir'/wagepan_data.csv", replace datafmt

display _newline
display "Vignette datasets exported successfully."
