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
export delimited using "`outdir'/mroz_data.csv", replace

// --- Wooldridge wagepan panel data ---
capture bcuse wagepan, clear
if _rc != 0 {
    display as error "Could not load wagepan dataset via bcuse."
    display as error "Install bcuse (ssc install bcuse) and rerun."
    exit 601
}
export delimited using "`outdir'/wagepan_data.csv", replace

display _newline
display "Vignette datasets exported successfully."
