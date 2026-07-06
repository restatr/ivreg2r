/*===========================================================================
  generate-stored-matrices-fixtures.do
  ------------------------------------
  Generates e(S) and e(W) matrix fixtures for the fit$S / fit$W export
  (ticket F1). One block per canonical configuration; canonical base specs
  with a single option varying (fixture-design rule).

  Data: all datasets load from the cached .dta files of record under ../validation/data/ from which the bundled data() objects are built bit-identically (M-27 re-base 2026-07-06 retired the griliches_data.csv, phillips_data.csv, and stockwatson_data.csv snapshots). stockwatson's derived columns are constructed in Stata from macrodat.dta EXACTLY as pkg/data-raw/stockwatson.R builds them (the same construction generate-cue-fixtures.do uses), so only the R data-raw script and one Stata construction recipe remain, both from the same float32 macrodat source. All loads use `use` (not bcuse), so they are safe after this file's program define blocks.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-stored-matrices-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

local outdir "tests/stata-benchmarks/fixtures"

/*---------------------------------------------------------------------------
  Helper: save a matrix to CSV (same as generate-wmatrix-fixtures.do)
---------------------------------------------------------------------------*/
capture program drop save_matrix
program define save_matrix
    syntax, mat(string) outfile(string)
    matrix M = `mat'
    local mr = rowsof(M)
    local mc = colsof(M)
    quietly {
        preserve
        clear
        set obs `mr'
        forvalues j = 1/`mc' {
            gen double v`j' = .
            forvalues i = 1/`mr' {
                replace v`j' = M[`i', `j'] in `i'
            }
        }
        export delimited using "`outfile'", replace
        restore
    }
end

/*---------------------------------------------------------------------------
  Helper: save a matrix's column names to a one-column CSV.
  The R tests map "_cons" -> "(Intercept)" and align by name, so no
  hand-coded order vectors are needed.
---------------------------------------------------------------------------*/
capture program drop save_matrix_names
program define save_matrix_names
    syntax, mat(string) outfile(string)
    matrix M = `mat'
    local cn : colnames M
    tempname fh
    file open `fh' using "`outfile'", write replace
    file write `fh' "name" _n
    foreach c of local cn {
        file write `fh' "`c'" _n
    }
    file close `fh'
end

/*===========================================================================
  1-3. Mroz: 2SLS iid / 2SLS robust / LIML (help.txt:1274 baseline)
       ivreg2 drops the 325 missing-lwage rows -> N = 428
===========================================================================*/
capture use "../validation/data/mroz.dta", clear
if _rc | _N == 0 {
    capture use http://fmwww.bc.edu/ec-p/data/wooldridge/mroz.dta, clear
}
quietly describe
if r(N) == 0 | r(k) == 0 {
    display as error "Could not load mroz (no local cache, bc.edu fetch failed)."
    display as error "Regenerate the cache by running the data-export chunk of validation/validate-helpfile.qmd (or: bcuse mroz, then save ../validation/data/mroz.dta)."
    exit 601
}

ivreg2 lwage exper expersq (educ = age kidslt6 kidsge6)
matrix S1 = e(S)
matrix W1 = e(W)
save_matrix, mat(S1) outfile("`outdir'/mroz_iv_eS_iid.csv")
save_matrix, mat(W1) outfile("`outdir'/mroz_iv_eW_iid.csv")
save_matrix_names, mat(S1) outfile("`outdir'/mroz_iv_inames.csv")

ivreg2 lwage exper expersq (educ = age kidslt6 kidsge6), robust
matrix S2 = e(S)
matrix W2 = e(W)
save_matrix, mat(S2) outfile("`outdir'/mroz_iv_eS_robust.csv")
save_matrix, mat(W2) outfile("`outdir'/mroz_iv_eW_robust.csv")

ivreg2 lwage exper expersq (educ = age kidslt6 kidsge6), liml
matrix S3 = e(S)
save_matrix, mat(S3) outfile("`outdir'/mroz_iv_eS_liml.csv")
* e(W) is not posted for LIML ("No weighting matrix defined")
capture matrix W3 = e(W)
display "LIML e(W) capture rc (expect nonzero): " _rc

/*===========================================================================
  4. Mroz: OLS robust
===========================================================================*/
ivreg2 lwage exper expersq, robust
matrix S4 = e(S)
matrix W4 = e(W)
save_matrix, mat(S4) outfile("`outdir'/mroz_ols_eS_robust.csv")
save_matrix, mat(W4) outfile("`outdir'/mroz_ols_eW_robust.csv")
save_matrix_names, mat(S4) outfile("`outdir'/mroz_ols_inames.csv")

/*===========================================================================
  5. Griliches: gmm2s robust (help.txt:1154, minus the xi year dummies --
     deviation disclosed; the year-dummy variant lives in the vignettes)
===========================================================================*/
capture use "../validation/data/griliches76.dta", clear
if _rc | _N == 0 {
    capture use http://fmwww.bc.edu/ec-p/data/hayashi/griliches76.dta, clear
}
quietly describe
if r(N) == 0 | r(k) == 0 {
    display as error "Could not load griliches76 (no local cache, bc.edu fetch failed)."
    exit 601
}

ivreg2 lw s expr tenure rns smsa (iq = med kww age mrt), gmm2s robust
matrix S5 = e(S)
matrix W5 = e(W)
save_matrix, mat(S5) outfile("`outdir'/griliches_eS_gmm2s_robust.csv")
save_matrix, mat(W5) outfile("`outdir'/griliches_eW_gmm2s_robust.csv")
save_matrix_names, mat(S5) outfile("`outdir'/griliches_inames.csv")

/*===========================================================================
  6-7. Phillips: AC (help.txt:1501) and HAC (help.txt:1511), Bartlett bw 3
===========================================================================*/
capture use "../validation/data/phillips.dta", clear
if _rc | _N == 0 {
    capture use http://fmwww.bc.edu/ec-p/data/wooldridge/phillips.dta, clear
}
quietly describe
if r(N) == 0 | r(k) == 0 {
    display as error "Could not load phillips (no local cache, bc.edu fetch failed)."
    exit 601
}
tsset year

ivreg2 cinf unem, bw(3) kernel(bartlett)
matrix S6 = e(S)
matrix W6 = e(W)
save_matrix, mat(S6) outfile("`outdir'/phillips_eS_ac_bw3.csv")
save_matrix, mat(W6) outfile("`outdir'/phillips_eW_ac_bw3.csv")
save_matrix_names, mat(S6) outfile("`outdir'/phillips_inames.csv")

ivreg2 cinf unem, robust bw(3) kernel(bartlett)
matrix S7 = e(S)
matrix W7 = e(W)
save_matrix, mat(S7) outfile("`outdir'/phillips_eS_hac_bw3.csv")
save_matrix, mat(W7) outfile("`outdir'/phillips_eW_hac_bw3.csv")

/*===========================================================================
  8. Stock-Watson: CUE, HAC Bartlett bw 5 (BSS 2007 SS4.3, p. 480)
===========================================================================*/
capture use "../validation/data/macrodat.dta", clear
if _rc | _N == 0 {
    capture use http://fmwww.bc.edu/ec-p/data/stockwatson/macrodat.dta, clear
}
quietly describe
if r(N) == 0 | r(k) == 0 {
    display as error "Could not load stockwatson (macrodat) (no local cache, bc.edu fetch failed)."
    exit 601
}
* Construct the derived columns EXACTLY as pkg/data-raw/stockwatson.R does, in
* double precision, so the Stata data is numerically coherent with the bundled
* data(stockwatson) the R test fits on (F4 precision-lesson discipline: both
* sides compute in double from the same float32 macrodat source).
sort date
tsset date
gen double inf = 100 * ln(CPI / L4.CPI)
gen double ggdp = 100 * ln(GDP / L4.GDP)
gen double dinf = inf - L.inf
gen double ggdp_2 = L2.ggdp
gen double TBILL_1 = L.TBILL
gen double ER_1 = L.ER
gen double TBON_1 = L.TBON

ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue robust kernel(bartlett) bw(5)
matrix S8 = e(S)
matrix W8 = e(W)
save_matrix, mat(S8) outfile("`outdir'/stockwatson_eS_cue.csv")
save_matrix, mat(W8) outfile("`outdir'/stockwatson_eW_cue.csv")
save_matrix_names, mat(S8) outfile("`outdir'/stockwatson_inames.csv")

display "generate-stored-matrices-fixtures.do complete"
