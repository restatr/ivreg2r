/*===========================================================================
  generate-stored-matrices-fixtures.do
  ------------------------------------
  Generates e(S) and e(W) matrix fixtures for the fit$S / fit$W export
  (ticket F1). One block per canonical configuration; canonical base specs
  with a single option varying (fixture-design rule).

  Data: CSV exports of the BUNDLED package datasets (written by R from the .rda files), so Stata and R compute on byte-identical inputs, except mroz: mroz_data.csv was retired 2026-07-06, and mroz now loads from ../validation/data/mroz.dta, the cached file of record from which the bundled data(mroz) is built bit-identically. No bcuse (the mroz load sits after this file's program define blocks).

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
import delimited "`outdir'/griliches_data.csv", clear case(preserve)

ivreg2 lw s expr tenure rns smsa (iq = med kww age mrt), gmm2s robust
matrix S5 = e(S)
matrix W5 = e(W)
save_matrix, mat(S5) outfile("`outdir'/griliches_eS_gmm2s_robust.csv")
save_matrix, mat(W5) outfile("`outdir'/griliches_eW_gmm2s_robust.csv")
save_matrix_names, mat(S5) outfile("`outdir'/griliches_inames.csv")

/*===========================================================================
  6-7. Phillips: AC (help.txt:1501) and HAC (help.txt:1511), Bartlett bw 3
===========================================================================*/
import delimited "`outdir'/phillips_data.csv", clear case(preserve)
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
import delimited "`outdir'/stockwatson_data.csv", clear case(preserve)
tsset date

ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1), cue robust kernel(bartlett) bw(5)
matrix S8 = e(S)
matrix W8 = e(W)
save_matrix, mat(S8) outfile("`outdir'/stockwatson_eS_cue.csv")
save_matrix, mat(W8) outfile("`outdir'/stockwatson_eW_cue.csv")
save_matrix_names, mat(S8) outfile("`outdir'/stockwatson_inames.csv")

display "generate-stored-matrices-fixtures.do complete"
