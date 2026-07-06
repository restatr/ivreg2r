/*===========================================================================
  generate-hols-fixtures.do
  -------------------------
  Generates CSV benchmark fixtures for the empty-endogenous / HOLS syntax
  (ticket F3): Stata's `(=z1 z2)` form, R's `y ~ exog | 0 | instruments`.

  Covers the four help-file acceptance examples (ivreg2-help.txt):
    H34 (~1290): ivreg2 lwage exper expersq educ (=age kidslt6 kidsge6), orthog(educ)
    H35 (~1296): ivreg2 lwage exper expersq educ (=age kidslt6 kidsge6), gmm2s
    H48 (~1356): ivreg2 educ exper expersq (=age kidslt6 kidsge6) if e(sample), robust
    H52 (~1375): ivreg2 educ exper expersq kidslt6 kidsge6 (=age) if e(sample), robust
  plus plain iid/robust variants and the companion IV fits whose idstat /
  redstat the H48/H52 J statistics reproduce.

  Probe blocks (logged, no fixtures): liml with (=z); saverf with (=z).

  Data: loads from ../validation/data/mroz.dta, the cached file of record from which the bundled data(mroz) is built bit-identically (753 rows incl. missing lwage; ivreg2 drops those rows -> N = 428). mroz_data.csv was retired 2026-07-06.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-hols-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

local outdir "tests/stata-benchmarks/fixtures"

/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV
  (schema follows generate-psd-fixtures.do, plus orthog/redundancy columns)
---------------------------------------------------------------------------*/
capture program drop save_ivreg2_results
program define save_ivreg2_results
    syntax, prefix(string) suffix(string) outdir(string)

    local N = e(N)
    local K = e(rankxx)
    local L = e(rankzz)

    // --- Coefficients and SEs ---
    matrix b = e(b)
    matrix V = e(V)
    local names : colnames b
    local ncols = colsof(b)

    quietly {
        preserve
        clear
        set obs `ncols'
        gen str32 term = ""
        gen double estimate = .
        gen double std_error = .

        forvalues i = 1/`ncols' {
            local nm : word `i' of `names'
            replace term = "`nm'" in `i'
            replace estimate = b[1, `i'] in `i'
            replace std_error = sqrt(V[`i', `i']) in `i'
        }
        export delimited using "`outdir'/`prefix'_coef_`suffix'.csv", replace
        restore
    }

    // --- Full VCV matrix ---
    matrix V = e(V)
    quietly {
        preserve
        clear
        local vr = rowsof(V)
        local vc = colsof(V)
        set obs `vr'
        forvalues j = 1/`vc' {
            gen double v`j' = .
            forvalues i = 1/`vr' {
                replace v`j' = V[`i', `j'] in `i'
            }
        }
        export delimited using "`outdir'/`prefix'_vcov_`suffix'.csv", replace
        restore
    }

    // --- e(S) matrix ---
    matrix S = e(S)
    quietly {
        preserve
        clear
        local sr = rowsof(S)
        set obs `sr'
        forvalues j = 1/`sr' {
            gen double v`j' = .
            forvalues i = 1/`sr' {
                replace v`j' = S[`i', `j'] in `i'
            }
        }
        export delimited using "`outdir'/`prefix'_eS_`suffix'.csv", replace
        restore
    }

    // --- Diagnostics ---
    quietly {
        preserve
        clear
        set obs 1

        // Overidentification
        gen double overid_stat = e(sargan)
        if overid_stat == . {
            replace overid_stat = e(j)
        }
        gen double overid_p = e(sarganp)
        if overid_p == . {
            replace overid_p = e(jp)
        }
        gen int overid_df = e(sargandf)
        if overid_df == . {
            replace overid_df = e(jdf)
        }

        // Underidentification (Anderson / KP) — expect missing for (=z)
        gen double underid_stat = e(idstat)
        gen double underid_p = e(idp)
        gen int underid_df = e(iddf)

        // Weak identification — expect missing for (=z)
        gen double weak_id_cd_f = e(cdf)
        gen double weak_id_kp_f = e(widstat)

        // Orthogonality C-stat
        gen double cstat = e(cstat)
        gen double cstat_p = e(cstatp)
        gen int cstat_df = e(cstatdf)

        // Redundancy
        gen double redstat = e(redstat)
        gen double redstat_p = e(redp)
        gen int redstat_df = e(reddf)

        // Stock-Wright S — expect missing for (=z)
        gen double sw_stat = e(sstat)
        gen double sw_p = e(sstatp)
        gen int sw_df = e(sstatdf)

        // LIML lambda and Anderson-Rubin LR overid (liml configs only)
        gen double lambda = e(lambda)
        gen double arubin = e(arubin)
        gen double arubin_p = e(arubinp)

        // Model F
        gen double model_f = e(F)
        gen double model_f_p = e(Fp)
        gen int model_f_df1 = e(Fdf1)
        gen int model_f_df2 = e(Fdf2)

        // Counts
        gen int N = e(N)
        gen int K = `K'
        gen int L = `L'

        export delimited using "`outdir'/`prefix'_diagnostics_`suffix'.csv", replace
        restore
    }
end

/*---------------------------------------------------------------------------
  Helper: save e(S) column names to a one-column CSV
---------------------------------------------------------------------------*/
capture program drop save_s_names
program define save_s_names
    syntax, outfile(string)
    matrix Smat = e(S)
    local cn : colnames Smat
    tempname fh
    file open `fh' using "`outfile'", write replace
    file write `fh' "name" _n
    foreach c of local cn {
        file write `fh' "`c'" _n
    }
    file close `fh'
end

/*===========================================================================
  Load mroz (cached source of record; also feeds the bundled data(mroz))
===========================================================================*/
capture use "../validation/data/mroz.dta", clear
if _rc | _N == 0 {
    display as error "Cannot load ../validation/data/mroz.dta (the cached source of record)."
    display as error "Regenerate the cache by running the data-export chunk of validation/validate-helpfile.qmd (or: bcuse mroz, then save ../validation/data/mroz.dta)."
    exit 601
}

/*===========================================================================
  Canonical HOLS base: lwage exper expersq educ (=age kidslt6 kidsge6)
===========================================================================*/

// --- 1. Plain iid (OLS with surplus instruments; Sargan posted) ---
ivreg2 lwage exper expersq educ (=age kidslt6 kidsge6)
save_ivreg2_results, prefix(mroz_hols) suffix(iid) outdir(`outdir')
save_s_names, outfile("`outdir'/mroz_hols_inames.csv")

// --- 2. Robust (Hansen J) ---
ivreg2 lwage exper expersq educ (=age kidslt6 kidsge6), robust
save_ivreg2_results, prefix(mroz_hols) suffix(robust) outdir(`outdir')

// --- 3. H34: orthog(educ) C-test of an included exogenous regressor ---
ivreg2 lwage exper expersq educ (=age kidslt6 kidsge6), orthog(educ)
save_ivreg2_results, prefix(mroz_hols) suffix(orthog) outdir(`outdir')

// --- 4a. H35 exact command: bare gmm2s. With no VCE option the 2-step
//     weighting matrix is proportional to (Z'Z)^-1, so the estimates equal
//     OLS exactly (help.txt:344-358: "IV/2SLS, SEs consistent under
//     homoskedasticity"; subtitle "efficient for homoskedasticity only") ---
ivreg2 lwage exper expersq educ (=age kidslt6 kidsge6), gmm2s
save_ivreg2_results, prefix(mroz_hols) suffix(gmm2s_iid) outdir(`outdir')

// --- 4b. HOLS proper: gmm2s + robust, where the surplus conditions buy
//     heteroskedasticity-efficiency and the estimates differ from OLS ---
ivreg2 lwage exper expersq educ (=age kidslt6 kidsge6), gmm2s robust
save_ivreg2_results, prefix(mroz_hols) suffix(gmm2s) outdir(`outdir')

/*===========================================================================
  H48: J of the empty-endo educ equation = KP idstat of the IV model
===========================================================================*/

// Companion IV model (defines e(sample) = 428 lwage rows)
ivreg2 lwage exper expersq (educ = age kidslt6 kidsge6), robust
save_ivreg2_results, prefix(mroz_hols) suffix(h48_iv) outdir(`outdir')

ivreg2 educ exper expersq (=age kidslt6 kidsge6) if e(sample), robust
save_ivreg2_results, prefix(mroz_hols) suffix(h48_j) outdir(`outdir')

/*===========================================================================
  H52: J of the empty-endo equation = redundancy LM stat
===========================================================================*/

ivreg2 lwage exper expersq (educ = age kidslt6 kidsge6), robust redundant(age)
save_ivreg2_results, prefix(mroz_hols) suffix(h52_iv) outdir(`outdir')

ivreg2 educ exper expersq kidslt6 kidsge6 (=age) if e(sample), robust
save_ivreg2_results, prefix(mroz_hols) suffix(h52_j) outdir(`outdir')

/*===========================================================================
  5. LIML with (=z): coefficients equal OLS exactly (X'M_Z = 0 makes the
     k-class normal equations collapse to X'X for any k); lambda is the
     1x1 eigenvalue RSS(y~X)/RSS(y~Z) and e(arubin) = (N)ln(lambda) is the
     LR overid statistic of the surplus conditions.
===========================================================================*/
ivreg2 lwage exper expersq educ (=age kidslt6 kidsge6), liml
save_ivreg2_results, prefix(mroz_hols) suffix(liml) outdir(`outdir')

/*===========================================================================
  Probe blocks — logged only, no fixtures
===========================================================================*/

display "===== PROBE: saverf with (=z) — stores nothing (e(rfeq) empty) ====="
capture noisily ivreg2 lwage exper expersq educ (=age kidslt6 kidsge6), saverf
display "PROBE saverf rc: " _rc
capture noisily display "PROBE saverf e(rfeq): [" e(rfeq) "]"
estimates dir

display "generate-hols-fixtures.do complete"
