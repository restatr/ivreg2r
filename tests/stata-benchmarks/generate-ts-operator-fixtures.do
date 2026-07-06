/*===========================================================================
  generate-ts-operator-fixtures.do
  --------------------------------
  Generates CSV benchmark fixtures for time-series operator (l()/d())
  testing in ivreg2r — the 9 help-file examples that use Stata ts operators
  in the varlist (help-file-specs.md H72-H76, H83-H85, H88).

  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Data sources: local cache in ../validation/data/ (built by
  validation/validate-helpfile.qmd); falls back to webuse / BC server.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-ts-operator-fixtures.do
===========================================================================*/

clear all
set more off
set sortseed 12345  // pin sort-tie order: Stata sorts place ties in random order and the sort RNG state persists across do-files in a session (see CLAUDE.md, Stata gotchas)
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV
  (Same as generate-liml-fixtures.do: includes kclass/lambda/fuller and
  C-stat columns, which the klein LIML-family and phillips orthog examples
  need. First-stage matrix output omitted — not consumed by these tests.)
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
        local dim = rowsof(V)
        svmat double V
        gen str32 term = ""
        forvalues i = 1/`dim' {
            local nm : word `i' of `names'
            replace term = "`nm'" in `i'
        }
        order term
        forvalues i = 1/`dim' {
            local nm : word `i' of `names'
            local cnm = subinstr("`nm'", ".", "_", .)
            rename V`i' vcov_`cnm'
        }
        export delimited using "`outdir'/`prefix'_vcov_`suffix'.csv", replace
        restore
    }

    // --- Diagnostics ---
    quietly {
        preserve
        clear
        set obs 1

        // Model stats
        gen double N = `N'
        gen double K = `K'
        gen double L = `L'
        gen double r2 = e(r2)
        gen double r2_a = e(r2_a)
        gen double r2u = e(r2u)
        gen double r2c = e(r2c)
        gen double rss = e(rss)
        gen double rmse = e(rmse)
        gen double F_stat = e(F)
        gen double F_p = e(Fp)
        gen double F_df1 = e(Fdf1)
        gen double F_df2 = e(Fdf2)

        // Initialize conditional columns to missing
        gen double sargan = .
        gen double sarganp = .
        gen double sargandf = .
        gen double j = .
        gen double jp = .
        gen double jdf = .
        gen double idstat = .
        gen double idp = .
        gen double iddf = .
        gen double cdf = .
        gen double widstat = .
        gen double arf = .
        gen double arfp = .
        gen double archi2 = .
        gen double archi2p = .
        gen double ardf = .
        gen double ardf_r = .
        gen double cstat = .
        gen double cstatp = .
        gen double cstatdf = .
        gen double sstat = .
        gen double sstatp = .
        gen double sstatdf = .
        gen double N_clust = .

        capture replace sargan = e(sargan)
        capture replace sarganp = e(sarganp)
        capture replace sargandf = e(sargandf)
        capture replace j = e(j)
        capture replace jp = e(jp)
        capture replace jdf = e(jdf)
        capture replace idstat = e(idstat)
        capture replace idp = e(idp)
        capture replace iddf = e(iddf)
        capture replace cdf = e(cdf)
        capture replace widstat = e(widstat)
        capture replace arf = e(arf)
        capture replace arfp = e(arfp)
        capture replace archi2 = e(archi2)
        capture replace archi2p = e(archi2p)
        capture replace ardf = e(ardf)
        capture replace ardf_r = e(ardf_r)
        capture replace cstat = e(cstat)
        capture replace cstatp = e(cstatp)
        capture replace cstatdf = e(cstatdf)
        capture replace sstat = e(sstat)
        capture replace sstatp = e(sstatp)
        capture replace sstatdf = e(sstatdf)
        capture replace N_clust = e(N_clust)

        // LIML / k-class metadata
        gen double kclass = .
        gen double lambda = .
        gen double fuller = .
        capture replace kclass = e(kclass)
        capture replace lambda = e(lambda)
        capture replace fuller = e(fuller)

        // Anderson-Rubin LIML overid
        gen double arubin = .
        gen double arubinp = .
        gen double arubin_lin = .
        gen double arubin_linp = .
        gen double arubindf = .
        capture replace arubin = e(arubin)
        capture replace arubinp = e(arubinp)
        capture replace arubin_lin = e(arubin_lin)
        capture replace arubin_linp = e(arubin_linp)
        capture replace arubindf = e(arubindf)

        export delimited using ///
            "`outdir'/`prefix'_diagnostics_`suffix'.csv", replace
        restore
    }
end

/*---------------------------------------------------------------------------
  Phillips data (H83-H85): ivreg2 cinf (unem = l(1/3).unem), ...
  Pure time series, tsset year. R side uses the bundled phillips dataset.
---------------------------------------------------------------------------*/
capture use "../validation/data/phillips.dta", clear
if _rc {
    use http://fmwww.bc.edu/ec-p/data/macro/phillips.dta, clear
}
tsset year

// --- H83 (help.txt:1520): AC, Bartlett (default kernel), bw(3) ---
ivreg2 cinf (unem = l(1/3).unem), bw(3)
save_ivreg2_results, prefix(tsop_phil) suffix(ac_bw3) outdir(`outdir')

// --- H84 (:1523): two-step GMM, HAC Tukey-Hanning, bw(3) ---
ivreg2 cinf (unem = l(1/3).unem), bw(3) gmm2s kernel(thann)
save_ivreg2_results, prefix(tsop_phil) suffix(gmm2s_thann) outdir(`outdir')

// --- H85 (:1526): two-step GMM, HAC QS, robust, orthog(l1.unem) ---
ivreg2 cinf (unem = l(1/3).unem), bw(3) gmm2s kernel(qs) robust orthog(l1.unem)
save_ivreg2_results, prefix(tsop_phil) suffix(gmm2s_qs_orthog) outdir(`outdir')

/*---------------------------------------------------------------------------
  Klein data (H72-H76): LIML-family estimators with L.profits / L.totinc
  tsset yr. Estimation columns exported for the R tests (klein is not
  bundled until Tranche 2).
---------------------------------------------------------------------------*/
capture use "../validation/data/klein.dta", clear
if _rc {
    webuse klein, clear
}
tsset yr

// Export estimation data for the R side. Full-precision export (%21.0g +
// datafmt): the default ~8-digit export does not round-trip float storage,
// which alone shifts the kclass(1.19) profits coefficient by 3e-6 relative
// (the validation F2 lesson, planning/13-validation-findings.md).
preserve
keep yr consump profits wagetot govt taxnetx year wagegovt capital1 totinc
format consump profits wagetot govt taxnetx year wagegovt capital1 totinc %21.0g
export delimited using "`outdir'/tsop_klein_data.csv", replace datafmt
restore

// --- H72 (:1462): LIML ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc), liml
save_ivreg2_results, prefix(tsop_klein) suffix(liml) outdir(`outdir')

// --- H73 (:1469): LIML + COVIV ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc), liml coviv
save_ivreg2_results, prefix(tsop_klein) suffix(liml_coviv) outdir(`outdir')

// --- H74 (:1473): CUE (LIML = CUE under iid demo) ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc), cue
save_ivreg2_results, prefix(tsop_klein) suffix(cue) outdir(`outdir')

// --- H75 (:1480): Fuller(1) ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc), fuller(1)
save_ivreg2_results, prefix(tsop_klein) suffix(fuller1) outdir(`outdir')

// --- H76 (:1487): k-class, k = 1.19 ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc), kclass(1.19)
save_ivreg2_results, prefix(tsop_klein) suffix(kclass119) outdir(`outdir')

/*---------------------------------------------------------------------------
  abdata (H88): panel differences as excluded instruments
  tsset id year. Estimation columns exported for the R tests.
---------------------------------------------------------------------------*/
capture use "../validation/data/abdata.dta", clear
if _rc {
    use http://fmwww.bc.edu/ec-p/data/macro/abdata.dta, clear
}
tsset id year

// Export estimation data for the R side (full precision; see klein note)
preserve
keep id year n w k ys
format n w k ys %21.0g
export delimited using "`outdir'/tsop_ab_data.csv", replace datafmt
restore

// --- H88 (:1541): two-step GMM, cluster(id), d. and d2. instruments ---
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), gmm2s cluster(id)
save_ivreg2_results, prefix(tsop_ab) suffix(gmm2s_cl) outdir(`outdir')

display "ts-operator fixture generation complete"
