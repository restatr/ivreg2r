/*===========================================================================
  generate-liml-fixtures.do
  -------------------------
  M-15 re-base: D5a option-variations on the klein H72-H76 LIML-family base and the abdata H88 base; canonical IID cells live in generate-ts-operator-fixtures.do as tsop_klein_*.

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-liml-fixtures.do

  Retired cells: Fuller(4) and kclass(0.5) were deleted as anti-pattern cells — they duplicated the option-variation shape of Fuller(1)/kclass(1.19) without adding new coverage.
  Just-identified LIML and kclass(1) cells are also retired: plain just-identified LIML short-circuits to 2SLS in R (a bit-identical identity test), and kclass(1) is a fixture-free identity with 2SLS itself; Stata parity for 2SLS is owned by the M-10 family, so neither needs a dedicated LIML-family fixture here.
  coviv x cluster, coviv x small, fuller x cluster, kclass x cluster, weighted x small, fuller x small, and kclass x small cells are retired compositionally: coviv only swaps the VCV bread and is already pinned at iid by the tsop H73 cell and at robust by klein_liml_coviv/hc1 here; fuller/kclass only change the k-class bread, and the k-class-bread x cluster-meat assembly is pinned by the ABDATA ab_liml cl/cl_small cells (cl_small is not a klein suffix); the small finite-sample correction is applied in the shared .vcov_from_omega with no method branch, so fuller x small and kclass x small are pinned by the klein_liml iid_small/hc1_small cells plus the method-uniformity identity test in test-liml.R.
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
  (Same as generate-ts-operator-fixtures.do: includes kclass/lambda/fuller and AR-LIML overid columns, plus a small-sample flag column since this generator's suffixes include iid_small/hc1_small/cl_small. First-stage matrix export omitted — not consumed by these tests.)
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

        // Small option flag
        gen double small = ("`suffix'" == "iid_small" | ///
                            "`suffix'" == "hc1_small" | ///
                            "`suffix'" == "cl_small")

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
  Klein data (H72-H76): LIML/Fuller/k-class option-variations
  tsset yr. R side uses the bundled klein dataset (no data export here).
---------------------------------------------------------------------------*/
capture use "../validation/data/klein.dta", clear
if _rc {
    webuse klein, clear
}
tsset yr

// Deterministic analytic weight (M-12/M-23 precedent) for the weighted cells
gen double awt = mod(yr, 5) + 1

// --- D5a option-variation on H72 (help.txt:1462): LIML + small ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc), liml small
save_ivreg2_results, prefix(klein_liml) suffix(iid_small) outdir(`outdir')

// --- D5a option-variation on H72 (help.txt:1462): LIML + robust ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc), liml robust
save_ivreg2_results, prefix(klein_liml) suffix(hc1) outdir(`outdir')

// --- D5a option-variation on H72 (help.txt:1462): LIML + robust + small ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc), liml robust small
save_ivreg2_results, prefix(klein_liml) suffix(hc1_small) outdir(`outdir')

// --- D5a option-variation on H73 (help.txt:1469): LIML + COVIV + robust ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc), liml coviv robust
save_ivreg2_results, prefix(klein_liml_coviv) suffix(hc1) outdir(`outdir')

// --- D5a option-variation on H75 (help.txt:1480): Fuller(1) + robust ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc), fuller(1) robust
save_ivreg2_results, prefix(klein_fuller1) suffix(hc1) outdir(`outdir')

// --- D5a option-variation on H76 (help.txt:1487): k-class 1.19 + robust ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc), kclass(1.19) robust
save_ivreg2_results, prefix(klein_kclass119) suffix(hc1) outdir(`outdir')

// --- D5a option-variation on H72 (help.txt:1462): weighted LIML, iid ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc) [aw=awt], liml
save_ivreg2_results, prefix(klein_liml_aw) suffix(iid) outdir(`outdir')

// --- D5a option-variation on H72 (help.txt:1462): weighted LIML, robust ---
ivreg2 consump L.profits (profits wagetot = govt taxnetx year ///
    wagegovt capital1 L.totinc) [aw=awt], liml robust
save_ivreg2_results, prefix(klein_liml_aw) suffix(hc1) outdir(`outdir')

/*---------------------------------------------------------------------------
  abdata (H88): LIML on the same panel-difference instrument set used by
  the gmm2s cell in generate-ts-operator-fixtures.do (help.txt:1541).
  tsset id year. R side uses the bundled abdata dataset (no data export here).
---------------------------------------------------------------------------*/
capture use "../validation/data/abdata.dta", clear
if _rc {
    use http://fmwww.bc.edu/ec-p/data/macro/abdata.dta, clear
}
tsset id year

// --- D5a option-variation on H88 (help.txt:1541): LIML + cluster(id) ---
ivreg2 n (w k ys = D.w D.k D.ys D2.w D2.k D2.ys), liml cluster(id)
save_ivreg2_results, prefix(ab_liml) suffix(cl) outdir(`outdir')

// --- D5a option-variation on H88 (help.txt:1541): LIML + cluster(id) + small ---
ivreg2 n (w k ys = D.w D.k D.ys D2.w D2.k D2.ys), liml cluster(id) small
save_ivreg2_results, prefix(ab_liml) suffix(cl_small) outdir(`outdir')

display "LIML fixture generation complete"
