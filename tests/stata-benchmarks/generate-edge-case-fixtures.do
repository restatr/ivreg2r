/*===========================================================================
  generate-edge-case-fixtures.do
  ------------------------------
  Generates CSV benchmark fixtures for edge-case models in ivreg2r testing.
  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Fixture: sim_intercept_only
    Intercept-only exogenous model: y ~ 1 | endo1 | z1 + z2
    No exogenous regressors beyond the constant.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (from the stata-benchmarks/ directory):
    cd /path/to/ivreg2r/pkg/tests/stata-benchmarks
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b generate-edge-case-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory (relative to CWD = stata-benchmarks/)
local outdir "fixtures"
capture mkdir "`outdir'"


/*---------------------------------------------------------------------------
  Helper program: extract all results from ivreg2 and save to CSV
  (Copied from generate-fixtures.do — bcuse calls clear all, so programs
   must be defined after any data loading that uses bcuse.)
---------------------------------------------------------------------------*/
capture program drop save_ivreg2_results
program define save_ivreg2_results
    syntax, prefix(string) suffix(string) outdir(string)

    // Number of observations
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
        // Rename columns
        forvalues i = 1/`dim' {
            local nm : word `i' of `names'
            // Clean name for variable naming (replace . with _)
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

        // Model stats (always returned by ivreg2)
        gen double N = `N'
        gen double K = `K'
        gen double L = `L'
        gen double r2 = e(r2)
        gen double r2_a = e(r2_a)
        gen double r2u = e(r2u)
        gen double r2c = e(r2c)
        gen double rss = e(rss)
        gen double rmse = e(rmse)
        gen double sigmasq = e(rmse)^2
        gen double F_stat = e(F)
        gen double F_p = e(Fp)
        gen double F_df1 = e(Fdf1)
        gen double F_df2 = e(Fdf2)

        // Initialize all conditional columns to missing (.) so every CSV
        // has the same schema regardless of model specification.
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
        gen double N_clust = .

        // Overidentification (conditionally returned)
        capture replace sargan = e(sargan)
        capture replace sarganp = e(sarganp)
        capture replace sargandf = e(sargandf)
        capture replace j = e(j)
        capture replace jp = e(jp)
        capture replace jdf = e(jdf)

        // Underidentification
        capture replace idstat = e(idstat)
        capture replace idp = e(idp)
        capture replace iddf = e(iddf)

        // Weak identification
        capture replace cdf = e(cdf)
        capture replace widstat = e(widstat)

        // Anderson-Rubin
        capture replace arf = e(arf)
        capture replace arfp = e(arfp)
        capture replace archi2 = e(archi2)
        capture replace archi2p = e(archi2p)
        capture replace ardf = e(ardf)
        capture replace ardf_r = e(ardf_r)

        // Endogeneity / C-statistic (orthog() option -> e(cstat))
        capture replace cstat = e(cstat)
        capture replace cstatp = e(cstatp)
        capture replace cstatdf = e(cstatdf)

        // Endogeneity test (endog() option -> e(estat))
        gen double estat = .
        gen double estatp = .
        gen double estatdf = .
        capture replace estat = e(estat)
        capture replace estatp = e(estatp)
        capture replace estatdf = e(estatdf)

        // Stock-Wright S statistic
        gen double sstat = .
        gen double sstatp = .
        gen double sstatdf = .
        capture replace sstat = e(sstat)
        capture replace sstatp = e(sstatp)
        capture replace sstatdf = e(sstatdf)

        // Cluster info (only when clustered)
        capture replace N_clust = e(N_clust)

        // Small option flag
        gen double small = ("`suffix'" == "iid_small" | ///
                            "`suffix'" == "hc1_small" | ///
                            "`suffix'" == "cl_small")

        export delimited using ///
            "`outdir'/`prefix'_diagnostics_`suffix'.csv", replace
        restore
    }

    // --- First-stage diagnostics (if IV model) ---
    capture confirm matrix e(first)
    if _rc == 0 {
        quietly {
            preserve
            clear

            // e(first) matrix: rows are statistics, columns are endogenous vars
            matrix F = e(first)
            local frows = rowsof(F)
            local fcols = colsof(F)
            local fnames : colnames F
            local rnames : rownames F

            set obs `frows'
            gen str32 statistic = ""
            forvalues i = 1/`frows' {
                local rn : word `i' of `rnames'
                replace statistic = "`rn'" in `i'
            }

            forvalues j = 1/`fcols' {
                local cn : word `j' of `fnames'
                local cnclean = subinstr("`cn'", ".", "_", .)
                gen double `cnclean' = .
                forvalues i = 1/`frows' {
                    replace `cnclean' = F[`i', `j'] in `i'
                }
            }

            export delimited using ///
                "`outdir'/`prefix'_firststage_`suffix'.csv", replace
            restore
        }
    }
end


/*---------------------------------------------------------------------------
  Helper: run one ivreg2 specification across all VCE combos
  Uses globals $ivreg2_depvar $ivreg2_exog $ivreg2_endo $ivreg2_iv
---------------------------------------------------------------------------*/
capture program drop run_all_vce_combos
program define run_all_vce_combos
    syntax, prefix(string) outdir(string) [clustervar(string) ///
            endogopt(string) modelopts(string)]

    // Build endog() option string
    local endogstr ""
    if "`endogopt'" != "" {
        local endogstr "endog(`endogopt')"
    }

    // Build the ivreg2 command from globals.
    local cmd "ivreg2 $ivreg2_depvar $ivreg2_exog"
    if "$ivreg2_endo" != "" {
        local cmd "`cmd' ($ivreg2_endo = $ivreg2_iv)"
    }

    // --- IID, small=FALSE ---
    `cmd', first `modelopts' `endogstr'
    save_ivreg2_results, prefix(`prefix') suffix(iid) outdir(`outdir')

    // --- IID, small=TRUE ---
    `cmd', first `modelopts' small `endogstr'
    save_ivreg2_results, prefix(`prefix') suffix(iid_small) outdir(`outdir')

    // --- HC1, small=FALSE ---
    `cmd', first `modelopts' robust `endogstr'
    save_ivreg2_results, prefix(`prefix') suffix(hc1) outdir(`outdir')

    // --- HC1, small=TRUE ---
    `cmd', first `modelopts' robust small `endogstr'
    save_ivreg2_results, prefix(`prefix') suffix(hc1_small) outdir(`outdir')

    // --- Cluster (if cluster variable specified) ---
    if "`clustervar'" != "" {
        // Cluster, small=FALSE
        `cmd', first `modelopts' cluster(`clustervar') `endogstr'
        save_ivreg2_results, prefix(`prefix') suffix(cl) outdir(`outdir')

        // Cluster, small=TRUE
        `cmd', first `modelopts' cluster(`clustervar') small `endogstr'
        save_ivreg2_results, prefix(`prefix') suffix(cl_small) outdir(`outdir')
    }
end


/*===========================================================================
  FIXTURE: sim_intercept_only
  Intercept-only exogenous model: y ~ 1 | endo1 | z1 + z2
  K2=1 (intercept only), K1=1, L1=2, overid_df=1
===========================================================================*/
display _newline(2) "=== FIXTURE: sim_intercept_only ==="

clear
set seed 20260221
set obs 500

// Excluded instruments
gen z1 = rnormal()
gen z2 = rnormal()

// Structural errors
gen v = rnormal()
gen u = rnormal()

// Endogenous variable (correlated with u via v)
gen endo1 = 0.6*z1 + 0.4*z2 + 0.7*v + 0.3*u

// Outcome: y = 3.0 + 1.5*endo1 + u (no exogenous regressors beyond intercept)
gen y = 3.0 + 1.5*endo1 + u

// Save simulated data
export delimited using "`outdir'/sim_intercept_only_data.csv", replace

global ivreg2_depvar "y"
global ivreg2_exog ""
global ivreg2_endo "endo1"
global ivreg2_iv "z1 z2"

run_all_vce_combos, ///
    prefix(sim_intercept_only) ///
    outdir(`outdir') ///
    endogopt(endo1)


/*===========================================================================
  Clean up globals and done
===========================================================================*/
macro drop ivreg2_depvar ivreg2_exog ivreg2_endo ivreg2_iv

display _newline(2) "=== Edge-case fixtures generated ==="
display "Output directory: `outdir'"
