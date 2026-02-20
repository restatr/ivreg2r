/*===========================================================================
  generate-twoway-fixtures.do
  ---------------------------
  Generates CSV benchmark fixtures for two-way cluster-robust VCE testing.
  Uses Cameron-Gelbach-Miller (2006) formula via Stata's ivreg2 cluster2().

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-twoway-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"


/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results with two-way cluster metadata
---------------------------------------------------------------------------*/
capture program drop save_twoway_results
program define save_twoway_results
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

    // --- Diagnostics (with N_clust1 and N_clust2) ---
    quietly {
        preserve
        clear
        set obs 1

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

        // Initialize all conditional columns
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
        gen double N_clust1 = .
        gen double N_clust2 = .
        gen double sstat = .
        gen double sstatp = .
        gen double sstatdf = .

        // Overidentification
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

        // Endogeneity / C-statistic
        capture replace cstat = e(cstat)
        capture replace cstatp = e(cstatp)
        capture replace cstatdf = e(cstatdf)

        // Endogeneity test (endog() option)
        gen double estat = .
        gen double estatp = .
        gen double estatdf = .
        capture replace estat = e(estat)
        capture replace estatp = e(estatp)
        capture replace estatdf = e(estatdf)

        // Stock-Wright S statistic
        capture replace sstat = e(sstat)
        capture replace sstatp = e(sstatp)
        capture replace sstatdf = e(sstatdf)

        // Cluster counts
        capture replace N_clust = e(N_clust)
        capture replace N_clust1 = e(N_clust1)
        capture replace N_clust2 = e(N_clust2)

        // Small option flag
        local is_small = (strpos("`suffix'", "small") > 0)
        gen double small = `is_small'

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


/*===========================================================================
  Simulated DGP: N=500 balanced panel, 25 firms x 20 obs (10 years x 2)
  Cluster-correlated errors so two-way clustering differs from one-way.
===========================================================================*/
display _newline(2) "=== Generating two-way cluster data ==="

clear
set seed 20260220
set obs 500

// Panel structure: 25 firms x 20 obs each
gen firm_id = ceil(_n / 20)

// 10 years, cycling within firm
gen year_id = mod(_n - 1, 10) + 1

// Firm-level effects
gen alpha_firm = rnormal() if mod(_n - 1, 20) == 0
bysort firm_id: replace alpha_firm = alpha_firm[1]

// Year-level effects
sort year_id
gen alpha_year = rnormal() if _n == 1 | year_id != year_id[_n-1]
bysort year_id: replace alpha_year = alpha_year[1]

// Individual-level variables
gen x1 = rnormal() + 0.3 * alpha_firm + 0.2 * alpha_year
gen z1 = rnormal() + 0.2 * alpha_firm + 0.1 * alpha_year
gen z2 = rnormal() + 0.1 * alpha_firm + 0.15 * alpha_year
gen v = rnormal()
gen u = rnormal() + 0.5 * alpha_firm + 0.4 * alpha_year

// Endogenous variable
gen endo1 = 0.5 * z1 + 0.3 * z2 + 0.2 * x1 + 0.6 * v + 0.3 * u

// Outcome
gen y = 1.0 + 0.8 * x1 + 1.5 * endo1 + u

// Weight variable (strictly positive)
gen wt = runiform() * 3 + 0.5

// Save simulated data
sort firm_id year_id
export delimited using "`outdir'/sim_twoway_data.csv", replace


/*===========================================================================
  FIXTURE 1: IV overid, two-way cluster (small=FALSE)
===========================================================================*/
display _newline(2) "=== Fixture: cl2 ==="
ivreg2 y x1 (endo1 = z1 z2), cluster(firm_id year_id) first endog(endo1)
save_twoway_results, prefix(sim_twoway) suffix(cl2) outdir(`outdir')


/*===========================================================================
  FIXTURE 2: IV overid, two-way cluster (small=TRUE)
===========================================================================*/
display _newline(2) "=== Fixture: cl2_small ==="
ivreg2 y x1 (endo1 = z1 z2), cluster(firm_id year_id) first small endog(endo1)
save_twoway_results, prefix(sim_twoway) suffix(cl2_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 3: OLS, two-way cluster (small=FALSE)
===========================================================================*/
display _newline(2) "=== Fixture: cl2_ols ==="
ivreg2 y x1 endo1, cluster(firm_id year_id)
save_twoway_results, prefix(sim_twoway) suffix(cl2_ols) outdir(`outdir')


/*===========================================================================
  FIXTURE 4: OLS, two-way cluster (small=TRUE)
===========================================================================*/
display _newline(2) "=== Fixture: cl2_ols_small ==="
ivreg2 y x1 endo1, cluster(firm_id year_id) small
save_twoway_results, prefix(sim_twoway) suffix(cl2_ols_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 5: IV overid weighted, two-way cluster (small=FALSE)
===========================================================================*/
display _newline(2) "=== Fixture: cl2_wt ==="
ivreg2 y x1 (endo1 = z1 z2) [aw=wt], cluster(firm_id year_id) first endog(endo1)
save_twoway_results, prefix(sim_twoway) suffix(cl2_wt) outdir(`outdir')


/*===========================================================================
  FIXTURE 6: IV overid weighted, two-way cluster (small=TRUE)
===========================================================================*/
display _newline(2) "=== Fixture: cl2_wt_small ==="
ivreg2 y x1 (endo1 = z1 z2) [aw=wt], cluster(firm_id year_id) first small endog(endo1)
save_twoway_results, prefix(sim_twoway) suffix(cl2_wt_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 7: IV overid dofminus(1) sdofminus(1), two-way cluster (small=FALSE)
===========================================================================*/
display _newline(2) "=== Fixture: cl2_dof ==="
ivreg2 y x1 (endo1 = z1 z2), cluster(firm_id year_id) first endog(endo1) ///
    dofminus(1) sdofminus(1)
save_twoway_results, prefix(sim_twoway) suffix(cl2_dof) outdir(`outdir')


/*===========================================================================
  FIXTURE 8: IV overid dofminus(1) sdofminus(1), two-way cluster (small=TRUE)
===========================================================================*/
display _newline(2) "=== Fixture: cl2_dof_small ==="
ivreg2 y x1 (endo1 = z1 z2), cluster(firm_id year_id) first small endog(endo1) ///
    dofminus(1) sdofminus(1)
save_twoway_results, prefix(sim_twoway) suffix(cl2_dof_small) outdir(`outdir')


/*===========================================================================
  Clean up and done
===========================================================================*/
display _newline(2) "=== All two-way cluster fixtures generated ==="
display "Output directory: `outdir'"
