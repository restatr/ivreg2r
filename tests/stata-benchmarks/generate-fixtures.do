/*===========================================================================
  generate-fixtures.do
  --------------------
  Generates CSV benchmark fixtures for ivreg2r R package testing.
  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)
  Each fixture produces:
    - *_data.csv          (dataset, for simulated data only)
    - *_coef_*.csv        (coefficients and SEs per VCE x small combo)
    - *_vcov_*.csv        (full VCV matrix per combo)
    - *_diagnostics_*.csv (test statistics per combo)
    - *_firststage_*.csv  (first-stage diagnostics per combo)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

// Pre-load Card data and save before defining programs.
// bcuse internally calls "clear all" which drops user-defined programs,
// so we must call it before any program definitions.
capture bcuse card, clear
if _rc != 0 {
    display as error "Could not load Card dataset via bcuse."
    display as error "Install bcuse (ssc install bcuse) and rerun."
    exit 601
}
save "`outdir'/_card_temp.dta", replace
export delimited using "`outdir'/card_data.csv", replace

/*---------------------------------------------------------------------------
  Helper program: extract all results from ivreg2 and save to CSV
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

        // Endogeneity / C-statistic (orthog() option → e(cstat))
        capture replace cstat = e(cstat)
        capture replace cstatp = e(cstatp)
        capture replace cstatdf = e(cstatdf)

        // Endogeneity test (endog() option → e(estat))
        gen double estat = .
        gen double estatp = .
        gen double estatdf = .
        capture replace estat = e(estat)
        capture replace estatp = e(estatp)
        capture replace estatdf = e(estatdf)

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
  Helper: run one ivreg2 specification, save results
  Uses globals $ivreg2_depvar $ivreg2_exog $ivreg2_endo $ivreg2_iv
  set by the caller before invoking run_all_vce_combos.
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
    // Globals are used because Stata's syntax parser cannot handle
    // nested parentheses inside string options.
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
  FIXTURE 1: card_just_id
  Dataset: Card (1995) — returns to education
  Model: lwage ~ exper expersq black south | educ | nearc4
  Purpose: Simplest real-data IV; exact identification (L1 = K1 = 1)
===========================================================================*/
display _newline(2) "=== FIXTURE 1: card_just_id ==="

// Reload Card data (pre-loaded before program definitions)
use "`outdir'/_card_temp.dta", clear

// Set model components via globals
global ivreg2_depvar "lwage"
global ivreg2_exog "exper expersq black south"
global ivreg2_endo "educ"
global ivreg2_iv "nearc4"

run_all_vce_combos, ///
    prefix(card_just_id) ///
    outdir(`outdir') ///
    endogopt(educ)


/*===========================================================================
  FIXTURE 2: card_overid
  Dataset: Card (1995)
  Model: lwage ~ exper expersq black south | educ | nearc2 nearc4
  Purpose: Overidentified (L1=2, K1=1); tests Sargan/Hansen J
===========================================================================*/
display _newline(2) "=== FIXTURE 2: card_overid ==="

// Reload Card data
use "`outdir'/_card_temp.dta", clear
global ivreg2_depvar "lwage"
global ivreg2_exog "exper expersq black south"
global ivreg2_endo "educ"
global ivreg2_iv "nearc2 nearc4"

run_all_vce_combos, ///
    prefix(card_overid) ///
    outdir(`outdir') ///
    endogopt(educ)


/*===========================================================================
  FIXTURE 3: sim_multi_endo
  Simulated data with 2 endogenous variables and 4 excluded IVs
  Purpose: Tests SW/AP, Shea partial R2, multi-endo eigenvalue paths
===========================================================================*/
display _newline(2) "=== FIXTURE 3: sim_multi_endo ==="

clear
set seed 20260214
set obs 1000

// Exogenous variables
gen x1 = rnormal()
gen x2 = rnormal()

// Excluded instruments (4 IVs for 2 endogenous)
gen z1 = rnormal()
gen z2 = rnormal()
gen z3 = rnormal()
gen z4 = rnormal()

// Structural errors (correlated with endogenous vars)
gen u = rnormal()
gen v1 = rnormal()
gen v2 = rnormal()

// Endogenous variables (correlated with u via v1, v2)
gen endo1 = 0.5*z1 + 0.3*z2 + 0.2*x1 + 0.8*v1 + 0.3*u
gen endo2 = 0.4*z2 + 0.3*z3 + 0.2*z4 + 0.1*x2 + 0.7*v2 + 0.2*u

// Outcome
gen y = 1.0 + 0.5*x1 - 0.3*x2 + 1.2*endo1 - 0.8*endo2 + u

// Save simulated data
export delimited using "`outdir'/sim_multi_endo_data.csv", replace

global ivreg2_depvar "y"
global ivreg2_exog "x1 x2"
global ivreg2_endo "endo1 endo2"
global ivreg2_iv "z1 z2 z3 z4"

run_all_vce_combos, ///
    prefix(sim_multi_endo) ///
    outdir(`outdir') ///
    endogopt(endo1 endo2)


/*===========================================================================
  FIXTURE 4: sim_no_constant
  Simulated data, model without intercept
  Purpose: Tests uncentered R-squared, no-constant code path
===========================================================================*/
display _newline(2) "=== FIXTURE 4: sim_no_constant ==="

clear
set seed 20260215
set obs 500

gen z1 = rnormal()
gen z2 = rnormal()
gen x1 = rnormal()
gen v = rnormal()
gen u = rnormal()

gen endo1 = 0.6*z1 + 0.4*z2 + 0.3*x1 + 0.7*v + 0.2*u
gen y = 2.0*x1 + 1.5*endo1 + u

// Save simulated data
export delimited using "`outdir'/sim_no_constant_data.csv", replace

global ivreg2_depvar "y"
global ivreg2_exog "x1"
global ivreg2_endo "endo1"
global ivreg2_iv "z1 z2"

// noconstant model
run_all_vce_combos, ///
    prefix(sim_no_constant) ///
    outdir(`outdir') ///
    endogopt(endo1) ///
    modelopts(noconstant)


/*===========================================================================
  FIXTURE 5: sim_cluster
  Simulated data with group structure (50 clusters)
  Purpose: Tests cluster VCV, cluster-robust diagnostics
===========================================================================*/
display _newline(2) "=== FIXTURE 5: sim_cluster ==="

clear
set seed 20260216
set obs 1000

// 50 clusters, 20 obs each
gen cluster_id = ceil(_n / 20)

// Cluster-level effects
gen alpha = rnormal() if mod(_n - 1, 20) == 0
bysort cluster_id: replace alpha = alpha[1]

// Individual-level variables
gen x1 = rnormal() + 0.3*alpha
gen z1 = rnormal() + 0.2*alpha
gen z2 = rnormal() + 0.1*alpha
gen v = rnormal()
gen u = rnormal() + 0.5*alpha

// Endogenous variable
gen endo1 = 0.5*z1 + 0.3*z2 + 0.2*x1 + 0.6*v + 0.3*u

// Outcome
gen y = 1.0 + 0.8*x1 + 1.5*endo1 + u

// Save simulated data
export delimited using "`outdir'/sim_cluster_data.csv", replace

global ivreg2_depvar "y"
global ivreg2_exog "x1"
global ivreg2_endo "endo1"
global ivreg2_iv "z1 z2"

// Run all combos including cluster
run_all_vce_combos, ///
    prefix(sim_cluster) ///
    outdir(`outdir') ///
    clustervar(cluster_id) ///
    endogopt(endo1)


/*===========================================================================
  FIXTURE 6: card_just_id_weighted
  Dataset: Card (1995) — returns to education
  Model: lwage ~ exper expersq black south | educ | nearc4 [aw=weight]
  Purpose: Weighted 2SLS with analytic weights; just-identified
===========================================================================*/
display _newline(2) "=== FIXTURE 6: card_just_id_weighted ==="

// Reload Card data
use "`outdir'/_card_temp.dta", clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [aw=weight], first endog(educ)
save_ivreg2_results, prefix(card_just_id_weighted) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [aw=weight], first small endog(educ)
save_ivreg2_results, prefix(card_just_id_weighted) suffix(iid_small) outdir(`outdir')

// --- HC1, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [aw=weight], first robust endog(educ)
save_ivreg2_results, prefix(card_just_id_weighted) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [aw=weight], first robust small endog(educ)
save_ivreg2_results, prefix(card_just_id_weighted) suffix(hc1_small) outdir(`outdir')

// --- Cluster on smsa66, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [aw=weight], first cluster(smsa66) endog(educ)
save_ivreg2_results, prefix(card_just_id_weighted) suffix(cl) outdir(`outdir')

// --- Cluster on smsa66, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [aw=weight], first cluster(smsa66) small endog(educ)
save_ivreg2_results, prefix(card_just_id_weighted) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  Clean up globals and done
===========================================================================*/
macro drop ivreg2_depvar ivreg2_exog ivreg2_endo ivreg2_iv
capture erase "`outdir'/_card_temp.dta"

display _newline(2) "=== All fixtures generated ==="
display "Output directory: `outdir'"
display "Run the R test suite to compare against these fixtures."
