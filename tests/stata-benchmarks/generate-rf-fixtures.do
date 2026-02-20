/*===========================================================================
  generate-rf-fixtures.do
  -----------------------
  Generates CSV benchmark fixtures for reduced-form regression tests.
  Run in Stata with ivreg2 installed (ssc install ivreg2).

  Uses ivreg2's saverf and savesfirst options to extract:
    - Coefficients and SEs from the reduced-form y ~ Z regression
    - F-stat of excluded instruments
    - RMSE

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-rf-fixtures.do
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
save "`outdir'/_card_temp_rf.dta", replace


/*---------------------------------------------------------------------------
  Helper program: extract reduced-form results from stored RF estimate

  After ivreg2 with saverf, the RF equation is stored as an estimate.
  We restore it, extract b, V, and compute F-test of excluded IVs.
---------------------------------------------------------------------------*/
capture program drop save_rf_results
program define save_rf_results
    syntax, prefix(string) suffix(string) outdir(string) ///
           rfeq(string) [excludedivs(string)]

    // Restore the RF estimate
    estimates restore `rfeq'

    local N = e(N)

    // Extract coefficients and SEs
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

        // RMSE and F-stat (from e() results of stored RF)
        gen double rmse = e(rmse)
        gen double N_obs = e(N)

        // F-test of excluded instruments
        gen double rf_f = .
        gen double rf_f_p = .
        gen double rf_f_df1 = .
        gen double rf_f_df2 = .

        // Test excluded instruments = 0
        if "`excludedivs'" != "" {
            capture test `excludedivs'
            if _rc == 0 {
                replace rf_f = r(F) in 1
                replace rf_f_p = r(p) in 1
                replace rf_f_df1 = r(df) in 1
                replace rf_f_df2 = r(df_r) in 1
            }
        }

        export delimited using ///
            "`outdir'/`prefix'_rf_`suffix'.csv", replace
        restore
    }
end


/*---------------------------------------------------------------------------
  Helper program: extract system (savesfirst) results
---------------------------------------------------------------------------*/
capture program drop save_system_results
program define save_system_results
    syntax, prefix(string) suffix(string) outdir(string) ///
           sfirsteq(string)

    // Restore the system estimate
    estimates restore `sfirsteq'

    matrix b = e(b)
    matrix V = e(V)
    local names : colnames b
    local eqnames : coleq b
    local ncols = colsof(b)

    quietly {
        preserve
        clear
        set obs `ncols'
        gen str64 eq = ""
        gen str32 term = ""
        gen double estimate = .
        gen double std_error = .

        forvalues i = 1/`ncols' {
            local nm : word `i' of `names'
            local eq : word `i' of `eqnames'
            replace eq = "`eq'" in `i'
            replace term = "`nm'" in `i'
            replace estimate = b[1, `i'] in `i'
            replace std_error = sqrt(V[`i', `i']) in `i'
        }

        gen double rmse = e(rmse)
        gen double N_obs = e(N)

        export delimited using ///
            "`outdir'/`prefix'_system_`suffix'.csv", replace
        restore
    }
end


/*---------------------------------------------------------------------------
  Helper: run one ivreg2 + saverf + savesfirst, save both RF and system
  Uses globals $ivreg2_cmd for the full command (before comma).
  This avoids problems with [aw=weight] placement.
---------------------------------------------------------------------------*/
capture program drop run_rf_save
program define run_rf_save
    syntax, prefix(string) suffix(string) outdir(string) ///
           rfeq(string) sfirsteq(string) excludedivs(string) ///
           [opts(string)]

    $ivreg2_cmd, saverf savesfirst `opts'

    save_rf_results, prefix(`prefix') suffix(`suffix') outdir(`outdir') ///
        rfeq(`rfeq') excludedivs(`excludedivs')
    save_system_results, prefix(`prefix') suffix(`suffix') outdir(`outdir') ///
        sfirsteq(`sfirsteq')
end


/*===========================================================================
  FIXTURE 1: card_just_id
  Dataset: Card (1995) — returns to education
  Model: lwage ~ exper expersq black south | educ | nearc4
  Just-identified (1 endogenous, 1 excluded IV)
===========================================================================*/
display _newline(2) "=== FIXTURE 1: card_just_id ==="

use "`outdir'/_card_temp_rf.dta", clear

global ivreg2_cmd "ivreg2 lwage exper expersq black south (educ = nearc4)"
local rfeq "_ivreg2_lwage"
local sfirsteq "_ivreg2_sfirst_lwage"
local excludedivs "nearc4"
local pre "card_just_id"

// IID
run_rf_save, prefix(`pre') suffix(iid) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs')

// IID small
run_rf_save, prefix(`pre') suffix(iid_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(small)

// HC1
run_rf_save, prefix(`pre') suffix(hc1) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(robust)

// HC1 small
run_rf_save, prefix(`pre') suffix(hc1_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(robust small)

// Cluster
run_rf_save, prefix(`pre') suffix(cl) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(cluster(smsa))

// Cluster small
run_rf_save, prefix(`pre') suffix(cl_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(cluster(smsa) small)


/*===========================================================================
  FIXTURE 2: card_overid
  Dataset: Card (1995) — overidentified model
  Model: lwage ~ exper expersq black south | educ | nearc2 nearc4
===========================================================================*/
display _newline(2) "=== FIXTURE 2: card_overid ==="

use "`outdir'/_card_temp_rf.dta", clear

global ivreg2_cmd "ivreg2 lwage exper expersq black south (educ = nearc2 nearc4)"
local rfeq "_ivreg2_lwage"
local sfirsteq "_ivreg2_sfirst_lwage"
local excludedivs "nearc2 nearc4"
local pre "card_overid"

// IID
run_rf_save, prefix(`pre') suffix(iid) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs')

// IID small
run_rf_save, prefix(`pre') suffix(iid_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(small)

// HC1
run_rf_save, prefix(`pre') suffix(hc1) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(robust)

// HC1 small
run_rf_save, prefix(`pre') suffix(hc1_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(robust small)

// Cluster
run_rf_save, prefix(`pre') suffix(cl) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(cluster(smsa))

// Cluster small
run_rf_save, prefix(`pre') suffix(cl_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(cluster(smsa) small)


/*===========================================================================
  FIXTURE 3: card_just_id_weighted
  Dataset: Card (1995) — weighted model
  Model: lwage ~ exper expersq black south | educ | nearc4 [aw=weight]
  Weights go before the comma, so we set them in the global command.
===========================================================================*/
display _newline(2) "=== FIXTURE 3: card_just_id_weighted ==="

use "`outdir'/_card_temp_rf.dta", clear

global ivreg2_cmd "ivreg2 lwage exper expersq black south (educ = nearc4) [aw=weight]"
local rfeq "_ivreg2_lwage"
local sfirsteq "_ivreg2_sfirst_lwage"
local excludedivs "nearc4"
local pre "card_just_id_weighted"

// IID
run_rf_save, prefix(`pre') suffix(iid) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs')

// IID small
run_rf_save, prefix(`pre') suffix(iid_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(small)

// HC1
run_rf_save, prefix(`pre') suffix(hc1) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(robust)

// HC1 small
run_rf_save, prefix(`pre') suffix(hc1_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(robust small)

// Cluster
run_rf_save, prefix(`pre') suffix(cl) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(cluster(smsa))

// Cluster small
run_rf_save, prefix(`pre') suffix(cl_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(cluster(smsa) small)


/*===========================================================================
  FIXTURE 4: sim_multi_endo
  Simulated data with 2 endogenous variables and 4 excluded IVs
  Purpose: Tests system mode with K1 > 1
===========================================================================*/
display _newline(2) "=== FIXTURE 4: sim_multi_endo ==="

import delimited "`outdir'/sim_multi_endo_data.csv", clear

global ivreg2_cmd "ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4)"
local rfeq "_ivreg2_y"
local sfirsteq "_ivreg2_sfirst_y"
local excludedivs "z1 z2 z3 z4"
local pre "sim_multi_endo"

// IID
run_rf_save, prefix(`pre') suffix(iid) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs')

// IID small
run_rf_save, prefix(`pre') suffix(iid_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(small)

// HC1
run_rf_save, prefix(`pre') suffix(hc1) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(robust)

// HC1 small
run_rf_save, prefix(`pre') suffix(hc1_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(robust small)


/*===========================================================================
  FIXTURE 5: sim_cluster
  Simulated cluster data
  Purpose: Tests cluster-robust RF
===========================================================================*/
display _newline(2) "=== FIXTURE 5: sim_cluster ==="

import delimited "`outdir'/sim_cluster_data.csv", clear

global ivreg2_cmd "ivreg2 y x1 (endo1 = z1 z2)"
local rfeq "_ivreg2_y"
local sfirsteq "_ivreg2_sfirst_y"
local excludedivs "z1 z2"
local pre "sim_cluster"

// Cluster
run_rf_save, prefix(`pre') suffix(cl) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(cluster(cluster_id))

// Cluster small
run_rf_save, prefix(`pre') suffix(cl_small) outdir(`outdir') ///
    rfeq(`rfeq') sfirsteq(`sfirsteq') excludedivs(`excludedivs') ///
    opts(cluster(cluster_id) small)


/*===========================================================================
  Cleanup
===========================================================================*/
capture erase "`outdir'/_card_temp_rf.dta"

display _newline(2) "=== All RF fixtures generated successfully ==="
