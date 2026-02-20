/*===========================================================================
  generate-liml-fixtures.do
  -------------------------
  Generates CSV benchmark fixtures for LIML, Fuller, and k-class estimation
  in the ivreg2r R package.

  Reuses the helper programs from generate-fixtures.do.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-liml-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

// Pre-load Card data before defining programs (bcuse calls clear all)
capture bcuse card, clear
if _rc != 0 {
    display as error "Could not load Card dataset via bcuse."
    display as error "Install bcuse (ssc install bcuse) and rerun."
    exit 601
}
save "`outdir'/_card_liml_temp.dta", replace

/*---------------------------------------------------------------------------
  Helper program: extract all results from ivreg2 and save to CSV
  (Same as generate-fixtures.do but adds kclass/lambda/fuller columns)
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
        forvalues i = 1/`dim' {
            local nm : word `i' of `names'
            local cnm = subinstr("`nm'", ".", "_", .)
            rename V`i' vcov_`cnm'
        }
        export delimited using "`outdir'/`prefix'_vcov_`suffix'.csv", replace
        restore
    }

    // --- Diagnostics (extended with LIML metadata) ---
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
        gen double sigmasq = e(rmse)^2
        gen double F_stat = e(F)
        gen double F_p = e(Fp)
        gen double F_df1 = e(Fdf1)
        gen double F_df2 = e(Fdf2)

        // Initialize all conditional columns to missing
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
        gen double sstat = .
        gen double sstatp = .
        gen double sstatdf = .
        capture replace sstat = e(sstat)
        capture replace sstatp = e(sstatp)
        capture replace sstatdf = e(sstatdf)

        // Cluster info
        capture replace N_clust = e(N_clust)

        // Small option flag
        gen double small = ("`suffix'" == "iid_small" | ///
                            "`suffix'" == "hc1_small" | ///
                            "`suffix'" == "cl_small")

        // --- LIML/k-class metadata ---
        gen double kclass = .
        gen double lambda = .
        gen double fuller = .
        capture replace kclass = e(kclass)
        capture replace lambda = e(lambda)
        capture replace fuller = e(fuller)

        // --- AR LIML overidentification (H3) ---
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
  FIXTURE 1: card_liml_overid
  Dataset: Card (1995)
  Model: lwage ~ exper expersq black south (educ = nearc2 nearc4), liml
  Purpose: LIML overidentified, IID + small
===========================================================================*/
display _newline(2) "=== FIXTURE 1: card_liml_overid ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml first
save_ivreg2_results, prefix(card_liml_overid) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml first small
save_ivreg2_results, prefix(card_liml_overid) suffix(iid_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 2: card_liml_justid
  Dataset: Card (1995)
  Model: lwage ~ exper expersq black south (educ = nearc4), liml
  Purpose: Exactly-identified LIML → lambda=1, should equal 2SLS
===========================================================================*/
display _newline(2) "=== FIXTURE 2: card_liml_justid ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4), liml first
save_ivreg2_results, prefix(card_liml_justid) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4), liml first small
save_ivreg2_results, prefix(card_liml_justid) suffix(iid_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 3: card_fuller1_overid
  Dataset: Card (1995)
  Model: lwage ~ exper expersq black south (educ = nearc2 nearc4), fuller(1)
  Purpose: Fuller(1) modification
===========================================================================*/
display _newline(2) "=== FIXTURE 3: card_fuller1_overid ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), fuller(1) first
save_ivreg2_results, prefix(card_fuller1_overid) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), fuller(1) first small
save_ivreg2_results, prefix(card_fuller1_overid) suffix(iid_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 4: card_fuller4_overid
  Dataset: Card (1995)
  Model: lwage ~ exper expersq black south (educ = nearc2 nearc4), fuller(4)
  Purpose: Larger Fuller alpha
===========================================================================*/
display _newline(2) "=== FIXTURE 4: card_fuller4_overid ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), fuller(4) first
save_ivreg2_results, prefix(card_fuller4_overid) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), fuller(4) first small
save_ivreg2_results, prefix(card_fuller4_overid) suffix(iid_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 5: card_kclass_half
  Dataset: Card (1995)
  Model: lwage ~ exper expersq black south (educ = nearc2 nearc4), kclass(0.5)
  Purpose: Arbitrary k value
===========================================================================*/
display _newline(2) "=== FIXTURE 5: card_kclass_half ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), kclass(0.5) first
save_ivreg2_results, prefix(card_kclass_half) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), kclass(0.5) first small
save_ivreg2_results, prefix(card_kclass_half) suffix(iid_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 6: card_kclass_1
  Dataset: Card (1995)
  Model: lwage ~ exper expersq black south (educ = nearc2 nearc4), kclass(1)
  Purpose: k=1 must equal 2SLS
===========================================================================*/
display _newline(2) "=== FIXTURE 6: card_kclass_1 ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), kclass(1) first
save_ivreg2_results, prefix(card_kclass_1) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), kclass(1) first small
save_ivreg2_results, prefix(card_kclass_1) suffix(iid_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 7: card_liml_weighted
  Dataset: Card (1995)
  Model: lwage ~ exper expersq black south (educ = nearc2 nearc4) [aw=weight], liml
  Purpose: Weighted LIML
===========================================================================*/
display _newline(2) "=== FIXTURE 7: card_liml_weighted ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4) [aw=weight], liml first
save_ivreg2_results, prefix(card_liml_weighted) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4) [aw=weight], liml first small
save_ivreg2_results, prefix(card_liml_weighted) suffix(iid_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 8: sim_multi_endo_liml
  Simulated data with 2 endogenous variables and 4 excluded IVs
  Purpose: Multi-endogenous LIML (K1>1 eigenvalue)
===========================================================================*/
display _newline(2) "=== FIXTURE 8: sim_multi_endo_liml ==="

import delimited using "`outdir'/sim_multi_endo_data.csv", clear

// --- IID, small=FALSE ---
ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4), liml first
save_ivreg2_results, prefix(sim_multi_endo_liml) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4), liml first small
save_ivreg2_results, prefix(sim_multi_endo_liml) suffix(iid_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 9: sim_multi_endo_fuller1
  Multi-endogenous + Fuller(1)
===========================================================================*/
display _newline(2) "=== FIXTURE 9: sim_multi_endo_fuller1 ==="

import delimited using "`outdir'/sim_multi_endo_data.csv", clear

// --- IID, small=FALSE ---
ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4), fuller(1) first
save_ivreg2_results, prefix(sim_multi_endo_fuller1) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4), fuller(1) first small
save_ivreg2_results, prefix(sim_multi_endo_fuller1) suffix(iid_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 10: card_liml_overid — robust/cluster
  Purpose: LIML overid with HC and cluster VCE
===========================================================================*/
display _newline(2) "=== FIXTURE 10: card_liml_overid robust/cluster ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- HC1 (Stata robust), small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml robust first
save_ivreg2_results, prefix(card_liml_overid) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml robust small first
save_ivreg2_results, prefix(card_liml_overid) suffix(hc1_small) outdir(`outdir')

// --- Cluster on smsa66, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml cluster(smsa66) first
save_ivreg2_results, prefix(card_liml_overid) suffix(cl) outdir(`outdir')

// --- Cluster on smsa66, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml cluster(smsa66) small first
save_ivreg2_results, prefix(card_liml_overid) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 11: card_liml_overid_coviv
  Purpose: LIML overid with COVIV (coviv option), all VCE types
===========================================================================*/
display _newline(2) "=== FIXTURE 11: card_liml_overid_coviv ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- IID, small=FALSE, coviv ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml coviv first
save_ivreg2_results, prefix(card_liml_overid_coviv) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE, coviv ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml coviv small first
save_ivreg2_results, prefix(card_liml_overid_coviv) suffix(iid_small) outdir(`outdir')

// --- HC1, small=FALSE, coviv ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml robust coviv first
save_ivreg2_results, prefix(card_liml_overid_coviv) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE, coviv ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml robust small coviv first
save_ivreg2_results, prefix(card_liml_overid_coviv) suffix(hc1_small) outdir(`outdir')

// --- Cluster, small=FALSE, coviv ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml cluster(smsa66) coviv first
save_ivreg2_results, prefix(card_liml_overid_coviv) suffix(cl) outdir(`outdir')

// --- Cluster, small=TRUE, coviv ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), liml cluster(smsa66) small coviv first
save_ivreg2_results, prefix(card_liml_overid_coviv) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 12: card_liml_justid — robust/cluster
  Purpose: Just-identified LIML (k=1), verifies LIML=2SLS under robust/cluster
===========================================================================*/
display _newline(2) "=== FIXTURE 12: card_liml_justid robust/cluster ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- HC1, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4), liml robust first
save_ivreg2_results, prefix(card_liml_justid) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4), liml robust small first
save_ivreg2_results, prefix(card_liml_justid) suffix(hc1_small) outdir(`outdir')

// --- Cluster, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4), liml cluster(smsa66) first
save_ivreg2_results, prefix(card_liml_justid) suffix(cl) outdir(`outdir')

// --- Cluster, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4), liml cluster(smsa66) small first
save_ivreg2_results, prefix(card_liml_justid) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 13: card_fuller1_overid — robust/cluster
  Purpose: Fuller(1) + robust/cluster VCE
===========================================================================*/
display _newline(2) "=== FIXTURE 13: card_fuller1_overid robust/cluster ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- HC1, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), fuller(1) robust first
save_ivreg2_results, prefix(card_fuller1_overid) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), fuller(1) robust small first
save_ivreg2_results, prefix(card_fuller1_overid) suffix(hc1_small) outdir(`outdir')

// --- Cluster, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), fuller(1) cluster(smsa66) first
save_ivreg2_results, prefix(card_fuller1_overid) suffix(cl) outdir(`outdir')

// --- Cluster, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), fuller(1) cluster(smsa66) small first
save_ivreg2_results, prefix(card_fuller1_overid) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 14: card_kclass_half — robust
  Purpose: kclass(0.5) + robust VCE
===========================================================================*/
display _newline(2) "=== FIXTURE 14: card_kclass_half robust ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- HC1, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), kclass(0.5) robust first
save_ivreg2_results, prefix(card_kclass_half) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4), kclass(0.5) robust small first
save_ivreg2_results, prefix(card_kclass_half) suffix(hc1_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 15: card_liml_weighted — robust
  Purpose: Weighted LIML + robust VCE
===========================================================================*/
display _newline(2) "=== FIXTURE 15: card_liml_weighted robust ==="

use "`outdir'/_card_liml_temp.dta", clear

// --- HC1, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4) [aw=weight], liml robust first
save_ivreg2_results, prefix(card_liml_weighted) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc2 nearc4) [aw=weight], liml robust small first
save_ivreg2_results, prefix(card_liml_weighted) suffix(hc1_small) outdir(`outdir')


/*===========================================================================
  FIXTURE 16: sim_multi_endo_liml — robust
  Purpose: Multi-endogenous LIML + robust VCE
===========================================================================*/
display _newline(2) "=== FIXTURE 16: sim_multi_endo_liml robust ==="

import delimited using "`outdir'/sim_multi_endo_data.csv", clear

// --- HC1, small=FALSE ---
ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4), liml robust first
save_ivreg2_results, prefix(sim_multi_endo_liml) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE ---
ivreg2 y x1 x2 (endo1 endo2 = z1 z2 z3 z4), liml robust small first
save_ivreg2_results, prefix(sim_multi_endo_liml) suffix(hc1_small) outdir(`outdir')


/*===========================================================================
  Clean up
===========================================================================*/
capture erase "`outdir'/_card_liml_temp.dta"

display _newline(2) "=== All LIML fixtures generated ==="
display "Output directory: `outdir'"
