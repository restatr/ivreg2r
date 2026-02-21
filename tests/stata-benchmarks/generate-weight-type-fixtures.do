/*===========================================================================
  generate-weight-type-fixtures.do
  --------------------------------
  Generates CSV benchmark fixtures for frequency weights (fweight),
  probability weights (pweight), and weighted overidentified aweight tests.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be pkg/):
    cd /path/to/ivreg2r/pkg
    do tests/stata-benchmarks/generate-weight-type-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

// Pre-load Card data and save before defining programs.
// bcuse calls "clear all" which drops user-defined programs.
capture bcuse card, clear
if _rc != 0 {
    display as error "Could not load Card dataset via bcuse."
    display as error "Install bcuse (ssc install bcuse) and rerun."
    exit 601
}
save "`outdir'/_card_temp.dta", replace

/*---------------------------------------------------------------------------
  Helper program: extract all results from ivreg2 and save to CSV
  (Identical to the one in generate-fixtures.do)
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

        // Endogeneity test
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
                            "`suffix'" == "cl_small" | ///
                            "`suffix'" == "hc0_small")

        export delimited using ///
            "`outdir'/`prefix'_diagnostics_`suffix'.csv", replace
        restore
    }

    // --- First-stage diagnostics ---
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
  FIXTURE SET 1: Frequency weights (fweight) — just identified
  Model: lwage ~ exper expersq black south | educ | nearc4  [fw=weight]
===========================================================================*/
use "`outdir'/_card_temp.dta", clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [fw=weight], first
save_ivreg2_results, prefix(card_fweight_just_id) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [fw=weight], first small
save_ivreg2_results, prefix(card_fweight_just_id) suffix(iid_small) outdir(`outdir')

// --- HC1, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [fw=weight], first robust
save_ivreg2_results, prefix(card_fweight_just_id) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [fw=weight], first robust small
save_ivreg2_results, prefix(card_fweight_just_id) suffix(hc1_small) outdir(`outdir')

// --- Cluster, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [fw=weight], first cluster(smsa)
save_ivreg2_results, prefix(card_fweight_just_id) suffix(cl) outdir(`outdir')

// --- Cluster, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [fw=weight], first cluster(smsa) small
save_ivreg2_results, prefix(card_fweight_just_id) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 2: Frequency weights (fweight) — overidentified
  Model: lwage ~ exper expersq black south | educ | nearc4 nearc2  [fw=weight]
===========================================================================*/
use "`outdir'/_card_temp.dta", clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=weight], first
save_ivreg2_results, prefix(card_fweight_overid) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=weight], first small
save_ivreg2_results, prefix(card_fweight_overid) suffix(iid_small) outdir(`outdir')

// --- HC1, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=weight], first robust
save_ivreg2_results, prefix(card_fweight_overid) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=weight], first robust small
save_ivreg2_results, prefix(card_fweight_overid) suffix(hc1_small) outdir(`outdir')

// --- Cluster, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=weight], first cluster(smsa)
save_ivreg2_results, prefix(card_fweight_overid) suffix(cl) outdir(`outdir')

// --- Cluster, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=weight], first cluster(smsa) small
save_ivreg2_results, prefix(card_fweight_overid) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 3: Probability weights (pweight) — just identified
  Model: lwage ~ exper expersq black south | educ | nearc4  [pw=weight]
  Note: pweight forces robust VCE, so no iid configs
===========================================================================*/
use "`outdir'/_card_temp.dta", clear

// --- HC0 (Stata's default for pweight), small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [pw=weight], first
save_ivreg2_results, prefix(card_pweight_just_id) suffix(hc0) outdir(`outdir')

// --- HC0, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [pw=weight], first small
save_ivreg2_results, prefix(card_pweight_just_id) suffix(hc0_small) outdir(`outdir')

// --- Cluster, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [pw=weight], first cluster(smsa)
save_ivreg2_results, prefix(card_pweight_just_id) suffix(cl) outdir(`outdir')

// --- Cluster, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4) [pw=weight], first cluster(smsa) small
save_ivreg2_results, prefix(card_pweight_just_id) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 4: Probability weights (pweight) — overidentified
  Model: lwage ~ exper expersq black south | educ | nearc4 nearc2  [pw=weight]
===========================================================================*/
use "`outdir'/_card_temp.dta", clear

// --- HC0, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [pw=weight], first
save_ivreg2_results, prefix(card_pweight_overid) suffix(hc0) outdir(`outdir')

// --- HC0, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [pw=weight], first small
save_ivreg2_results, prefix(card_pweight_overid) suffix(hc0_small) outdir(`outdir')

// --- Cluster, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [pw=weight], first cluster(smsa)
save_ivreg2_results, prefix(card_pweight_overid) suffix(cl) outdir(`outdir')

// --- Cluster, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [pw=weight], first cluster(smsa) small
save_ivreg2_results, prefix(card_pweight_overid) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  FIXTURE SET 5: Analytic weights (aweight) — overidentified
  Fills Tier 1 gap: no weighted overid fixtures exist yet.
  Model: lwage ~ exper expersq black south | educ | nearc4 nearc2  [aw=weight]
===========================================================================*/
use "`outdir'/_card_temp.dta", clear

// --- IID, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [aw=weight], first
save_ivreg2_results, prefix(card_aweight_overid) suffix(iid) outdir(`outdir')

// --- IID, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [aw=weight], first small
save_ivreg2_results, prefix(card_aweight_overid) suffix(iid_small) outdir(`outdir')

// --- HC1, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [aw=weight], first robust
save_ivreg2_results, prefix(card_aweight_overid) suffix(hc1) outdir(`outdir')

// --- HC1, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [aw=weight], first robust small
save_ivreg2_results, prefix(card_aweight_overid) suffix(hc1_small) outdir(`outdir')

// --- Cluster, small=FALSE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [aw=weight], first cluster(smsa)
save_ivreg2_results, prefix(card_aweight_overid) suffix(cl) outdir(`outdir')

// --- Cluster, small=TRUE ---
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [aw=weight], first cluster(smsa) small
save_ivreg2_results, prefix(card_aweight_overid) suffix(cl_small) outdir(`outdir')


/*===========================================================================
  Cleanup
===========================================================================*/
capture erase "`outdir'/_card_temp.dta"

display ""
display "=========================================="
display "Weight-type fixtures generated successfully"
display "=========================================="
