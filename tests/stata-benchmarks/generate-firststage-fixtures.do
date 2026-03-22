* generate-firststage-fixtures.do
* Generate Stata fixtures for first-stage model objects (Ticket S1)
* Writes CSV files to fixtures/ for automated comparison in R tests.

clear all
set more off

* --- Load Card data ---
bcuse card, clear

* Helper program: save first-stage coefs + VCV to CSV files
* Must be defined AFTER bcuse (which calls clear all)
capture program drop save_fs_fixture
program define save_fs_fixture
    syntax , prefix(string) eqname(string)

    estimates restore `eqname'

    * --- Coefficients ---
    matrix b = e(b)
    local K = colsof(b)
    local names : colnames b
    file open fh using "fixtures/`prefix'_coef.csv", write replace
    file write fh "term,estimate" _n
    forval i = 1/`K' {
        local nm : word `i' of `names'
        file write fh "`nm'," %20.12f (b[1,`i']) _n
    }
    file close fh

    * --- VCV matrix ---
    matrix V = e(V)
    file open fh using "fixtures/`prefix'_vcov.csv", write replace
    * Header row
    file write fh "term"
    forval j = 1/`K' {
        local nm : word `j' of `names'
        file write fh ",vcov_`nm'"
    }
    file write fh _n
    * Data rows
    forval i = 1/`K' {
        local rn : word `i' of `names'
        file write fh "`rn'"
        forval j = 1/`K' {
            file write fh "," %20.12f (V[`i',`j'])
        }
        file write fh _n
    }
    file close fh

    * --- Scalars ---
    file open fh using "fixtures/`prefix'_scalars.csv", write replace
    file write fh "quantity,value" _n
    file write fh "rmse," %20.12f (e(rmse)) _n
    file write fh "df_r," %10.0f (e(df_r)) _n
    file write fh "N," %10.0f (e(N)) _n
    if !missing(e(N_clust)) {
        file write fh "N_clust," %10.0f (e(N_clust)) _n
    }
    file close fh

end


* =========================================================================
* 1. Overidentified, IID VCE
* =========================================================================
ivreg2 lwage exper expersq black south smsa (educ = nearc4 nearc2), savefirst
save_fs_fixture, prefix(fs_overid_iid) eqname(_ivreg2_educ)

* =========================================================================
* 2. Overidentified, Robust VCE
* =========================================================================
ivreg2 lwage exper expersq black south smsa (educ = nearc4 nearc2), savefirst robust
save_fs_fixture, prefix(fs_overid_robust) eqname(_ivreg2_educ)

* =========================================================================
* 3. Overidentified, Cluster VCE (age, M=11)
* =========================================================================
ivreg2 lwage exper expersq black south smsa (educ = nearc4 nearc2), savefirst cluster(age)
save_fs_fixture, prefix(fs_overid_cluster) eqname(_ivreg2_educ)

* =========================================================================
* 4. Just-identified, Robust VCE
* =========================================================================
ivreg2 lwage exper expersq black south smsa (educ = nearc4), savefirst robust
save_fs_fixture, prefix(fs_just_robust) eqname(_ivreg2_educ)

* =========================================================================
* 5. LIML (first-stage is same OLS regardless of method)
* =========================================================================
ivreg2 lwage exper expersq black south smsa (educ = nearc4 nearc2), savefirst liml
save_fs_fixture, prefix(fs_overid_liml) eqname(_ivreg2_educ)

* =========================================================================
* 6. Weighted (aweight), IID
* =========================================================================
gen aw = _n
ivreg2 lwage exper expersq black south smsa (educ = nearc4 nearc2) [aw=aw], savefirst
save_fs_fixture, prefix(fs_overid_aw_iid) eqname(_ivreg2_educ)

* =========================================================================
* 7. Weighted (aweight), Robust
* =========================================================================
ivreg2 lwage exper expersq black south smsa (educ = nearc4 nearc2) [aw=aw], savefirst robust
save_fs_fixture, prefix(fs_overid_aw_robust) eqname(_ivreg2_educ)
