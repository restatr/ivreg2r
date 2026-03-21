* ==========================================================================
* Generate Stata fixtures for R2 (Additional Stored Results)
* ==========================================================================
* Quantities: e(yy), e(yyc), e(rankxx), e(rankzz), e(condxx), e(condzz),
*             e(ll), e(ccev), e(cdev)

clear all
set more off

* --- Load Card data ---
bcuse card, clear nodesc

cap gen lwage = log(wage)

* --- Create output directory ---
cap mkdir fixtures

* ==========================================================================
* Config 1: 2SLS just-identified (IID)
* ==========================================================================
ivreg2 lwage exper expersq black south (educ = nearc4)

mat ccev_vec = e(ccev)
mat cdev_vec = e(cdev)

file open fh using "fixtures/stored_results_justid.csv", write replace
file write fh "quantity,value" _n
file write fh "yy," %20.12f (e(yy)) _n
file write fh "yyc," %20.12f (e(yyc)) _n
file write fh "rankxx," %5.0f (e(rankxx)) _n
file write fh "rankzz," %5.0f (e(rankzz)) _n
file write fh "condxx," %20.12f (e(condxx)) _n
file write fh "condzz," %20.12f (e(condzz)) _n
file write fh "ll," %20.12f (e(ll)) _n
file write fh "N," %10.0f (e(N)) _n
file write fh "rss," %20.12f (e(rss)) _n
local nccev = colsof(ccev_vec)
forval i = 1/`nccev' {
    file write fh "ccev`i'," %20.12f (ccev_vec[1,`i']) _n
    file write fh "cdev`i'," %20.12f (cdev_vec[1,`i']) _n
}
file close fh

* ==========================================================================
* Config 2: 2SLS overidentified (robust)
* ==========================================================================
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2), robust

mat ccev_vec = e(ccev)
mat cdev_vec = e(cdev)

file open fh using "fixtures/stored_results_overid_robust.csv", write replace
file write fh "quantity,value" _n
file write fh "yy," %20.12f (e(yy)) _n
file write fh "yyc," %20.12f (e(yyc)) _n
file write fh "rankxx," %5.0f (e(rankxx)) _n
file write fh "rankzz," %5.0f (e(rankzz)) _n
file write fh "condxx," %20.12f (e(condxx)) _n
file write fh "condzz," %20.12f (e(condzz)) _n
file write fh "ll," %20.12f (e(ll)) _n
file write fh "N," %10.0f (e(N)) _n
file write fh "rss," %20.12f (e(rss)) _n
local nccev = colsof(ccev_vec)
forval i = 1/`nccev' {
    file write fh "ccev`i'," %20.12f (ccev_vec[1,`i']) _n
    file write fh "cdev`i'," %20.12f (cdev_vec[1,`i']) _n
}
file close fh

* ==========================================================================
* Config 3: 2SLS overidentified with clustering
* ==========================================================================
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2), cluster(smsa)

mat ccev_vec = e(ccev)
mat cdev_vec = e(cdev)

file open fh using "fixtures/stored_results_overid_cluster.csv", write replace
file write fh "quantity,value" _n
file write fh "yy," %20.12f (e(yy)) _n
file write fh "yyc," %20.12f (e(yyc)) _n
file write fh "rankxx," %5.0f (e(rankxx)) _n
file write fh "rankzz," %5.0f (e(rankzz)) _n
file write fh "condxx," %20.12f (e(condxx)) _n
file write fh "condzz," %20.12f (e(condzz)) _n
file write fh "ll," %20.12f (e(ll)) _n
file write fh "N," %10.0f (e(N)) _n
file write fh "rss," %20.12f (e(rss)) _n
local nccev = colsof(ccev_vec)
forval i = 1/`nccev' {
    file write fh "ccev`i'," %20.12f (ccev_vec[1,`i']) _n
    file write fh "cdev`i'," %20.12f (cdev_vec[1,`i']) _n
}
file close fh

* ==========================================================================
* Config 4: Aweight
* ==========================================================================
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [aw=age]

mat ccev_vec = e(ccev)
mat cdev_vec = e(cdev)

file open fh using "fixtures/stored_results_aweight.csv", write replace
file write fh "quantity,value" _n
file write fh "yy," %20.12f (e(yy)) _n
file write fh "yyc," %20.12f (e(yyc)) _n
file write fh "rankxx," %5.0f (e(rankxx)) _n
file write fh "rankzz," %5.0f (e(rankzz)) _n
file write fh "condxx," %20.12f (e(condxx)) _n
file write fh "condzz," %20.12f (e(condzz)) _n
file write fh "ll," %20.12f (e(ll)) _n
file write fh "N," %10.0f (e(N)) _n
file write fh "rss," %20.12f (e(rss)) _n
local nccev = colsof(ccev_vec)
forval i = 1/`nccev' {
    file write fh "ccev`i'," %20.12f (ccev_vec[1,`i']) _n
    file write fh "cdev`i'," %20.12f (cdev_vec[1,`i']) _n
}
file close fh

* ==========================================================================
* Config 5: LIML
* ==========================================================================
ivreg2 lwage exper expersq black south (educ = nearc4 nearc2), liml

mat ccev_vec = e(ccev)
mat cdev_vec = e(cdev)

file open fh using "fixtures/stored_results_liml.csv", write replace
file write fh "quantity,value" _n
file write fh "yy," %20.12f (e(yy)) _n
file write fh "yyc," %20.12f (e(yyc)) _n
file write fh "rankxx," %5.0f (e(rankxx)) _n
file write fh "rankzz," %5.0f (e(rankzz)) _n
file write fh "condxx," %20.12f (e(condxx)) _n
file write fh "condzz," %20.12f (e(condzz)) _n
file write fh "ll," %20.12f (e(ll)) _n
file write fh "N," %10.0f (e(N)) _n
file write fh "rss," %20.12f (e(rss)) _n
local nccev = colsof(ccev_vec)
forval i = 1/`nccev' {
    file write fh "ccev`i'," %20.12f (ccev_vec[1,`i']) _n
    file write fh "cdev`i'," %20.12f (cdev_vec[1,`i']) _n
}
file close fh

* ==========================================================================
* Config 6: Fweight
* ==========================================================================
bcuse card, clear nodesc
cap gen lwage = log(wage)

ivreg2 lwage exper expersq black south (educ = nearc4 nearc2) [fw=age]

mat ccev_vec = e(ccev)
mat cdev_vec = e(cdev)

file open fh using "fixtures/stored_results_fweight.csv", write replace
file write fh "quantity,value" _n
file write fh "yy," %20.12f (e(yy)) _n
file write fh "yyc," %20.12f (e(yyc)) _n
file write fh "rankxx," %5.0f (e(rankxx)) _n
file write fh "rankzz," %5.0f (e(rankzz)) _n
file write fh "condxx," %20.12f (e(condxx)) _n
file write fh "condzz," %20.12f (e(condzz)) _n
file write fh "ll," %20.12f (e(ll)) _n
file write fh "N," %10.0f (e(N)) _n
file write fh "rss," %20.12f (e(rss)) _n
local nccev = colsof(ccev_vec)
forval i = 1/`nccev' {
    file write fh "ccev`i'," %20.12f (ccev_vec[1,`i']) _n
    file write fh "cdev`i'," %20.12f (cdev_vec[1,`i']) _n
}
file close fh

* ==========================================================================
* Config 7: OLS (no instruments — rankzz/condzz should be missing)
* ==========================================================================
bcuse card, clear nodesc
cap gen lwage = log(wage)

ivreg2 lwage exper expersq black south educ

file open fh using "fixtures/stored_results_ols.csv", write replace
file write fh "quantity,value" _n
file write fh "yy," %20.12f (e(yy)) _n
file write fh "yyc," %20.12f (e(yyc)) _n
file write fh "rankxx," %5.0f (e(rankxx)) _n
* rankzz may equal rankxx for OLS in Stata
cap file write fh "rankzz," %5.0f (e(rankzz)) _n
file write fh "condxx," %20.12f (e(condxx)) _n
* condzz may equal condxx for OLS in Stata
cap file write fh "condzz," %20.12f (e(condzz)) _n
file write fh "ll," %20.12f (e(ll)) _n
file write fh "N," %10.0f (e(N)) _n
file write fh "rss," %20.12f (e(rss)) _n
file close fh

di "=== Fixture generation complete ==="
