/*===========================================================================
  generate-helpfile-fixtures.do
  ------------------------------
  Generates CSV benchmark fixtures for the in-package help-file example
  suite (planning/24-execution-roadmap.md, decision D1): reproduce every
  `ivreg2` invocation in the Examples section of the Stata ivreg2 help file
  (planning/help-file-specs.md, H01-H105) against the bundled ivreg2r
  datasets, so test-helpfile-examples.R can assert Stata parity on the
  package's own worked examples rather than only on synthetic/Card fixtures.

  Coverage: griliches76, mroz, phillips, abdata, grunfeld, nlswork.
  klein is intentionally OMITTED here -- the 5 klein commands (H72-H76) are
  already fixture-verified in generate-ts-operator-fixtures.do (time-series
  operator examples), and HOLS commands (H34/H35/H48/H52) are covered by
  generate-hols-fixtures.do. See the coverage ledger at the top of
  test-helpfile-examples.R for the full H01-H105 accounting.

  Data sources: local cache in ../validation/data/ (built by
  validation/validate-helpfile.qmd); falls back to the BC.edu / webuse
  sources cited in the help file if the cache is absent.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-helpfile-fixtures.do
===========================================================================*/

clear all
set more off
version 14

local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"

/*---------------------------------------------------------------------------
  Helper program: extract ivreg2 results and save to CSV.
  Cloned from generate-ts-operator-fixtures.do's fullest variant: coef/vcov/
  diagnostics CSVs including kclass/lambda/fuller, cstat, sstat, arubin,
  N_clust (needed across the griliches CUE/orthog, klein-style LIML probes
  we do NOT run here, mroz endog/redundant/AR, and abdata/grunfeld/nlswork
  cluster+kernel examples).
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
        gen double estat = .
        gen double estatp = .
        gen double estatdf = .
        gen double redstat = .
        gen double redp = .
        gen double reddf = .
        gen double N_clust = .
        gen double N_clust1 = .
        gen double N_clust2 = .

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
        capture replace estat = e(estat)
        capture replace estatp = e(estatp)
        capture replace estatdf = e(estatdf)
        capture replace redstat = e(redstat)
        capture replace redp = e(redp)
        capture replace reddf = e(reddf)
        capture replace N_clust = e(N_clust)
        capture replace N_clust1 = e(N_clust1)
        capture replace N_clust2 = e(N_clust2)

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


/*===========================================================================
  Griliches76 (H01-H29): help.txt:1128-1263
===========================================================================*/

capture use "../validation/data/griliches76.dta", clear
if _rc {
    use http://fmwww.bc.edu/ec-p/data/hayashi/griliches76.dta, clear
}
quietly xi i.year

// --- H03/H23 (help.txt:1138/1223): baseline 2SLS ---
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt)
save_ivreg2_results, prefix(hf_gril) suffix(H03) outdir(`outdir')

// --- H04 (help.txt:1141): small + ffirst ---
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), small ffirst
save_ivreg2_results, prefix(hf_gril) suffix(H04) outdir(`outdir')

// --- H06 (help.txt:1154): GMM2S, robust ---
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), gmm2s robust
save_ivreg2_results, prefix(hf_gril) suffix(H06) outdir(`outdir')

// --- H07 (help.txt:1162): robust IV (also base model for the S/W-matrix
//     demos below) ---
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), robust
save_ivreg2_results, prefix(hf_gril) suffix(H07) outdir(`outdir')

// --- H08-H11 (help.txt:1165-1174): S-matrix demo. Rerun the H07 fit,
//     predict residuals, build S = (1/N) Z'diag(uhat^2) Z, refit GMM2S
//     with smatrix(S) ---
quietly {
    predict double uhat if e(sample), resid
    mat accum S = `e(insts)' [iw=uhat^2]
    mat S = 1/`e(N)' * S
}
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), gmm2s robust smatrix(S)
save_ivreg2_results, prefix(hf_gril) suffix(H11) outdir(`outdir')

// --- H12-H13 (help.txt:1178-1182): W-matrix demo. W = invsym(S) from
//     the S-matrix demo above ---
quietly mat W = invsym(S)
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), gmm2s robust wmatrix(W)
save_ivreg2_results, prefix(hf_gril) suffix(H13) outdir(`outdir')

// --- H14-H21 (help.txt:1189-1210): Ahn (1997) J = Wald equivalence demo.
//     Base GMM2S fit with only iq endogenous (med kww age excluded);
//     S0 = e(S); refit three times with the SAME instrument set but
//     different included/excluded partitions, reusing smatrix(S0). ---
capture use "../validation/data/griliches76.dta", clear
if _rc {
    use http://fmwww.bc.edu/ec-p/data/hayashi/griliches76.dta, clear
}

// --- H14 (help.txt:1189): base fit, instruments = _cons med kww age ---
ivreg2 lw (iq=med kww age), gmm2s
save_ivreg2_results, prefix(hf_gril) suffix(H14) outdir(`outdir')
quietly mat S0 = e(S)

// --- H16 (help.txt:1195): kww excluded, med age included ---
qui ivreg2 lw (iq=kww) med age, gmm2s smatrix(S0)
save_ivreg2_results, prefix(hf_gril) suffix(H16) outdir(`outdir')

// --- H18 (help.txt:1201): med excluded, kww age included ---
qui ivreg2 lw (iq=med) kww age, gmm2s smatrix(S0)
save_ivreg2_results, prefix(hf_gril) suffix(H18) outdir(`outdir')

// --- H20 (help.txt:1207): age excluded, med kww included ---
qui ivreg2 lw (iq=age) med kww, gmm2s smatrix(S0)
save_ivreg2_results, prefix(hf_gril) suffix(H20) outdir(`outdir')

// --- H22 (help.txt:1217): CUE, robust. Stata's optimizer converges to a
//     KNOWN DEGENERATE BASIN on this 13-regressor model (R^2 = -54.4,
//     J = 40.1) -- a documented CUE ratio-objective pathology, not a bug
//     in either implementation. See planning/14-cue-optimizer-research.md
//     and STATUS.md. Stata's answer is recorded here DELIBERATELY so the
//     R test can assert the divergence explicitly (a staleness tripwire),
//     not because we expect R to reproduce it. ---
capture use "../validation/data/griliches76.dta", clear
if _rc {
    use http://fmwww.bc.edu/ec-p/data/hayashi/griliches76.dta, clear
}
quietly xi i.year
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), cue robust
save_ivreg2_results, prefix(hf_gril) suffix(H22) outdir(`outdir')
// H22 is a divergence tripwire: the R test compares COEFFICIENTS only
// (asserting they do NOT match). Erase the vcov/diagnostics files the
// shared saver wrote so no orphaned fixtures ship (planning/18 rule).
erase "`outdir'/hf_gril_vcov_H22.csv"
erase "`outdir'/hf_gril_diagnostics_H22.csv"

// --- H25 (help.txt:1235): GMM2S, orthog on an included regressor (s) ---
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), gmm2s orthog(s)
save_ivreg2_results, prefix(hf_gril) suffix(H25) outdir(`outdir')

// --- H26 (help.txt:1242): GMM2S, orthog on two excluded IVs (age mrt) ---
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), gmm2s orthog(age mrt)
save_ivreg2_results, prefix(hf_gril) suffix(H26) outdir(`outdir')

// --- H27 (help.txt:1250): cluster(year); note only 3 excluded IVs (no mrt) ---
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), cluster(year)
save_ivreg2_results, prefix(hf_gril) suffix(H27) outdir(`outdir')

// --- H28 (help.txt:1253): cluster(year) + partial(_I*) (FWL demo) ---
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), cluster(year) partial(_I*)
save_ivreg2_results, prefix(hf_gril) suffix(H28) outdir(`outdir')

// --- H29 (help.txt:1261): DOCUMENTED HELP-FILE BUG. The help file claims
//     partialling out the year dummies makes GMM2S feasible with
//     #clusters < #instruments. It does not: after partialling, 5 included
//     + 3 excluded = 8 moment conditions remain against M = 7 year
//     clusters, so the cluster-robust S matrix is singular. Stata itself
//     rejects this command (r(506): "optimal GMM weighting matrix not
//     unique"). Encoded as an assert, not a fixture. ---
capture noisily ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), ///
    cluster(year) partial(_I*) gmm2s
local rc29 = _rc
assert `rc29' == 506
display "H29 (help.txt:1261): confirmed Stata rejects this command, rc=`rc29'"


/*===========================================================================
  Mroz (H30-H69): help.txt:1268-1449
===========================================================================*/

capture use "../validation/data/mroz.dta", clear
if _rc {
    use http://fmwww.bc.edu/ec-p/data/wooldridge/mroz.dta, clear
}

// --- H31 (help.txt:1274): baseline 2SLS ---
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6)
save_ivreg2_results, prefix(hf_mroz) suffix(H31) outdir(`outdir')

// --- H33 (help.txt:1283): endog(educ) ---
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6), endog(educ)
save_ivreg2_results, prefix(hf_mroz) suffix(H33) outdir(`outdir')

// --- H41/H46 (help.txt:1325/1350): robust IV. H41 adds savefirst (Stata
//     bookkeeping only, does not change the main-equation estimates); H46
//     is the identical spec without savefirst. One fixture covers both. ---
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6), robust savefirst
save_ivreg2_results, prefix(hf_mroz) suffix(H41) outdir(`outdir')

// --- H50 (help.txt:1368): robust + redundant(age) ---
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6), robust redundant(age)
save_ivreg2_results, prefix(hf_mroz) suffix(H50) outdir(`outdir')

// --- H54 (help.txt:1385): robust + ffirst + saverf (AR / Stock-Wright demo;
//     ffirst/saverf are Stata bookkeeping options) ---
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6), robust ffirst saverf
save_ivreg2_results, prefix(hf_mroz) suffix(H54) outdir(`outdir')

// --- H61 (help.txt:1416): reduced-form OLS, robust (AR chi-sq demo) ---
ivreg2 lwage exper expersq age kidslt6 kidsge6, robust
save_ivreg2_results, prefix(hf_mroz) suffix(H61) outdir(`outdir')

// --- H66 (help.txt:1440): bare OLS (b0 demo, coefficient source) ---
qui ivreg2 lwage exper expersq
save_ivreg2_results, prefix(hf_mroz) suffix(H66) outdir(`outdir')

// --- H64-H68 (help.txt:1434-1446): b0 demo. Build b = (0, OLS coefs) in
//     the IV model's e(sample), exactly as help.txt lines 1434-1443 /
//     validate-helpfile.qmd #mroz-b0 do. ---
quietly {
    ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6)
    gen byte mroz_iv_sample = e(sample)
    ivreg2 lwage exper expersq if mroz_iv_sample
    mat b_ols = e(b)
    mat b0 = (0)
    mat colnames b0 = educ
    mat b0 = b0, b_ols
}
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6) if mroz_iv_sample, robust b0(b0)
save_ivreg2_results, prefix(hf_mroz) suffix(H68) outdir(`outdir')


/*===========================================================================
  Phillips (H77-H85): help.txt:1494-1527
  Klein (H70-H76) is covered by generate-ts-operator-fixtures.do -- OMITTED.
  Only the pure-OLS AC/HAC examples (H79-H81) are fixtured here; the IV +
  ts-operator examples (H83-H85) are also covered by
  generate-ts-operator-fixtures.do.
===========================================================================*/

capture use "../validation/data/phillips.dta", clear
if _rc {
    use http://fmwww.bc.edu/ec-p/data/wooldridge/phillips.dta, clear
}
quietly tsset year, yearly

// --- H79 (help.txt:1501): AC, default Bartlett kernel, bw(3) ---
ivreg2 cinf unem, bw(3)
save_ivreg2_results, prefix(hf_phil) suffix(H79) outdir(`outdir')

// --- H80 (help.txt:1504): AC, Quadratic Spectral kernel, auto-bandwidth ---
ivreg2 cinf unem, kernel(qs) bw(auto)
save_ivreg2_results, prefix(hf_phil) suffix(H80) outdir(`outdir')

// --- H81 (help.txt:1511): HAC, Bartlett bw(3), robust, small ---
ivreg2 cinf unem, bw(3) kernel(bartlett) robust small
save_ivreg2_results, prefix(hf_phil) suffix(H81) outdir(`outdir')


/*===========================================================================
  Abdata (H86-H92): help.txt:1533-1561
  H88 (gmm2s + d./d2. instruments) is covered by
  generate-ts-operator-fixtures.do -- OMITTED here.
===========================================================================*/

capture use "../validation/data/abdata.dta", clear
if _rc {
    use http://fmwww.bc.edu/ec-p/data/macro/abdata.dta, clear
}
quietly tsset id year

// --- H89 (help.txt:1548): Kiefer SEs ---
ivreg2 n w k, kiefer
save_ivreg2_results, prefix(hf_ab) suffix(H89) outdir(`outdir')

// --- H90 (help.txt:1551): DOCUMENTED HELP-FILE BUG. bw(9) exceeds the
//     T=9 timespan (max lag = T-1 = 8); Stata rejects it even though
//     kiefer internally uses bw = T_span. bw(8) is numerically identical
//     (truncated kernel gives weight 1 to every available lag). ---
capture noisily ivreg2 n w k, bw(9) kernel(tru)
local rc90 = _rc
assert `rc90' == 198
display "H90 (help.txt:1551): confirmed Stata rejects bw(9), rc=`rc90'"
ivreg2 n w k, bw(8) kernel(tru)
save_ivreg2_results, prefix(hf_ab) suffix(H90) outdir(`outdir')

// --- H91 (help.txt:1558): one-way cluster on id ---
ivreg2 n w k, cluster(id)
save_ivreg2_results, prefix(hf_ab) suffix(H91) outdir(`outdir')

// --- H92 (help.txt:1561): SAME help-file bug as H90 (bw(9) -> bw(8)),
//     HAC variant (+ robust) ---
capture noisily ivreg2 n w k, bw(9) kernel(tru) robust
local rc92 = _rc
assert `rc92' == 198
display "H92 (help.txt:1561): confirmed Stata rejects bw(9), rc=`rc92'"
ivreg2 n w k, bw(8) kernel(tru) robust
save_ivreg2_results, prefix(hf_ab) suffix(H92) outdir(`outdir')


/*===========================================================================
  Grunfeld (H93-H100): help.txt:1569-1601
===========================================================================*/

capture use "../validation/data/grunfeld.dta", clear
if _rc {
    webuse grunfeld, clear
}
quietly tsset

// --- H95 (help.txt:1576): AC, Truncated kernel, bw(1) ---
ivreg2 invest mvalue kstock, bw(1) kernel(tru)
save_ivreg2_results, prefix(hf_grun) suffix(H95) outdir(`outdir')

// --- H96 (help.txt:1582): HAC, robust, Truncated kernel, bw(1) ---
ivreg2 invest mvalue kstock, robust bw(1) kernel(tru)
save_ivreg2_results, prefix(hf_grun) suffix(H96) outdir(`outdir')

// --- H97 (help.txt:1588): cluster(year) + robust + Truncated kernel bw(1) ---
ivreg2 invest mvalue kstock, robust cluster(year) bw(1) kernel(tru)
save_ivreg2_results, prefix(hf_grun) suffix(H97) outdir(`outdir')

// --- H98 (help.txt:1595): Driscoll-Kraay, bw(2), small ---
ivreg2 invest mvalue kstock, dkraay(2) small
save_ivreg2_results, prefix(hf_grun) suffix(H98) outdir(`outdir')

// --- H99 (help.txt:1598): DK equivalence via cluster(year) + Bartlett
//     bw(2), small ---
ivreg2 invest mvalue kstock, cluster(year) bw(2) small
save_ivreg2_results, prefix(hf_grun) suffix(H99) outdir(`outdir')


/*===========================================================================
  nlswork (H101-H105): help.txt:1604-1638
===========================================================================*/

capture use "../validation/data/nlswork.dta", clear
if _rc {
    webuse nlswork, clear
}
quietly tsset

// --- H103 (help.txt:1618): one-way cluster on idcode ---
ivreg2 ln_w grade age ttl_exp tenure, cluster(idcode)
save_ivreg2_results, prefix(hf_nls) suffix(H103) outdir(`outdir')

// --- H104 (help.txt:1628): two-way cluster on idcode x year ---
ivreg2 ln_w grade age ttl_exp tenure, cluster(idcode year)
save_ivreg2_results, prefix(hf_nls) suffix(H104) outdir(`outdir')

// --- H105 (help.txt:1637): two-way cluster + Truncated kernel bw(2) ---
ivreg2 ln_w grade age ttl_exp tenure, cluster(idcode year) bw(2) kernel(tru)
save_ivreg2_results, prefix(hf_nls) suffix(H105) outdir(`outdir')


/*===========================================================================
  Stored scalar results (family M-28, absorbed into the hf suite 2026-07-04
  per planning/22-spec-matrix.md): e(yy), e(yyc), e(rankxx), e(rankzz),
  e(condxx), e(condzz), e(ll), e(ccev), e(cdev) on the canonical mroz
  bases H31 (2SLS, help.txt:1274) and H61 (OLS via ivreg2, help.txt:1416).
  These are bookkeeping outputs with no worked example in the help file
  (decision D5a); they do not vary with the VCE type, so two cells
  (IV overidentified + OLS) suffice -- estimator/weight variation is
  compositionally covered by each family's rss/r2/N parity fixtures.
  Consumed by test-stored-results.R.
===========================================================================*/

capture use "../validation/data/mroz.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/wooldridge/mroz.dta, clear
    if _rc {
        display as error "Could not load mroz dataset (no local cache, no network)."
        exit 601
    }
}

// --- H31 (help.txt:1274): baseline 2SLS -- stored scalars + ccev/cdev ---
ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6)

mat ccev_vec = e(ccev)
mat cdev_vec = e(cdev)

file open fh using "`outdir'/hf_mroz_stored_H31.csv", write replace
file write fh "quantity,value" _n
file write fh "yy," %21.0g (e(yy)) _n
file write fh "yyc," %21.0g (e(yyc)) _n
file write fh "rankxx," %21.0g (e(rankxx)) _n
file write fh "rankzz," %21.0g (e(rankzz)) _n
file write fh "condxx," %21.0g (e(condxx)) _n
file write fh "condzz," %21.0g (e(condzz)) _n
file write fh "ll," %21.0g (e(ll)) _n
file write fh "N," %21.0g (e(N)) _n
file write fh "rss," %21.0g (e(rss)) _n
local nccev = colsof(ccev_vec)
forval i = 1/`nccev' {
    file write fh "ccev`i'," %21.0g (ccev_vec[1,`i']) _n
    file write fh "cdev`i'," %21.0g (cdev_vec[1,`i']) _n
}
file close fh

// --- H61 (help.txt:1416): reduced-form OLS via ivreg2 -- stored scalars.
//     The robust option is part of the verbatim H61 command but does not
//     affect any of these scalars. No ccev/cdev block (OLS has none). ---
ivreg2 lwage exper expersq age kidslt6 kidsge6, robust

file open fh using "`outdir'/hf_mroz_stored_H61.csv", write replace
file write fh "quantity,value" _n
file write fh "yy," %21.0g (e(yy)) _n
file write fh "yyc," %21.0g (e(yyc)) _n
file write fh "rankxx," %21.0g (e(rankxx)) _n
* Stata sets rankzz = rankxx for OLS; an unconditional write fails loudly
* if that ever changes.
file write fh "rankzz," %21.0g (e(rankzz)) _n
file write fh "condxx," %21.0g (e(condxx)) _n
* Stata sets condzz = condxx for OLS; an unconditional write fails loudly
* if that ever changes.
file write fh "condzz," %21.0g (e(condzz)) _n
file write fh "ll," %21.0g (e(ll)) _n
file write fh "N," %21.0g (e(N)) _n
file write fh "rss," %21.0g (e(rss)) _n
file close fh


display "generate-helpfile-fixtures.do complete"
