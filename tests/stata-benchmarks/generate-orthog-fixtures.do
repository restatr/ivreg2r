/*===========================================================================
  generate-orthog-fixtures.do
  ---------------------------
  Generates CSV benchmark fixtures for orthogonality (instrument-subset)
  C-statistic testing (M-22, planning/22-spec-matrix.md), re-based off the
  Card/sim ad hoc data onto the canonical help-file bases:

    - griliches76 (H03/H25/H26, help.txt:1138/1235/1242): the schooling-IQ
      wage model, orthog(s) [an included regressor] and orthog(age mrt)
      [two excluded IVs], both taken minus gmm2s since the H25/H26 worked
      examples use gmm2s but the orthog C-stat itself is estimator-agnostic
      and 2SLS keeps this family independent of the GMM2S fixture family.
    - abdata (H88, help.txt:1541): the Arellano-Bond employment model with
      first- and second-difference instruments, minus gmm2s, with orthog()
      added as a D5a option-variation to exercise one-way cluster-robust
      orthog.
    - nlswork (H104, help.txt:1628): the two-way-cluster wage model, with
      tenure endogenized against geography instruments as a D5a
      option-variation (WITH DISCLOSURE -- see the nlswork cell below) to
      exercise the two-way-cluster orthog code path, since nlswork has no
      overidentified worked example in the help file.

  `small` variants are deliberately NOT generated for any cell in this file:
  Stata's e(cstat)/e(cstatp)/e(cstatdf) are invariant to the `small` option
  (verified byte-identical in the retired card/sim fixtures this family
  replaces), so a small/non-small pair would just duplicate the same numbers
  under a different suffix. The R tests fit both small = FALSE and
  small = TRUE against each fixture here and assert both reproduce the same
  C-statistic.

  Output directory: tests/stata-benchmarks/fixtures/ (relative to pkg/)

  Usage (CWD must be the package root, i.e. pkg/):
    /Applications/StataNow/StataSE.app/Contents/MacOS/stata-se -b \
      tests/stata-benchmarks/generate-orthog-fixtures.do
===========================================================================*/

clear all
set more off
version 14

// Output directory
local outdir "tests/stata-benchmarks/fixtures"
capture mkdir "`outdir'"


/*===========================================================================
  Load all three datasets FIRST, before the program is defined (bcuse-style
  loaders can issue "clear all" internally; keeping every use/save block
  ahead of "program define" avoids any ordering hazard).
===========================================================================*/

// --- griliches76: local cache, then BC.edu web fallback (per the H03 base
//     model, help.txt:1138) ---
capture use "../validation/data/griliches76.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/hayashi/griliches76.dta, clear
    if _rc {
        display as error "Could not load griliches76 dataset (no local cache, no network)."
        exit 601
    }
}
quietly xi i.year
save "`outdir'/_orthog_gril_temp.dta", replace

// --- abdata: local cache, then BC.edu web fallback (per the H88 base
//     model, help.txt:1541) ---
capture use "../validation/data/abdata.dta", clear
if _rc {
    capture use http://fmwww.bc.edu/ec-p/data/macro/abdata.dta, clear
    if _rc {
        display as error "Could not load abdata dataset (no local cache, no network)."
        exit 601
    }
}
save "`outdir'/_orthog_ab_temp.dta", replace

// --- nlswork: local cache, then webuse fallback (per the H104 base model,
//     help.txt:1628) ---
capture use "../validation/data/nlswork.dta", clear
if _rc {
    capture webuse nlswork, clear
    if _rc {
        display as error "Could not load nlswork dataset (no local cache, no network)."
        exit 601
    }
}
save "`outdir'/_orthog_nls_temp.dta", replace


/*---------------------------------------------------------------------------
  Helper program: save orthog-specific diagnostics to CSV.
---------------------------------------------------------------------------*/
capture program drop save_orthog_results
program define save_orthog_results
    syntax, prefix(string) suffix(string) outdir(string)

    quietly {
        preserve
        clear
        set obs 1

        gen double cstat = .
        gen double cstatp = .
        gen double cstatdf = .
        gen str80 clist = ""
        gen double j = .
        gen double jp = .
        gen double jdf = .

        // Orthogonality C-statistic
        capture replace cstat = e(cstat)
        capture replace cstatp = e(cstatp)
        capture replace cstatdf = e(cstatdf)
        capture replace clist = "`e(clist)'"

        // Overidentification (for reference / sanity checks)
        capture replace j = e(j)
        capture replace jp = e(jp)
        capture replace jdf = e(jdf)
        // IID path uses sargan instead of j
        capture {
            if missing(j) {
                replace j = e(sargan)
                replace jp = e(sarganp)
                replace jdf = e(sargandf)
            }
        }

        export delimited using ///
            "`outdir'/`prefix'_orthog_`suffix'.csv", replace
        restore
    }
end


/*===========================================================================
  BLOCK 1: griliches76
  Base model = help-file H03 (help.txt:1138):
    ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt)
===========================================================================*/
use "`outdir'/_orthog_gril_temp.dta", clear

// --- D5a: H26 (help.txt:1242) minus gmm2s -- 2SLS IID orthog of two
//     excluded IVs; restricted model stays overidentified (overid df
//     3 -> 1), nondegenerate difference-of-Sargan ---
display _newline(2) "=== griliches orthog(age mrt), IID ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), orthog(age mrt)
save_orthog_results, prefix(gril_orthog1) suffix(iid) outdir(`outdir')

// --- D5a: H26 minus gmm2s, robust VCE -- difference-of-J ---
display _newline(2) "=== griliches orthog(age mrt), HC1 ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), orthog(age mrt) robust
save_orthog_results, prefix(gril_orthog1) suffix(hc1) outdir(`outdir')

// --- D5a: H25 (help.txt:1235) minus gmm2s -- orthog of an INCLUDED
//     regressor under 2SLS ---
display _newline(2) "=== griliches orthog(s), IID ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), orthog(s)
save_orthog_results, prefix(gril_orthogs) suffix(iid) outdir(`outdir')

// --- D5a option-variation: dofminus threading into the C-stat ---
display _newline(2) "=== griliches orthog(age mrt), IID, dofminus(1) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), orthog(age mrt) dofminus(1)
save_orthog_results, prefix(gril_orthog1_dof) suffix(iid) outdir(`outdir')

display _newline(2) "=== griliches orthog(age mrt), HC1, dofminus(1) ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), robust orthog(age mrt) dofminus(1)
save_orthog_results, prefix(gril_orthog1_dof) suffix(hc1) outdir(`outdir')

// --- All 4 excluded IVs tested: restricted model underidentified ->
//     Stata reports cstat=0 (guard cell, same role as the retired
//     card_orthog2) ---
display _newline(2) "=== griliches orthog(med kww age mrt), IID ==="
ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age mrt), orthog(med kww age mrt)
save_orthog_results, prefix(gril_orthog2) suffix(iid) outdir(`outdir')


/*===========================================================================
  BLOCK 2: abdata
  Base model = help-file H88 (help.txt:1541), minus gmm2s:
    ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), cluster(id)
===========================================================================*/
use "`outdir'/_orthog_ab_temp.dta", clear
quietly tsset id year

// --- D5a: H88 minus gmm2s plus orthog -- one-way-cluster orthog;
//     restricted overid 3 -> 2 ---
display _newline(2) "=== abdata orthog(d2.w), cluster(id) ==="
ivreg2 n (w k ys = d.w d.k d.ys d2.w d2.k d2.ys), cluster(id) orthog(d2.w)
save_orthog_results, prefix(ab_orthog1) suffix(cl) outdir(`outdir')


/*===========================================================================
  BLOCK 3: nlswork
  Base model = help-file H104 (help.txt:1628):
    ivreg2 ln_w grade age ttl_exp tenure, cluster(idcode year)

  D5a WITH DISCLOSURE: nlswork has no overidentified worked example in the
  help file; this spec endogenizes tenure with geography instruments purely
  to exercise the two-way-cluster orthog code path on the H104 base model.
  It is NOT an economically motivated model.
===========================================================================*/
use "`outdir'/_orthog_nls_temp.dta", clear

display _newline(2) "=== nlswork orthog(south), cluster(idcode year) ==="
ivreg2 ln_wage grade age ttl_exp (tenure = south c_city not_smsa), cluster(idcode year) orthog(south)
save_orthog_results, prefix(nls_orthog1) suffix(cl2) outdir(`outdir')


/*===========================================================================
  Clean up temp files and done
===========================================================================*/
capture erase "`outdir'/_orthog_gril_temp.dta"
capture erase "`outdir'/_orthog_ab_temp.dta"
capture erase "`outdir'/_orthog_nls_temp.dta"
display _newline(2) "=== All orthog fixtures generated ==="
display "Output directory: `outdir'"
