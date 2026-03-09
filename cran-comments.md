## R CMD check results

0 errors | 0 warnings | 0 notes

## Test methodology

This package is inspired by Stata's `ivreg2` command (Baum, Schaffer &
Stillman, 2003, 2007). Correctness is verified by comparing R output
against Stata `ivreg2` output across 892 fixture files covering all
estimation methods, VCE types, and diagnostic tests. All fixture
comparisons pass within their specified tolerances (8731 tests, 0
failures).

## First submission

This is the first CRAN submission of this package.
