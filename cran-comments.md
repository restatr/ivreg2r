## R CMD check results

0 errors | 0 warnings | 1 note

The one note is the standard "New submission" note (this is the first
CRAN submission of this package).

## Test methodology

This package is inspired by Stata's `ivreg2` command (Baum, Schaffer &
Stillman, 2003, 2007). Correctness is verified by comparing R output
against Stata `ivreg2` output across an extensive benchmark suite
covering all estimation methods, VCE types, and diagnostic tests. Every
comparison passes within documented numerical tolerances, and the test
suite runs with no failures, warnings, or skips.

## First submission

This is the first CRAN submission of this package.
