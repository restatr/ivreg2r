## R CMD check results

0 errors | 0 warnings | 1 note

The one note is the standard "New submission" note (this is the first
CRAN submission of this package).

The two Stata Journal DOIs cited in CITATION and the README
(10.1177/1536867X0300300101 and 10.1177/1536867X0800700402) return
HTTP 403 to automated URL checkers because the publisher (SAGE) blocks
non-browser clients; both resolve normally in a browser.

## Test methodology

This package is inspired by Stata's `ivreg2` command (Baum, Schaffer &
Stillman, 2003, 2007). Correctness is verified by comparing R output
against Stata `ivreg2` output across an extensive benchmark suite
covering all estimation methods, VCE types, and diagnostic tests. Every
comparison passes within documented numerical tolerances, and the test
suite runs with no failures, warnings, or skips.

## First submission

This is the first CRAN submission of this package.
