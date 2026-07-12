## Resubmission

This is a resubmission. In this version I have:

* Single-quoted 'Stata' at both of its occurrences in the Description field, as requested.

Regarding the reviewer's note on writing function names as foo(): the Title and Description contain no R function names. The quoted 'ivreg2' in the Description refers to the 'Stata' command of that name, not to a function.

## R CMD check results

0 errors | 0 warnings | 1 note

Checked with `R CMD check --as-cran` on the resubmission tarball (R 4.6.1, macOS arm64, 2026-07-12) and, for the originally submitted tarball (identical apart from the DESCRIPTION quoting fix above), on win-builder under both R-release (4.6.1) and R-devel (2026-07-10 r90234), each returning this single note. The note is the standard "New submission" note (this is the first CRAN submission of this package); the "possibly misspelled" words it lists (Driscoll, Kraay, GMM, HAC, LIML, SLS) are estimator acronyms and author surnames. The package also checks clean on the GitHub Actions matrix: Linux (R release and oldrel-1), macOS (release and oldrel-1), and Windows (release).

The two Stata Journal DOIs cited in DESCRIPTION, CITATION, and the README (10.1177/1536867X0300300101 and 10.1177/1536867X0800700402) return HTTP 403 to automated URL checkers because the publisher (SAGE) blocks non-browser clients; both resolve normally in a browser. Verified 2026-07-11 with `curl -I -L https://doi.org/<doi>` from a residential network: doi.org answers 302 and the SAGE landing page answers 403 to the automated HEAD request.

## Test methodology

This package is inspired by Stata's `ivreg2` command (Baum, Schaffer & Stillman, 2003, 2007), whose authors are credited with ctb and cph roles in Authors@R. Correctness is verified by comparing R output against Stata `ivreg2` output across an extensive benchmark suite covering all estimation methods, VCE types, and diagnostic tests. Every comparison passes within documented numerical tolerances, and the test suite (15,491 tests, including a machine-checked tolerance-policy gate) runs with no failures, warnings, or skips.
