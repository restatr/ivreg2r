# ivreg2r 0.1.0

Initial CRAN release. `ivreg2r` provides extended instrumental variables and GMM estimation with automatic diagnostics, inspired by Stata's `ivreg2` (Baum, Schaffer & Stillman).

* Estimators: 2SLS, LIML, Fuller, k-class, two-step efficient GMM, and the continuously-updated estimator (CUE).
* Variance estimators: classical, heteroskedasticity-robust, one- and two-way cluster-robust, HAC/AC with eight kernels and automatic bandwidth selection, Kiefer, and Driscoll-Kraay.
* Diagnostics reported at estimation time: weak identification (Kleibergen-Paap, Cragg-Donald, Stock-Yogo critical values), underidentification, overidentification (Sargan, Hansen J, Stock-Wright S), first-stage tests (Sanderson-Windmeijer, Angrist-Pischke), endogeneity, and orthogonality.
* Tidyverse integration: a three-part formula interface and `tidy()`, `glance()`, and `augment()` methods, plus `diagnostics()` and `first_stage()` accessors that return the specification tests and first-stage results as tidy objects.
* Bundled datasets used throughout the help-file examples and vignettes: `card`, `mroz`, `wagepan`, `griliches`, `klein`, `grunfeld`, `abdata`, `nlswork`, `phillips`, `stockwatson`, and `cigar`.
* All outputs verified against Stata `ivreg2` within tight numerical tolerances.
