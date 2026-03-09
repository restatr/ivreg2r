# ivreg2r 0.1.0

Initial CRAN release.

## Estimation methods

* Two-stage least squares (2SLS) via three-part formula:
  `y ~ exog | endo | instruments`.
* OLS as a special case via one-part formula: `y ~ x1 + x2`.
* LIML, Fuller, and arbitrary k-class estimators (`method`, `fuller`,
  `kclass` arguments).
* Two-step efficient GMM (`method = "gmm2s"`).
* Continuously-updated estimator (`method = "cue"`).
* User-supplied weighting and covariance matrices (`wmatrix`, `smatrix`).
* Frisch-Waugh-Lovell projection (`partial`).
* Fixed-coefficient evaluation (`b0`).

## Variance-covariance estimators

* Classical (IID), heteroskedasticity-robust, one-way and two-way
  cluster-robust standard errors (`vcov`, `clusters`).
* HAC and AC kernels: Bartlett, Parzen, Quadratic Spectral,
  Tukey-Hanning, Tukey-Hamming, Truncated, Daniell, Tent
  (`kernel`, `bw`).
* Automatic bandwidth selection via Newey-West (1994) (`bw = "auto"`).
* Kiefer (1980) robust VCE (`kiefer`).
* Driscoll-Kraay (1998) standard errors (`dkraay`).
* Thompson (2009) two-way cluster + kernel (`clusters` + `kernel`).
* COVIV option for LIML/k-class models (`coviv`).
* Finite-sample corrections controlled uniformly by `small`.
* Centered moment conditions (`center`).
* PSD corrections for non-positive-semidefinite VCV matrices (`psd`).

## Diagnostics (computed automatically)

* Weak identification: Cragg-Donald Wald F, Kleibergen-Paap rk Wald F,
  Stock-Yogo (2005) critical values for IV/LIML/Fuller.
* Underidentification: Anderson canonical correlations LM,
  Kleibergen-Paap rk LM.
* Overidentification: Sargan (IID), Hansen J (robust/cluster/HAC).
* Anderson-Rubin (1949) test (F and chi-squared forms).
* Stock-Wright (2000) S statistic.
* Endogeneity test (C statistic / difference-in-Sargan).
* First-stage diagnostics: F-test, partial R-squared, Shea partial
  R-squared, Sanderson-Windmeijer conditional F, Angrist-Pischke
  conditional F/chi-squared.
* Anderson-Rubin LIML overidentification (LR and linearized forms).
* Model F-test across all VCE types.
* Orthogonality test for instrument subsets (`orthog`).
* Redundancy test for excluded instruments (`redundant`).
* Reduced-form regression (`reduced_form`).

## Other features

* Analytic, frequency, and probability weights (`weights`, `weight_type`).
* Factor variable support with automatic expansion.
* Degree-of-freedom adjustments (`dofminus`, `sdofminus`).
* Three-part collinearity detection with endogenous reclassification.
* S3 methods: `summary`, `print`, `coef`, `vcov`, `residuals`, `fitted`,
  `predict`, `confint`, `nobs`, `formula`, `model.matrix`, `terms`,
  `update`.
* Broom methods: `tidy()`, `glance()`, `augment()`.
* All outputs verified against Stata's `ivreg2` within tight numerical
  tolerances (892 fixture files, 8700+ test expectations).
