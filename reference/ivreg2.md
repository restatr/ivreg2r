# Extended Instrumental Variables Estimation

Estimate models by OLS, two-stage least squares (2SLS), LIML, Fuller, or
k-class with automatic diagnostic tests. Uses a three-part formula for
IV: `y ~ exog | endo | instruments`.

## Usage

``` r
ivreg2(
  formula,
  data,
  weights,
  subset,
  na.action = stats::na.omit,
  vcov = "iid",
  clusters = NULL,
  endog = NULL,
  orthog = NULL,
  redundant = NULL,
  method = "2sls",
  kclass = NULL,
  fuller = 0,
  coviv = FALSE,
  small = FALSE,
  dofminus = 0L,
  sdofminus = 0L,
  weight_type = "aweight",
  kernel = NULL,
  bw = NULL,
  tvar = NULL,
  ivar = NULL,
  kiefer = FALSE,
  dkraay = NULL,
  wmatrix = NULL,
  smatrix = NULL,
  b0 = NULL,
  noid = FALSE,
  partial = NULL,
  nopartialsmall = FALSE,
  center = FALSE,
  psd = NULL,
  sw = FALSE,
  reduced_form = "none",
  first_stage = FALSE,
  model = TRUE,
  x = FALSE,
  y = TRUE
)
```

## Arguments

- formula:

  A formula: `y ~ exog` (OLS), `y ~ exog | endo | instruments` (IV), or
  `y ~ exog | 0 | instruments` (no endogenous regressors, with the
  excluded instruments contributing surplus moment conditions — Stata's
  `(=z1 z2)` form). The latter is estimated by OLS; the surplus
  conditions drive the Sargan/Hansen J overidentification test and
  `orthog` C-tests, and `method = "gmm2s"` with `vcov = "robust"` gives
  Cragg's (1983) heteroskedastic OLS (HOLS) estimator (with the default
  iid VCE the two-step weighting matrix is proportional to
  \\(Z'Z)^{-1}\\ and the estimates equal OLS exactly). The response must
  be numeric; a factor or character response is rejected with an error.

- data:

  A data frame containing the variables in the formula.

- weights:

  Optional analytic weights expression (evaluated in `data`), equivalent
  to Stata's `[aw=varname]`. Must be strictly positive. For aweights and
  pweights, the weights are normalized internally to sum to N, following
  Stata's convention. This makes sigma (RMSE) scale-invariant for those
  types: multiplying all weights by a constant changes no coefficients,
  standard errors, or test statistics. Frequency and importance weights
  are *not* normalized — they redefine N as `sum(weights)` (see
  `weight_type`), so their scale matters by construction. See
  [ivreg2r-conventions](https://restatr.com/ivreg2r/reference/ivreg2r-conventions.md)
  for how weights interact with sigma and for the comparison to
  [`lm()`](https://rdrr.io/r/stats/lm.html) (coefficients match
  [`lm()`](https://rdrr.io/r/stats/lm.html) for OLS; SEs and the VCV
  match only when `small = TRUE`).

- subset:

  Optional subset expression (evaluated in `data`).

- na.action:

  Function for handling `NA`s (default
  [na.omit](https://rdrr.io/r/stats/na.fail.html)). Use
  [na.exclude](https://rdrr.io/r/stats/na.fail.html) to drop incomplete
  cases during estimation but reinsert `NA`s at the omitted positions in
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html), and
  [`predict()`](https://rdrr.io/r/stats/predict.html) output, so results
  align with the original data frame. Stata drops missing observations
  silently (equivalent to
  [na.omit](https://rdrr.io/r/stats/na.fail.html)) and has no analogue
  of [na.exclude](https://rdrr.io/r/stats/na.fail.html).

- vcov:

  Character: covariance type. One of `"iid"` (classical), `"robust"`
  (White heteroskedasticity-robust), `"HAC"` (heteroskedasticity and
  autocorrelation consistent), or `"AC"` (autocorrelation consistent).
  Matching is case-insensitive (`"hac"`, `"Hac"`, and `"HAC"` are all
  accepted), but the value is always normalized to and printed as the
  canonical spelling shown above. Use `small = TRUE` to apply
  finite-sample corrections. To match Stata's `ivreg2, robust`: use
  `vcov = "robust"`. To match Stata's `ivreg2, robust small`: use
  `vcov = "robust", small = TRUE`.

- clusters:

  One-sided formula specifying one or two cluster variables (e.g.
  `~ firmid` for one-way, `~ firmid + year` for two-way). Two-way
  clustering uses the Cameron-Gelbach-Miller (2006) formula. The
  effective cluster count is `min(M1, M2)` per Stata convention. The
  `small` argument controls whether the finite-sample correction
  `(N-1)/(N-K) * M/(M-1)` is applied (matching Stata's `cluster() small`
  combination). Unlike Stata's `ivreg2`, which silently drops
  observations with missing cluster values from the estimation sample,
  ivreg2r raises an error; the user should drop or handle those rows
  explicitly (e.g., by filtering them out) so the estimation sample is
  never changed silently.

- endog:

  Character vector of endogenous regressor names to test for exogeneity
  (endogeneity test / C-statistic). If `NULL` (default), tests all
  endogenous regressors. Names must match variables in the endogenous
  part of the formula. Ignored for OLS models.

  **Note:** Unlike Stata's `ivreg2`, which only computes the endogeneity
  test when the `endog()` option is explicitly specified, `ivreg2r`
  computes it automatically for all IV models.

- orthog:

  Character vector of instrument names to test for orthogonality
  (instrument-subset C-statistic). Names must be included or excluded
  instruments (not endogenous regressors or the intercept). If `NULL`
  (default), no orthogonality test is computed. Ignored for OLS models.
  Equivalent to Stata's `orthog()` option.

- redundant:

  Character vector of excluded instrument names to test for redundancy
  (zero first-stage explanatory power). The test is a KP rk LM test of
  H0: rank=0 on the first-stage coefficient matrix for the tested
  instruments, conditional on maintained instruments. If `NULL`
  (default), no redundancy test is computed. On a model with no
  endogenous regressors (a one-part OLS formula, or the IV form
  `y ~ exog | 0 | instruments`), a warning reports that the test is
  skipped; Stata's `redundant()` is silently ignored in those cases.
  Equivalent to Stata's `redundant()` option.

- method:

  Character: estimation method. One of `"2sls"` (default), `"liml"`,
  `"kclass"`, `"gmm2s"` (two-step efficient GMM), or `"cue"`
  (continuously updated GMM estimator). For OLS models (1-part formula),
  this is ignored. When `fuller > 0` is specified, method is
  automatically promoted to `"liml"`. `"gmm2s"` uses the inverse of the
  moment covariance matrix as the optimal weighting matrix, yielding
  efficient estimates under the specified error structure. Incompatible
  with `fuller` and `kclass`. `"cue"` continuously updates both moment
  conditions and moment covariance at each candidate beta during
  optimization. Under iid errors (no robust/cluster/kernel VCE), CUE
  with no other options is equivalent to LIML with `coviv = TRUE`. Under
  non-iid errors, CUE generalizes LIML to produce efficient estimates.
  Incompatible with `fuller`, `kclass`, `wmatrix`, and `smatrix`.

  **CUE optimizer note:** This package uses R's
  [`optim()`](https://rdrr.io/r/stats/optim.html) (BFGS + Nelder-Mead,
  scaled by the starting values) while Stata uses Mata's
  [`optimize()`](https://rdrr.io/r/stats/optimize.html) (modified
  Newton-Raphson). The CUE objective is non-convex and can have multiple
  local minima, and it can possess economically degenerate minima (e.g.,
  negative R-squared) – a documented property of the CUE ratio objective
  (Hausman, Lewis, Menzel & Newey, 2011), not an optimizer failure. On
  every benchmarked model the two implementations agree on the minimum,
  including one such pathological case, where both converge to the same
  degenerate basin; because those basins are extremely flat, coefficient
  agreement across implementations there is limited to a few digits.
  Inspect the J-statistic and R-squared to judge whether a CUE solution
  is economically sensible.

- kclass:

  Numeric scalar: user-supplied k value for k-class estimation. When
  supplied, `method` is automatically set to `"kclass"`. Must be
  non-negative. Cannot be combined with `method = "liml"` or `fuller`.

- fuller:

  Numeric scalar: Fuller (1977) modification parameter. Non-negative.
  `0` (the default) applies no Fuller adjustment; a positive value
  activates the Fuller-modified LIML estimator, which automatically sets
  `method` to `"liml"` and `k = lambda - fuller / (N - L)`. `fuller = 1`
  gives the bias-corrected LIML estimator; `fuller = 4` targets MSE.
  Cannot be combined with `kclass`.

- coviv:

  Logical: if `TRUE`, use the 2SLS bread `(X_hat'X_hat)^{-1}` instead of
  the k-class bread for VCV computation in LIML/k-class estimation. This
  gives the "COVIV" (covariance at the IV estimates) VCV that is robust
  to misspecification of the LIML model. Silently ignored for OLS and
  2SLS. Default `FALSE`.

- small:

  Logical: if `TRUE`, use small-sample corrections (t/F instead of
  z/chi-squared, `N-K` denominator for sigma).

- dofminus:

  Non-negative integer: large-sample degrees-of-freedom adjustment.
  Subtracted from N in large-sample variance formulas (e.g., sigma =
  rss/(N-dofminus)). Useful when fixed effects have been partialled out.
  Equivalent to Stata's `dofminus()` option.

- sdofminus:

  Non-negative integer: small-sample degrees-of-freedom adjustment.
  Subtracted from the residual degrees of freedom alongside K (e.g.,
  df.residual = N - K - dofminus - sdofminus). Useful when partialling
  out regressors. Equivalent to Stata's `sdofminus()` option.

- weight_type:

  Character: type of weights. One of `"aweight"` (analytic weights,
  default), `"fweight"` (frequency weights), `"pweight"`
  (probability/sampling weights), or `"iweight"` (importance weights).

  **aweight**: Normalized to sum to N. Standard for WLS.

  **fweight**: Integer-valued; N is redefined as `sum(weights)`. HC meat
  uses `w * e^2` (linear, not quadratic).

  **pweight**: Normalized to sum to N. Forces robust VCE (overrides
  `vcov = "iid"` to `"robust"` and `vcov = "AC"` — explicit or
  kiefer-implied — to `"HAC"`).

  **iweight**: Not normalized; N is redefined as `sum(weights)` (float).
  All intermediate calculations use float N; posted `nobs` and
  `df.residual` use `floor(N)`. Restricted to IID VCE only —
  incompatible with robust, cluster, HAC, and GMM estimation.

- kernel:

  Character: kernel function for HAC/AC standard errors. One of
  `"bartlett"`, `"parzen"`, `"truncated"`, `"tukey-hanning"`,
  `"tukey-hamming"`, `"qs"` (quadratic spectral), `"daniell"`, or
  `"tent"`. When specified without an explicit `vcov` change, `vcov` is
  automatically set: `"iid"` becomes `"AC"`, `"robust"` becomes `"HAC"`.
  Requires `bw` and `tvar`.

  **Note:** Kleibergen-Paap identification tests (underidentification
  and weak identification) always use the Bartlett kernel internally,
  regardless of the user-specified kernel. This matches a known behavior
  in Stata's `ivreg2`, where `ranktest` hard-codes the Bartlett kernel.

- bw:

  Numeric or `"auto"`: bandwidth for kernel estimation. Must be positive
  numeric, or `"auto"` for automatic selection via Newey-West (1994).
  Auto selection is available for Bartlett, Parzen, and Quadratic
  Spectral kernels only, and is not supported for panel data (`ivar`).
  Required when `kernel` is specified.

- tvar:

  Character: name of the time variable in `data`. Required for HAC/AC
  estimation and for time-series operators
  [`l()`](https://restatr.com/ivreg2r/reference/ts-operators.md)/[`d()`](https://restatr.com/ivreg2r/reference/ts-operators.md)
  in the formula (see
  [`ts-operators`](https://restatr.com/ivreg2r/reference/ts-operators.md)).
  Lags are resolved by time *value* (gap-aware), mirroring Stata's
  `tsset`.

- ivar:

  Character: name of the panel identifier variable in `data`. Optional;
  only needed for panel data (HAC/AC estimation, or so that time-series
  operators lag within panel units).

- kiefer:

  Logical: if `TRUE`, use the Kiefer (1980) VCE —
  autocorrelation-consistent with kernel = Truncated and bandwidth = T
  (the full time span). Requires panel data (`tvar` + `ivar`).
  Incompatible with robust VCE, clustering, explicit kernel or
  bandwidth. Equivalent to Stata's `ivreg2 ..., kiefer`.

- dkraay:

  Positive numeric scalar: bandwidth for Driscoll-Kraay (1998) VCE. When
  specified, clusters on the time variable and applies kernel smoothing
  across time lags, producing standard errors robust to cross-sectional
  dependence. Requires panel data (`tvar` + `ivar`). Incompatible with
  explicit `bw`. If `kernel` is not specified, defaults to Bartlett.
  Equivalent to Stata's `ivreg2 ..., dkraay(3)`.

- wmatrix:

  Numeric matrix: user-supplied L x L weighting matrix for GMM
  estimation, where L is the number of instruments (including exogenous
  regressors). When supplied without `method = "gmm2s"`, produces
  inefficient GMM with full sandwich VCV (`method` becomes `"gmmw"`).
  When combined with `method = "gmm2s"`, used as the first-step
  weighting matrix (second step uses the optimal weighting matrix). Must
  be symmetric. Incompatible with `method = "liml"`, `kclass`, and
  `fuller`. When `method` is not `"gmm2s"`, ignored without robust VCE,
  clustering, or HAC (with a warning). Equivalent to Stata's `wmatrix()`
  option. When used (not ignored), it is echoed back as `fit$W`,
  matching Stata's `e(W)`. A matrix with dimnames is matched to the
  instrument columns by name (reordered as needed, mirroring Stata's
  matsort; names that do not cover the instrument list are an error); an
  unnamed matrix is used positionally, in instrument column order.

- smatrix:

  Numeric matrix: user-supplied L x L moment covariance matrix for GMM
  estimation. When supplied, the GMM estimation uses this matrix instead
  of computing Omega from residuals. Implies efficient GMM (`method` is
  promoted to `"gmm2s"` if `"2sls"`). Must be symmetric. Incompatible
  with `method = "liml"`, `kclass`, and `fuller`. Equivalent to Stata's
  `smatrix()` option. Echoed back as `fit$S` (never recomputed),
  matching Stata's `e(S)`. A matrix with dimnames is matched to the
  instrument columns by name (reordered as needed, mirroring Stata's
  matsort); an unnamed matrix is used positionally, in instrument column
  order.

- b0:

  Numeric vector: evaluate the CUE objective at this fixed parameter
  vector without optimization. When supplied, `method` is promoted to
  `"cue"` and the J(b0) statistic is computed and stored. Identification
  diagnostics (underid, weak-id, first-stage, AR, SW, endogeneity,
  orthog) are suppressed (matching Stata's `b0()` option). Length must
  equal the number of regressors K. Named vectors are reordered to match
  model matrix columns; unnamed vectors are used as-is in model matrix
  column order. Incompatible with `method = "gmm2s"`, `"liml"`,
  `kclass`, `fuller`, and `wmatrix`.

- noid:

  Logical: if `TRUE`, suppress computation of underidentification, weak
  identification, and redundancy test statistics. Overidentification,
  Anderson-Rubin, Stock-Wright, endogeneity, and first-stage diagnostics
  are still computed. Default `FALSE`. Useful for speeding up simulation
  loops where rank-based identification tests are expensive. Equivalent
  to Stata's `noid` option.

- partial:

  Character vector: exogenous regressors to partial out via
  Frisch-Waugh-Lovell projection before estimation. Coefficients on
  partialled variables are not recoverable and are not reported. Special
  values:

  - `"_cons"` or `"(Intercept)"`: partial only the constant (demean).

  - `"_all"`: partial all included exogenous regressors.

  - `NULL` (default): no partialling.

  By default, `sdofminus` is automatically incremented by the number of
  partialled variables (including the constant). Use
  `nopartialsmall = TRUE` to suppress this adjustment.

  **Note:** FWL invariance holds for OLS, 2SLS, LIML, and two-step GMM,
  but **not** for CUE. The CUE objective recomputes the moment
  covariance at each iteration, so partialled residuals produce a
  different optimization surface. A warning is issued when `partial` is
  used with `method = "cue"`.

  After partialling, [`predict()`](https://rdrr.io/r/stats/predict.html)
  is restricted to residuals only (no `newdata`). Summary output notes
  that total SS, model F, and R-squared are partial-model values.

  Equivalent to Stata's `partial()` option.

- nopartialsmall:

  Logical: if `TRUE`, suppress the automatic `sdofminus` adjustment from
  partialling. Default `FALSE`. Equivalent to Stata's `nopartialsmall`
  option.

- center:

  Logical: if `TRUE`, subtract the mean of moment conditions (scores)
  before computing the S matrix (meat of the sandwich VCE). Default
  `FALSE`. Centering only affects non-homoskedastic VCE types (HC,
  cluster, HAC); a warning is issued if used with IID or AC VCE.

  **Note:** Centering is applied to the main model's VCE and to
  diagnostic tests that use the main model's S matrix
  (overidentification, orthogonality). However, the endogeneity test
  (`endog`) computes its own S matrix from the restricted model
  *without* centering, even when `center = TRUE`. This matches Stata's
  `ivreg2`, where `center` is not forwarded to the recursive call for
  the endogeneity test.

- psd:

  Character or NULL: PSD correction for the moment covariance matrix S.
  `NULL` (default) applies no correction. `"psd0"` zeroes negative
  eigenvalues (Politis 2007). `"psda"` replaces negative eigenvalues
  with their absolute values (Stock & Watson 2008). The correction is
  applied to S *before* the variance-covariance matrix is assembled from
  it, matching Stata's `m_omega`; the corrected S is stored as `$S`. The
  conventional (`vcov = "iid"`) VCV and the Kleibergen-Paap
  identification statistics are never corrected, also matching Stata
  (its `ranktest` does not receive the psd option). A warning is emitted
  whenever negative eigenvalues are detected and corrected (Stata
  corrects silently). Ignored when a user-supplied `smatrix` is given:
  the user matrix is used as-is for estimation, the VCV, and `$S`,
  matching Stata (whose `m_omega` is never invoked for a supplied S).
  Equivalent to Stata's `psd0` and `psda` options.

- sw:

  Logical: if `TRUE`, compute the Stock-Watson (2008, Econometrica)
  panel-robust VCE. Requires panel data (`ivar`). Incompatible with
  clustering, HAC kernels, kiefer, dkraay, fweights, and iweights.
  Forces `vcov = "robust"` automatically. Equivalent to Stata's `sw`
  option (labeled "BETA VERSION" in Stata).

  **Precondition:** the Stock-Watson (2008) correction is derived for a
  within-transformed (fixed-effects) panel regression with fixed `T` and
  large `N`. `ivreg2()` does not perform the within transformation for
  you: within-transform the data yourself (demean every variable by
  panel unit) before fitting, and pass `dofminus` equal to the number of
  panel units, which charges the absorbed unit means to the degrees of
  freedom. See the worked example in the "Fixed-effects panels: the
  Stock-Watson correction" section of
  [`vignette("time-series-gmm")`](https://restatr.com/ivreg2r/articles/time-series-gmm.md).

- reduced_form:

  Character: what reduced-form output to store. `"none"` (default)
  stores nothing. `"rf"` stores the y ~ Z regression (equivalent to
  Stata's `saverf`; Stata's `rf` option *displays* the reduced form
  without storing it). `"system"` stores the full system of y + all
  endogenous variables regressed on Z, with cross-equation VCV
  (equivalent to Stata's `savesfirst`, an option present in `ivreg2.ado`
  but absent from its help file; Stata's `sfirst` displays the system
  without storing it). On a model with no endogenous regressors (a
  one-part OLS formula, or the IV form `y ~ exog | 0 | instruments`), no
  reduced form is computed — matching Stata, which skips reduced-form
  estimation whenever the endogenous list is empty — and a warning
  reports that the request was ignored.

- first_stage:

  Logical: if `TRUE`, store extractable first-stage regression objects
  on the fitted model. Access them via
  [`first_stage()`](https://restatr.com/ivreg2r/reference/first_stage.md);
  see
  [`first_stage()`](https://restatr.com/ivreg2r/reference/first_stage.md)
  for the full list of supported S3 methods
  ([`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html), and
  others). Equivalent to Stata's `savefirst` option. Default `FALSE`. On
  a model with no endogenous regressors (a one-part OLS formula, or the
  IV form `y ~ exog | 0 | instruments`), nothing is stored and a warning
  reports that the request was ignored.

- model:

  Logical: if `TRUE` (default), store the model frame in the return
  object.

- x:

  Logical: if `TRUE`, store model matrices (`X`, `Z`) in the return
  object.

- y:

  Logical: if `TRUE` (default), store the response vector in the return
  object.

## Value

An object of class `"ivreg2"`. Beyond the usual components
(`coefficients`, `residuals`, `vcov`, `diagnostics`, ...), the object
stores the estimated moment-condition covariance matrix as `$S` and the
GMM weighting matrix as `$W`, corresponding to Stata's `e(S)` and
`e(W)`:

- `$S` is present for every fit. Normalization matches Stata's
  `m_omega`: iid `sigma^2 (Z'Z)/N` with `sigma^2 = RSS/(N - dofminus)`;
  robust and HAC, meat divided by `N - dofminus`; cluster, meat divided
  by `N`; psd-corrected when `psd` is set; built from centered moments
  when `center = TRUE`. A user-supplied `smatrix` is echoed, never
  recomputed and never psd-corrected.

- `$W` is `NULL` for `method = "liml"` and `"kclass"` (as in Stata, no
  weighting matrix is defined for k-class estimators). For two-step GMM,
  CUE, and `b0` fits it is the inverse of the defining `$S`; for
  one-step OLS/2SLS it is `(1/sigma^2) (Z'Z/N)^{-1}`; a user-supplied
  `wmatrix` is echoed (it is the first-step weighting matrix under
  `method = "gmm2s"`).

Rows and columns are named by the R instrument columns — `(Intercept)`
rather than Stata's `_cons`, and R column order (exogenous regressors
before excluded instruments) rather than Stata's. Passing `fit$S` back
via `smatrix` reproduces the corresponding efficient-GMM fit, mirroring
the `e(S)` reuse workflow in Stata's `ivreg2` help file.

## Overidentification and endogeneity tests

When the model carries more instruments than regressors, `ivreg2()`
automatically reports an overidentification test: the Sargan (1958)
statistic under `vcov = "iid"`, and the Hansen J statistic under robust,
cluster, or HAC errors. The test evaluates the `L - K` overidentifying
restrictions jointly, where `L` is the number of instruments and `K` the
number of regressors. A rejection indicates that at least one instrument
violates its orthogonality condition, without identifying which one. A
non-rejection is necessary but not sufficient for instrument validity.
The test has no power against a bias shared by every instrument, so
orthogonality conditions that fail in the same direction can still
appear jointly satisfied. The test is also uninformative about the `K`
exactly-identifying moment conditions, whose validity must be defended
on substantive grounds. An exactly identified model has no
overidentifying restrictions, and therefore no overidentification test.

The endogeneity test, also called the C-statistic or difference test,
asks whether regressors treated as endogenous could instead be treated
as exogenous. It is the overidentification statistic of the model that
adds the suspect regressors to the instrument set, treating them as
exogenous, minus that of the model that treats them as endogenous. Under
the null of exogeneity, this difference follows a chi-square
distribution with degrees of freedom equal to the number of regressors
tested. Both overidentification statistics are formed from the same
estimated moment covariance, which guarantees that the difference is
non-negative (Hayashi 2000). `ivreg2()` computes this test automatically
for every IV model, whereas Stata computes it only when the `endog()`
option is given explicitly. The `endog` argument selects which
endogenous regressors are tested. The statistic takes the
difference-in-Sargan form under iid errors and the
difference-in-Hansen-J form under robust, cluster, or HAC errors.

## LIML, Fuller, and k-class estimators

LIML, Fuller's modification, and the general k-class estimators are
alternatives to 2SLS motivated by weak instruments. LIML is
approximately median-unbiased, and its Stock-Yogo weak-identification
thresholds are more forgiving than those for 2SLS. LIML has no
finite-sample moments of any order, however. Its sampling distribution
can therefore be dispersed, and individual estimates can lie far from
the truth. Fuller's (1977) modification sets
`k = lambda - fuller / (N - L)`, which yields an estimator with finite
moments. The choice `fuller = 1` is approximately best-unbiased, and
`fuller = 4` minimizes mean squared error, both within the standard
weak-instrument framework. Weak-identification pretesting for LIML and
Fuller uses estimator-specific critical values rather than the 2SLS
values. LIML coincides with 2SLS when the model is exactly identified,
because the LIML eigenvalue equals one in that case. Fuller does not,
since it subtracts `fuller / (N - L)` from that eigenvalue.

## See also

[ivreg2r-conventions](https://restatr.com/ivreg2r/reference/ivreg2r-conventions.md)
for the statistical conventions (sigma normalization, `small`, weights,
degrees of freedom) and how they differ from
[`lm()`](https://rdrr.io/r/stats/lm.html);
[ivreg2r-glossary](https://restatr.com/ivreg2r/reference/ivreg2r-glossary.md)
for the diagnostic-output acronyms;
[`summary.ivreg2()`](https://restatr.com/ivreg2r/reference/summary.ivreg2.md)
and
[`print.ivreg2()`](https://restatr.com/ivreg2r/reference/print.ivreg2.md)
for console output;
[`tidy.ivreg2()`](https://restatr.com/ivreg2r/reference/tidy.ivreg2.md),
[`glance.ivreg2()`](https://restatr.com/ivreg2r/reference/glance.ivreg2.md),
and
[`augment.ivreg2()`](https://restatr.com/ivreg2r/reference/augment.ivreg2.md)
for broom-style output;
[`first_stage()`](https://restatr.com/ivreg2r/reference/first_stage.md)
for first-stage regressions. The vignettes
[`vignette("introduction", "ivreg2r")`](https://restatr.com/ivreg2r/articles/introduction.md),
[`vignette("advanced-iv", "ivreg2r")`](https://restatr.com/ivreg2r/articles/advanced-iv.md),
and
[`vignette("time-series-gmm", "ivreg2r")`](https://restatr.com/ivreg2r/articles/time-series-gmm.md)
work through applied examples.

Other ivreg2 methods:
[`coef.ivreg2()`](https://restatr.com/ivreg2r/reference/coef.ivreg2.md),
[`confint.ivreg2()`](https://restatr.com/ivreg2r/reference/confint.ivreg2.md),
[`diagnostics()`](https://restatr.com/ivreg2r/reference/diagnostics.md),
[`first_stage()`](https://restatr.com/ivreg2r/reference/first_stage.md),
[`fitted.ivreg2()`](https://restatr.com/ivreg2r/reference/fitted.ivreg2.md),
[`formula.ivreg2()`](https://restatr.com/ivreg2r/reference/formula.ivreg2.md),
[`model.matrix.ivreg2()`](https://restatr.com/ivreg2r/reference/model.matrix.ivreg2.md),
[`nobs.ivreg2()`](https://restatr.com/ivreg2r/reference/nobs.ivreg2.md),
[`predict.ivreg2()`](https://restatr.com/ivreg2r/reference/predict.ivreg2.md),
[`print.ivreg2()`](https://restatr.com/ivreg2r/reference/print.ivreg2.md),
[`print.summary.ivreg2()`](https://restatr.com/ivreg2r/reference/print.summary.ivreg2.md),
[`residuals.ivreg2()`](https://restatr.com/ivreg2r/reference/residuals.ivreg2.md),
[`summary.ivreg2()`](https://restatr.com/ivreg2r/reference/summary.ivreg2.md),
[`terms.ivreg2()`](https://restatr.com/ivreg2r/reference/terms.ivreg2.md),
[`update.ivreg2()`](https://restatr.com/ivreg2r/reference/update.ivreg2.md),
[`vcov.ivreg2()`](https://restatr.com/ivreg2r/reference/vcov.ivreg2.md)

## Examples

``` r
data(mroz)
mroz_work <- subset(mroz, inlf == 1)

# --- Just-identified IV: return to schooling ---
# Wooldridge (2020), Example 15.1: father's education instruments
# education in a simple wage equation (replicates the published
# estimates: educ = 0.059, SE = 0.035). With one instrument for one
# endogenous variable, no overid test is possible; instrument validity
# must be defended on substantive grounds.
fit <- ivreg2(lwage ~ 1 | educ | fatheduc,
              data = mroz_work)
summary(fit)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ 1 | educ | fatheduc, data = mroz_work)
#> 
#> Observations: 428 
#> VCV type:     Classical (iid) 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  0.44110    0.44506   0.991   0.3216  
#> educ         0.05917    0.03506   1.688   0.0915 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.0934 
#> Adj. R-squared: 0.0913 
#> Wald chi2(1):  2.8 (p = 0.0929)
#> Root MSE:       0.6878 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(1) = 73.86 (p < 2.2e-16)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           88.84 
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       16.38 
#>      15%  maximal IV size       8.96 
#>      20%  maximal IV size       6.66 
#>      25%  maximal IV size       5.53 
#> 
#> Overidentification test (Sargan):  (equation exactly identified)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(1,426) = 2.59 (p = 0.1086)
#>   Anderson-Rubin Wald Chi-sq(1) = 2.60 (p = 0.1070)
#>   Stock-Wright LM S Chi-sq(1) = 2.58 (p = 0.1081)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 2.47 (p = 0.1158)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ               88.84    0.0000      0.1726    0.1726     88.84     88.84
#> 
#> Instrumented:          educ 
#> Excluded instruments:  fatheduc 
#> 

# --- Overidentified IV: testing instrument validity ---
# The Stata ivreg2 help file's illustrative baseline (line 1274):
# age, kidslt6, and kidsge6 give two overidentifying restrictions,
# so the Sargan test can detect misspecification. Note its weak
# first stage (Cragg-Donald F ~ 4.3) -- see the LIML example below.
fit_overid <- ivreg2(lwage ~ exper + expersq | educ |
                       age + kidslt6 + kidsge6, data = mroz_work)
summary(fit_overid)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq | educ | age + kidslt6 + 
#>     kidsge6, data = mroz_work)
#> 
#> Observations: 428 
#> VCV type:     Classical (iid) 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)   
#> (Intercept) -0.3848718  1.0115512  -0.380  0.70359   
#> educ         0.0964002  0.0814278   1.184  0.23646   
#> exper        0.0421930  0.0138831   3.039  0.00237 **
#> expersq     -0.0008323  0.0004204  -1.980  0.04773 * 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.1556 
#> Adj. R-squared: 0.1496 
#> Wald chi2(3):  7.5 (p = 0.0001)
#> Root MSE:       0.6638 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(3) = 12.82 (p = 0.0051)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           4.34 
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       22.30 
#>      15%  maximal IV size       12.83 
#>      20%  maximal IV size       9.54 
#>      25%  maximal IV size       7.80 
#>   Stock-Yogo critical values (IV relative bias):
#>      5%  maximal IV relative bias   13.91 
#>      10%  maximal IV relative bias   9.08 
#>      20%  maximal IV relative bias   6.46 
#>      30%  maximal IV relative bias   5.39 
#> 
#> Overidentification test (Sargan):
#>   Chi-sq(2) = 0.70 (p = 0.7042)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(3,422) = 0.61 (p = 0.6076)
#>   Anderson-Rubin Wald Chi-sq(3) = 1.86 (p = 0.6016)
#>   Stock-Wright LM S Chi-sq(3) = 1.85 (p = 0.6033)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 0.02 (p = 0.8899)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ                4.34    0.0050      0.0299    0.0299      4.34      4.34
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq 
#> Excluded instruments:  age, kidslt6, kidsge6 
#> 

# --- Robust standard errors ---
# With heteroskedasticity, Kleibergen-Paap diagnostics replace
# Anderson/Cragg-Donald, and Hansen J replaces Sargan.
fit_robust <- ivreg2(lwage ~ exper + expersq | educ |
                       age + kidslt6 + kidsge6, data = mroz_work,
                       vcov = "robust")
summary(fit_robust)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq | educ | age + kidslt6 + 
#>     kidsge6, data = mroz_work, vcov = "robust")
#> 
#> Observations: 428 
#> VCV type:     Robust 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)  
#> (Intercept) -0.3848718  1.0599330  -0.363   0.7165  
#> educ         0.0964002  0.0864626   1.115   0.2649  
#> exper        0.0421930  0.0166585   2.533   0.0113 *
#> expersq     -0.0008323  0.0004707  -1.768   0.0770 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.1556 
#> Adj. R-squared: 0.1496 
#> Wald chi2(3):  6.0 (p = 0.0005)
#> Root MSE:       0.6638 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(3) = 11.23 (p = 0.0105)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           4.34 
#>   Kleibergen-Paap rk Wald F:     5.02 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       22.30 
#>      15%  maximal IV size       12.83 
#>      20%  maximal IV size       9.54 
#>      25%  maximal IV size       7.80 
#>   Stock-Yogo critical values (IV relative bias):
#>      5%  maximal IV relative bias   13.91 
#>      10%  maximal IV relative bias   9.08 
#>      20%  maximal IV relative bias   6.46 
#>      30%  maximal IV relative bias   5.39 
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(2) = 0.51 (p = 0.7734)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(3,422) = 0.58 (p = 0.6287)
#>   Anderson-Rubin Wald Chi-sq(3) = 1.76 (p = 0.6229)
#>   Stock-Wright LM S Chi-sq(3) = 1.66 (p = 0.6454)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 0.00 (p = 0.9712)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ                5.02    0.0020      0.0299    0.0299      5.02      5.02
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq 
#> Excluded instruments:  age, kidslt6, kidsge6 
#> 

# \donttest{
# --- LIML ---
# Same baseline as above (a documented variation on the help-file
# spec). Its weak first stage (Cragg-Donald F ~ 4.3) is the setting
# where 2SLS bias toward OLS grows with the overidentification degree
# relative to instrument strength, and LIML's approximate
# median-unbiasedness -- established within the iid Stock-Yogo
# framework -- is attractive. LIML equals 2SLS when just-identified.
fit_liml <- ivreg2(lwage ~ exper + expersq | educ |
                     age + kidslt6 + kidsge6, data = mroz_work,
                     method = "liml")
summary(fit_liml)
#> 
#> LIML Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq | educ | age + kidslt6 + 
#>     kidsge6, data = mroz_work, method = "liml")
#> 
#> Observations: 428 
#> VCV type:     Classical (iid) 
#> lambda:       1.001642 
#> kclass:       1.001642 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)   
#> (Intercept) -0.3769294  1.0394246  -0.363  0.71688   
#> educ         0.0957581  0.0836906   1.144  0.25254   
#> exper        0.0422292  0.0139270   3.032  0.00243 **
#> expersq     -0.0008335  0.0004220  -1.975  0.04827 * 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.1555 
#> Adj. R-squared: 0.1495 
#> Wald chi2(3):  7.5 (p = 0.0001)
#> Root MSE:       0.6638 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(3) = 12.82 (p = 0.0051)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           4.34 
#>   Stock-Yogo critical values (LIML size):
#>      10%  maximal LIML size      6.46 
#>      15%  maximal LIML size      4.36 
#>      20%  maximal LIML size      3.69 
#>      25%  maximal LIML size      3.32 
#> 
#> Overidentification test (Sargan):
#>   Chi-sq(2) = 0.70 (p = 0.7042)
#> 
#> Anderson-Rubin overidentification:
#>   LR Chi-sq(2) = 0.702 (p = 0.7040)
#>   Linearized Chi-sq(2) = 0.703 (p = 0.7038)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(3,422) = 0.61 (p = 0.6076)
#>   Anderson-Rubin Wald Chi-sq(3) = 1.86 (p = 0.6016)
#>   Stock-Wright LM S Chi-sq(3) = 1.85 (p = 0.6033)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 0.02 (p = 0.8899)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ                4.34    0.0050      0.0299    0.0299      4.34      4.34
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq 
#> Excluded instruments:  age, kidslt6, kidsge6 
#> 

# --- Endogeneity test ---
# Is education actually endogenous? The C-statistic tests
# H0: OLS is consistent (educ is exogenous).
# ivreg2 help file line 1283: the endog(educ) example on this equation.
fit_endog <- ivreg2(lwage ~ exper + expersq | educ |
                      age + kidslt6 + kidsge6, data = mroz_work,
                      endog = "educ")
if (requireNamespace("dplyr", quietly = TRUE)) {
  diagnostics(fit_endog) |> dplyr::filter(test == "endogeneity")
}
#> # A tibble: 1 × 7
#>   test        test_name   statistic    df   df2 p_value tested_vars
#>   <chr>       <chr>           <dbl> <int> <int>   <dbl> <chr>      
#> 1 endogeneity Endogeneity    0.0191     1    NA   0.890 educ       

# --- Clustering ---
# Griliches (1976) wage equation, cluster on year.
# Adapted from Stata `ivreg2` help-file example (line 1250): the help
# command also includes year dummies (xi i.year), omitted here, and no
# small option. Either way the fit warns -- 7 year clusters < the
# instrument count makes the moment covariance singular, which is the
# help file's own lead-in to partialling; see
# vignette("time-series-gmm").
data(griliches)
fit_cl <- ivreg2(lw ~ s + expr + tenure + rns + smsa |
                   iq | med + kww + age,
                 data = griliches, clusters = ~year, small = TRUE)
#> Warning: Hansen J statistic not computed; the moment covariance or Hessian matrix is singular.
#> Warning: Stock-Wright: Omega is rank-deficient; S statistic not computed.
summary(fit_cl)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lw ~ s + expr + tenure + rns + smsa | iq | med + 
#>     kww + age, data = griliches, clusters = ~year, small = TRUE)
#> 
#> Observations: 758 
#> VCV type:     Cluster-robust, small-sample corrected 
#> Clusters:    7 (year)
#> 
#> Coefficients:
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  4.707083   0.461224  10.206 5.16e-05 ***
#> iq          -0.009490   0.005897  -1.609 0.158704    
#> s            0.131165   0.013092  10.019 5.73e-05 ***
#> expr         0.034573   0.016874   2.049 0.086375 .  
#> tenure       0.039682   0.011483   3.456 0.013539 *  
#> rns         -0.111709   0.049898  -2.239 0.066460 .  
#> smsa         0.148145   0.018974   7.808 0.000233 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.2415 
#> Adj. R-squared: 0.2354 
#> F(6, 6):     2398.4 (p = 0.0000)
#> Root MSE:       0.3751 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(3) = 5.86 (p = 0.1188)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           13.75 
#>   Kleibergen-Paap rk Wald F:     38.47 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       22.30 
#>      15%  maximal IV size       12.83 
#>      20%  maximal IV size       9.54 
#>      25%  maximal IV size       7.80 
#>   Stock-Yogo critical values (IV relative bias):
#>      5%  maximal IV relative bias   13.91 
#>      10%  maximal IV relative bias   9.08 
#>      20%  maximal IV relative bias   6.46 
#>      30%  maximal IV relative bias   5.39 
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(2) =  NA (p = NA)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(3,6) = 87.30 (p = 0.0000)
#>   Anderson-Rubin Wald Chi-sq(3) = 308.81 (p < 2.2e-16)
#> 
#> Endogeneity test:  not computed (rank-deficient S)
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   iq                 38.47    0.0003      0.0522    0.0522     38.47     38.47
#> 
#> Instrumented:          iq 
#> Included instruments:  s, expr, tenure, rns, smsa 
#> Excluded instruments:  med, kww, age 
#> 

# --- Two-step efficient GMM ---
# ivreg2 help file line 1154: efficient GMM with robust weighting.
# The optimal weighting matrix uses the first-step robust moment
# covariance; Hansen J replaces Sargan.
fit_gmm <- ivreg2(lw ~ s + expr + tenure + rns + smsa + factor(year) |
                    iq | med + kww + age + mrt,
                  data = griliches, method = "gmm2s", vcov = "robust")
summary(fit_gmm)
#> 
#> 2-Step GMM Estimation
#> 
#> Call:
#> ivreg2(formula = lw ~ s + expr + tenure + rns + smsa + factor(year) | 
#>     iq | med + kww + age + mrt, data = griliches, vcov = "robust", 
#>     method = "gmm2s")
#> 
#> Observations: 758 
#> VCV type:     Robust 
#>               Estimates efficient for arbitrary heteroskedasticity
#> 
#> Coefficients:
#>                 Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)     4.436784   0.289950  15.302  < 2e-16 ***
#> iq             -0.001401   0.004113  -0.341 0.733314    
#> s               0.076835   0.013186   5.827 5.64e-09 ***
#> expr            0.031234   0.006693   4.667 3.06e-06 ***
#> tenure          0.049000   0.007344   6.672 2.52e-11 ***
#> rns            -0.100681   0.029589  -3.403 0.000667 ***
#> smsa            0.133597   0.026325   5.075 3.87e-07 ***
#> factor(year)67 -0.021013   0.045543  -0.461 0.644515    
#> factor(year)68  0.089099   0.042702   2.087 0.036930 *  
#> factor(year)69  0.207248   0.040800   5.080 3.78e-07 ***
#> factor(year)70  0.233831   0.052851   4.424 9.67e-06 ***
#> factor(year)71  0.234552   0.042566   5.510 3.58e-08 ***
#> factor(year)73  0.336027   0.040410   8.315  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.4166 
#> Adj. R-squared: 0.4072 
#> Wald chi2(12):  49.7 (p < 2.2e-16)
#> Root MSE:       0.3274 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(4) = 41.54 (p = 0.0000)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           13.79 
#>   Kleibergen-Paap rk Wald F:     12.17 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       24.58 
#>      15%  maximal IV size       13.96 
#>      20%  maximal IV size       10.26 
#>      25%  maximal IV size       8.31 
#>   Stock-Yogo critical values (IV relative bias):
#>      5%  maximal IV relative bias   16.85 
#>      10%  maximal IV relative bias   10.27 
#>      20%  maximal IV relative bias   6.71 
#>      30%  maximal IV relative bias   5.34 
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(3) = 74.16 (p = 0.0000)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(4,742) = 25.77 (p < 2.2e-16)
#>   Anderson-Rubin Wald Chi-sq(4) = 105.31 (p < 2.2e-16)
#>   Stock-Wright LM S Chi-sq(4) = 76.06 (p = 0.0000)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 0.45 (p = 0.5021)
#>   Tested: iq
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   iq                 12.17    0.0000      0.0692    0.0692     12.17     12.17
#> 
#> Instrumented:          iq 
#> Included instruments:  s, expr, tenure, rns, smsa, factor(year)67, factor(year)68, factor(year)69, factor(year)70, factor(year)71, factor(year)73 
#> Excluded instruments:  med, kww, age, mrt 
#> 
# }
```
