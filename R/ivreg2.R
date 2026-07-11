# --------------------------------------------------------------------------
# ivreg2 — the public orchestrator
#
# The estimation pipeline is one phase per file, called in this order:
#   R/validate-args.R         .validate_and_normalize_args
#   R/parse-formula.R         .parse_formula
#   R/prepare-model.R         .prepare_model
#   R/fit-dispatch.R          .resolve_and_fit
#   R/vcov-dispatch.R         .compute_vcov
#   R/diagnostics-dispatch.R  .compute_diagnostics
#   R/stored-matrices.R       .compute_stored_sw
#   R/ivreg2-class.R          .new_ivreg2 (constructor, print, summary)
# --------------------------------------------------------------------------

#' Extended Instrumental Variables Estimation
#'
#' Estimate models by OLS, two-stage least squares (2SLS), LIML, Fuller, or
#' k-class with automatic diagnostic tests. Uses a three-part formula for IV:
#' `y ~ exog | endo | instruments`.
#'
#' @param formula A formula: `y ~ exog` (OLS),
#'   `y ~ exog | endo | instruments` (IV), or
#'   `y ~ exog | 0 | instruments` (no endogenous regressors, with the
#'   excluded instruments contributing surplus moment conditions —
#'   Stata's `(=z1 z2)` form). The latter is estimated by OLS; the surplus
#'   conditions drive the Sargan/Hansen J overidentification test and
#'   `orthog` C-tests, and `method = "gmm2s"` with `vcov = "robust"` gives
#'   Cragg's (1983) heteroskedastic OLS (HOLS) estimator (with the default
#'   iid VCE the two-step weighting matrix is proportional to
#'   \eqn{(Z'Z)^{-1}} and the estimates equal OLS exactly).
#'   The response must be numeric; a factor or character response is rejected
#'   with an error.
#' @param data A data frame containing the variables in the formula.
#' @param weights Optional analytic weights expression (evaluated in `data`),
#'   equivalent to Stata's `[aw=varname]`. Must be strictly positive.
#'   For aweights and pweights, the weights are normalized internally to sum
#'   to N, following Stata's convention. This makes sigma (RMSE)
#'   scale-invariant for those types: multiplying all weights by a constant
#'   changes no coefficients, standard errors, or test statistics. Frequency
#'   and importance weights are *not* normalized --- they redefine N as
#'   `sum(weights)` (see `weight_type`), so their scale matters by
#'   construction.
#'   See [ivreg2r-conventions] for how weights interact with sigma and for
#'   the comparison to [lm()] (coefficients match `lm()` for OLS; SEs and the
#'   VCV match only when `small = TRUE`).
#' @param subset Optional subset expression (evaluated in `data`).
#' @param na.action Function for handling `NA`s (default [na.omit]).
#'   Use [na.exclude] to drop incomplete cases during estimation but
#'   reinsert `NA`s at the omitted positions in [fitted()], [residuals()],
#'   and [predict()] output, so results align with the original data frame.
#'   Stata drops missing observations silently (equivalent to [na.omit])
#'   and has no analogue of [na.exclude].
#' @param vcov Character: covariance type. One of `"iid"` (classical),
#'   `"robust"` (White heteroskedasticity-robust), `"HAC"`
#'   (heteroskedasticity and autocorrelation consistent), or `"AC"`
#'   (autocorrelation consistent). Matching is case-insensitive (`"hac"`,
#'   `"Hac"`, and `"HAC"` are all accepted), but the value is always
#'   normalized to and printed as the canonical spelling shown above. Use
#'   `small = TRUE` to apply finite-sample corrections. To match Stata's
#'   `ivreg2, robust`: use `vcov = "robust"`. To match Stata's
#'   `ivreg2, robust small`: use `vcov = "robust", small = TRUE`.
#' @param clusters One-sided formula specifying one or two cluster variables
#'   (e.g. `~ firmid` for one-way, `~ firmid + year` for two-way).
#'   Two-way clustering uses the Cameron-Gelbach-Miller (2006) formula.
#'   The effective cluster count is `min(M1, M2)` per Stata convention.
#'   The `small` argument controls whether the finite-sample correction
#'   `(N-1)/(N-K) * M/(M-1)` is applied (matching Stata's
#'   `cluster() small` combination).
#'   Unlike Stata's `ivreg2`, which silently drops observations with missing
#'   cluster values from the estimation sample, ivreg2r raises an error;
#'   the user should drop or handle those rows explicitly (e.g., by
#'   filtering them out) so the estimation sample is never changed silently.
#' @param endog Character vector of endogenous regressor names to test for
#'   exogeneity (endogeneity test / C-statistic). If `NULL` (default), tests
#'   all endogenous regressors. Names must match variables in the endogenous
#'   part of the formula. Ignored for OLS models.
#'
#'   **Note:** Unlike Stata's `ivreg2`, which only computes the endogeneity
#'   test when the `endog()` option is explicitly specified, `ivreg2r` computes
#'   it automatically for all IV models.
#' @param orthog Character vector of instrument names to test for
#'   orthogonality (instrument-subset C-statistic). Names must be included
#'   or excluded instruments (not endogenous regressors or the intercept).
#'   If `NULL` (default), no orthogonality test is computed. Ignored for
#'   OLS models. Equivalent to Stata's `orthog()` option.
#' @param redundant Character vector of excluded instrument names to test
#'   for redundancy (zero first-stage explanatory power). The test is a
#'   KP rk LM test of H0: rank=0 on the first-stage coefficient matrix
#'   for the tested instruments, conditional on maintained instruments.
#'   If `NULL` (default), no redundancy test is computed. On a model with no
#'   endogenous regressors (a one-part OLS formula, or the IV form
#'   `y ~ exog | 0 | instruments`), a warning reports that the test is
#'   skipped; Stata's `redundant()` is silently ignored in those cases.
#'   `noid = TRUE` and `b0` also suppress the redundancy test; supplying
#'   `redundant` alongside either is dropped with a warning rather than
#'   silently ignored.
#'   Equivalent to Stata's `redundant()` option.
#' @param small Logical: if `TRUE`, use small-sample corrections
#'   (t/F instead of z/chi-squared, `N-K` denominator for sigma).
#' @param weight_type Character: type of weights. One of `"aweight"`
#'   (analytic weights, default), `"fweight"` (frequency weights),
#'   `"pweight"` (probability/sampling weights), or `"iweight"`
#'   (importance weights).
#'
#'   **aweight**: Normalized to sum to N. Standard for WLS.
#'
#'   **fweight**: Integer-valued; N is redefined as `sum(weights)`.
#'   HC meat uses `w * e^2` (linear, not quadratic).
#'
#'   **pweight**: Normalized to sum to N. Forces robust VCE, with a warning
#'   (overrides `vcov = "iid"` to `"robust"` and `vcov = "AC"` ---
#'   explicit or kiefer-implied --- to `"HAC"`).
#'
#'   **iweight**: Not normalized; N is redefined as `sum(weights)` (float).
#'   All intermediate calculations use float N; posted `nobs` and `df.residual`
#'   use `floor(N)`. Restricted to IID VCE only — incompatible with robust,
#'   cluster, HAC, and GMM estimation.
#' @param dofminus Non-negative integer: large-sample degrees-of-freedom
#'   adjustment. Subtracted from N in large-sample variance formulas
#'   (e.g., sigma = rss/(N-dofminus)). Useful when fixed effects have been
#'   partialled out. Equivalent to Stata's `dofminus()` option.
#' @param sdofminus Non-negative integer: small-sample degrees-of-freedom
#'   adjustment. Subtracted from the residual degrees of freedom alongside K
#'   (e.g., df.residual = N - K - dofminus - sdofminus). Useful when
#'   partialling out regressors. Equivalent to Stata's `sdofminus()` option.
#' @param method Character: estimation method. One of `"2sls"` (default),
#'   `"liml"`, `"kclass"`, `"gmm2s"` (two-step efficient GMM), or `"cue"`
#'   (continuously updated GMM estimator).
#'   For OLS models (1-part formula), this is ignored. When `fuller > 0`
#'   is specified, method is automatically promoted to `"liml"`.
#'   `"gmm2s"` uses the inverse of the moment covariance matrix as the
#'   optimal weighting matrix, yielding efficient estimates under the
#'   specified error structure. Incompatible with `fuller` and `kclass`.
#'   `"cue"` continuously updates both moment conditions and moment
#'   covariance at each candidate beta during optimization. Under iid errors
#'   (no robust/cluster/kernel VCE), CUE with no other options is equivalent
#'   to LIML with `coviv = TRUE`. Under non-iid errors, CUE generalizes
#'   LIML to produce efficient estimates.
#'   Incompatible with `fuller`, `kclass`, `wmatrix`, and `smatrix`.
#'
#'   **CUE optimizer note:** This package uses R's `optim()` (BFGS +
#'   Nelder-Mead, scaled by the starting values) while Stata uses Mata's
#'   `optimize()` (modified Newton-Raphson). The CUE objective is
#'   non-convex and can have multiple local minima, and it can possess
#'   economically degenerate minima (e.g., negative R-squared) -- a
#'   documented property of the CUE ratio objective (Hausman, Lewis,
#'   Menzel & Newey, 2011), not an optimizer failure. On every benchmarked
#'   model the two implementations agree on the minimum, including one such
#'   pathological case, where both converge to the same degenerate basin;
#'   because those basins are extremely flat, coefficient agreement across
#'   implementations there is limited to a few digits. Inspect the
#'   J-statistic and R-squared to judge whether a CUE solution is
#'   economically sensible.
#' @param kclass Numeric scalar: user-supplied k value for k-class
#'   estimation. When supplied, `method` is automatically set to `"kclass"`.
#'   Must be non-negative. Cannot be combined with `method = "liml"` or
#'   `fuller`.
#' @param fuller Numeric scalar: Fuller (1977) modification parameter.
#'   Non-negative. `0` (the default) applies no Fuller adjustment; a
#'   positive value activates the Fuller-modified LIML estimator, which
#'   automatically sets `method` to `"liml"` and
#'   `k = lambda - fuller / (N - L)`. `fuller = 1` gives the bias-corrected
#'   LIML estimator; `fuller = 4` targets MSE. Cannot be combined with
#'   `kclass`.
#' @param coviv Logical: if `TRUE`, use the 2SLS bread `(X_hat'X_hat)^{-1}`
#'   instead of the k-class bread for VCV computation in LIML/k-class
#'   estimation. This gives the "COVIV" (covariance at the IV estimates)
#'   VCV that is robust to misspecification of the LIML model. Applies only
#'   to `method = "liml"` and `method = "kclass"`; for every other method
#'   (OLS, 2SLS, GMM2S, CUE) it is reset to `FALSE` with a warning.
#'   Default `FALSE`.
#' @param kernel Character: kernel function for HAC/AC standard errors.
#'   One of `"bartlett"`, `"parzen"`, `"truncated"`, `"tukey-hanning"`,
#'   `"tukey-hamming"`, `"qs"` (quadratic spectral), `"daniell"`, or
#'   `"tent"`. When specified without an explicit `vcov` change, `vcov`
#'   is automatically set: `"iid"` becomes `"AC"`, `"robust"` becomes
#'   `"HAC"`. Requires `bw` and `tvar`.
#'
#'   **Note:** Kleibergen-Paap identification tests (underidentification and
#'   weak identification) always use the Bartlett kernel internally, regardless
#'   of the user-specified kernel. This matches a known behavior in Stata's
#'   `ivreg2`, where `ranktest` hard-codes the Bartlett kernel.
#' @param bw Numeric or `"auto"`: bandwidth for kernel estimation. Must be
#'   positive numeric, or `"auto"` for automatic selection via Newey-West
#'   (1994). Auto selection is available for Bartlett, Parzen, and Quadratic
#'   Spectral kernels only, and is not supported for panel data (`ivar`).
#'   Required when `kernel` is specified.
#' @param tvar Character: name of the time variable in `data`. Required
#'   for HAC/AC estimation and for time-series operators `l()`/`d()` in the
#'   formula (see [`ts-operators`][l]). Lags are resolved by time *value*
#'   (gap-aware), mirroring Stata's `tsset`.
#' @param ivar Character: name of the panel identifier variable in `data`.
#'   Optional; only needed for panel data (HAC/AC estimation, or so that
#'   time-series operators lag within panel units).
#' @param kiefer Logical: if `TRUE`, use the Kiefer (1980) VCE —
#'   autocorrelation-consistent with kernel = Truncated and bandwidth = T
#'   (the full time span). Requires panel data (`tvar` + `ivar`).
#'   Incompatible with robust VCE, clustering, explicit kernel or bandwidth.
#'   Equivalent to Stata's `ivreg2 ..., kiefer`.
#' @param dkraay Positive numeric scalar: bandwidth for Driscoll-Kraay (1998)
#'   VCE. When specified, clusters on the time variable and applies kernel
#'   smoothing across time lags, producing standard errors robust to
#'   cross-sectional dependence. Requires panel data (`tvar` + `ivar`).
#'   Incompatible with explicit `bw`. If `kernel` is not specified, defaults
#'   to Bartlett. Equivalent to Stata's `ivreg2 ..., dkraay(3)`.
#' @param wmatrix Numeric matrix: user-supplied L x L weighting matrix for
#'   GMM estimation, where L is the number of instruments (including exogenous
#'   regressors). When supplied without `method = "gmm2s"`, produces
#'   inefficient GMM with full sandwich VCV (`method` becomes `"gmmw"`).
#'   When combined with `method = "gmm2s"`, used as the first-step weighting
#'   matrix (second step uses the optimal weighting matrix). Must be symmetric.
#'   Incompatible with `method = "liml"`, `kclass`, and `fuller`.
#'   When `method` is not `"gmm2s"`, ignored without robust VCE, clustering,
#'   or HAC (with a warning).
#'   Equivalent to Stata's `wmatrix()` option. When used (not ignored), it is
#'   echoed back as `fit$W`, matching Stata's `e(W)`. A matrix with dimnames
#'   is matched to the instrument columns by name (reordered as needed,
#'   mirroring Stata's matsort; names that do not cover the instrument list
#'   are an error); an unnamed matrix is used positionally, in instrument
#'   column order.
#' @param smatrix Numeric matrix: user-supplied L x L moment covariance matrix
#'   for GMM estimation. When supplied, the GMM estimation uses this matrix
#'   instead of computing Omega from residuals. Implies efficient GMM
#'   (`method` is promoted to `"gmm2s"` if `"2sls"`). Must be symmetric.
#'   Incompatible with `method = "liml"`, `kclass`, and `fuller`.
#'   Equivalent to Stata's `smatrix()` option. Echoed back as `fit$S`
#'   (never recomputed), matching Stata's `e(S)`. A matrix with dimnames is
#'   matched to the instrument columns by name (reordered as needed,
#'   mirroring Stata's matsort); an unnamed matrix is used positionally, in
#'   instrument column order.
#' @param b0 Numeric vector: evaluate the CUE objective at this fixed
#'   parameter vector without optimization. When supplied, `method` is
#'   promoted to `"cue"` and the J(b0) statistic is computed and stored.
#'   Identification diagnostics (underid, weak-id, first-stage, AR, SW,
#'   endogeneity, orthog, and the redundancy test) are suppressed (matching
#'   Stata's `b0()` option); `reduced_form` is not suppressed. An explicit
#'   `endog`, `orthog`, `redundant`, or `first_stage = TRUE` request supplied
#'   alongside `b0` is dropped with a warning, rather than silently ignored.
#'   Length must equal the number of regressors K. Named vectors are
#'   reordered to match model matrix columns; unnamed vectors are used
#'   as-is in model matrix column order.
#'   Incompatible with `method = "gmm2s"`, `"liml"`, `kclass`, `fuller`,
#'   and `wmatrix`.
#' @param noid Logical: if `TRUE`, suppress computation of underidentification,
#'   weak identification, and redundancy test statistics. Overidentification,
#'   Anderson-Rubin, Stock-Wright, endogeneity, and first-stage diagnostics
#'   are still computed. Default `FALSE`. Useful for speeding up simulation
#'   loops where rank-based identification tests are expensive.
#'   An explicit `redundant` request supplied alongside `noid = TRUE` is
#'   dropped with a warning, rather than silently ignored.
#'   Equivalent to Stata's `noid` option.
#' @param partial Character vector: exogenous regressors to partial out
#'   via Frisch-Waugh-Lovell projection before estimation. Coefficients on
#'   partialled variables are not recoverable and are not reported.
#'   Special values:
#'   - `"_cons"` or `"(Intercept)"`: partial only the constant (demean).
#'   - `"_all"`: partial all included exogenous regressors.
#'   - `NULL` (default): no partialling.
#'
#'   By default, `sdofminus` is automatically incremented by the number of
#'   partialled variables (including the constant). Use `nopartialsmall = TRUE`
#'   to suppress this adjustment.
#'
#'   **Note:** FWL invariance holds for OLS, 2SLS, LIML, and two-step GMM, but
#'   **not** for CUE. The CUE objective recomputes the moment covariance at each
#'   iteration, so partialled residuals produce a different optimization surface.
#'   A warning is issued when `partial` is used with `method = "cue"`.
#'
#'   After partialling, `predict()` is restricted to residuals only (no
#'   `newdata`). Summary output notes that total SS, model F, and R-squared
#'   are partial-model values.
#'
#'   Equivalent to Stata's `partial()` option.
#' @param nopartialsmall Logical: if `TRUE`, suppress the automatic
#'   `sdofminus` adjustment from partialling. Default `FALSE`.
#'   Equivalent to Stata's `nopartialsmall` option.
#' @param center Logical: if `TRUE`, subtract the mean of moment conditions
#'   (scores) before computing the S matrix (meat of the sandwich VCE).
#'   Default `FALSE`. Centering only affects non-homoskedastic VCE types
#'   (HC, cluster, HAC); a warning is issued if used with IID or AC VCE.
#'
#'   **Note:** Centering is applied to the main model's VCE and to
#'   diagnostic tests that use the main model's S matrix (overidentification,
#'   orthogonality). However, the endogeneity test (`endog`) computes its own
#'   S matrix from the restricted model *without* centering, even when
#'   `center = TRUE`. This matches Stata's `ivreg2`, where `center` is not
#'   forwarded to the recursive call for the endogeneity test.
#'
#'   The Kleibergen-Paap identification statistics (underidentification,
#'   weak identification) and the instrument redundancy test are likewise
#'   always computed *without* centering, regardless of `center`. This
#'   matches Stata, whose `ranktest` calls (underid, weak-id, and
#'   redundancy) never receive the `center` option.
#'
#'   When a user-supplied `smatrix` is given, `center` has no effect on the
#'   estimator, the GMM weighting matrix, the stored `$S`, or any statistic
#'   computed from the supplied matrix (the overidentification and
#'   orthogonality tests reuse it): the supplied matrix is used as-is
#'   (mirroring the `psd` argument's identical interaction with `smatrix`,
#'   documented below). Centering still applies only to diagnostics that
#'   rebuild their own moment covariance, such as the Stock-Wright S
#'   statistic.
#' @param psd Character or NULL: PSD correction for the moment covariance
#'   matrix S. `NULL` (default) applies no correction. `"psd0"` zeroes
#'   negative eigenvalues (Politis 2007). `"psda"` replaces negative
#'   eigenvalues with their absolute values (Stock & Watson 2008). The
#'   correction is applied to S *before* the variance-covariance matrix is
#'   assembled from it, matching Stata's `m_omega`; the corrected S is
#'   stored as `$S`. The conventional (`vcov = "iid"`) VCV and the
#'   Kleibergen-Paap identification statistics are never corrected, also
#'   matching Stata (its `ranktest` does not receive the psd option). A
#'   warning is emitted whenever negative eigenvalues are detected and
#'   corrected (Stata corrects silently). Ignored when a user-supplied
#'   `smatrix` is given: the user matrix is used as-is for estimation, the
#'   VCV, and `$S`, matching Stata (whose `m_omega` is never invoked for a
#'   supplied S). Equivalent to Stata's `psd0` and `psda` options.
#'
#'   Under an exactly-singular `psd0`-corrected VCV, the model F statistic
#'   uses a conditioning-guarded (swept) inverse and can differ from Stata,
#'   which inverts the near-singular block directly; neither is a
#'   well-defined Wald statistic in that degenerate case.
#' @param sw Logical: if `TRUE`, compute the Stock-Watson (2008,
#'   Econometrica) panel-robust VCE. Requires panel data (`ivar`).
#'   Incompatible with clustering, HAC kernels, kiefer, dkraay,
#'   fweights, and iweights. Forces `vcov = "robust"` automatically, with a
#'   warning when `vcov = "iid"` is thereby overridden.
#'   Equivalent to Stata's `sw` option (labeled "BETA VERSION" in Stata).
#'
#'   **Precondition:** the Stock-Watson (2008) correction is derived for a
#'   within-transformed (fixed-effects) panel regression with fixed `T` and
#'   large `N`. `ivreg2()` does not perform the within transformation for
#'   you: within-transform the data yourself (demean every variable by
#'   panel unit) before fitting, and pass `dofminus` equal to the number of
#'   panel units, which charges the absorbed unit means to the degrees of
#'   freedom. See the worked example in the "Fixed-effects panels: the
#'   Stock-Watson correction" section of `vignette("time-series-gmm")`.
#' @param reduced_form Character: what reduced-form output to store.
#'   `"none"` (default) stores nothing. `"rf"` stores the y ~ Z regression
#'   (equivalent to Stata's `saverf`; Stata's `rf` option *displays* the
#'   reduced form without storing it). `"system"` stores the full system of
#'   y + all endogenous variables regressed on Z, with cross-equation VCV
#'   (equivalent to Stata's `savesfirst`, an option present in `ivreg2.ado`
#'   but absent from its help file; Stata's `sfirst` displays the system
#'   without storing it). On a model with no endogenous regressors (a
#'   one-part OLS formula, or the IV form `y ~ exog | 0 | instruments`), no
#'   reduced form is computed — matching Stata, which skips reduced-form
#'   estimation whenever the endogenous list is empty — and a warning
#'   reports that the request was ignored.
#'
#'   **Note:** the reduced-form variance-covariance matrix always applies
#'   OLS-style small-sample degrees-of-freedom scaling, regardless of the
#'   main model's `small` argument. This matches Stata's `ivreg2`
#'   (ivreg2.ado:3091-3097). See [ivreg2r-conventions] for the full statement
#'   of which corrections `small` controls.
#' @param first_stage Logical: if `TRUE`, store extractable first-stage
#'   regression objects on the fitted model. Access them via
#'   [first_stage()]; see [first_stage()] for the full list of supported
#'   S3 methods (`coef()`, `vcov()`, `summary()`, `tidy()`, `glance()`,
#'   and others). Equivalent to Stata's `savefirst` option. Default `FALSE`. On a model with no endogenous
#'   regressors (a one-part OLS formula, or the IV form
#'   `y ~ exog | 0 | instruments`), nothing is stored and a warning reports
#'   that the request was ignored.
#' @param model Logical: if `TRUE` (default), store the model frame in the
#'   return object.
#' @param x Logical: if `TRUE`, store model matrices (`X`, `Z`) in the
#'   return object.
#' @param y Logical: if `TRUE` (default), store the response vector in the
#'   return object.
#'
#' @return An object of class `"ivreg2"`. Beyond the usual components
#'   (`coefficients`, `residuals`, `vcov`, `diagnostics`, ...), the object
#'   stores the estimated moment-condition covariance matrix as `$S` and the
#'   GMM weighting matrix as `$W`, corresponding to Stata's `e(S)` and
#'   `e(W)`:
#'
#'   - `$S` is present for every fit. Normalization matches Stata's
#'     `m_omega`: iid `sigma^2 (Z'Z)/N` with `sigma^2 = RSS/(N - dofminus)`;
#'     robust and HAC, meat divided by `N - dofminus`; cluster, meat divided
#'     by `N`; psd-corrected when `psd` is set; built from centered moments
#'     when `center = TRUE`. A user-supplied `smatrix` is echoed, never
#'     recomputed and never psd-corrected.
#'   - `$W` is `NULL` for `method = "liml"` and `"kclass"` (as in Stata, no
#'     weighting matrix is defined for k-class estimators). For two-step
#'     GMM, CUE, and `b0` fits it is the inverse of the defining `$S`; for
#'     one-step OLS/2SLS it is `(1/sigma^2) (Z'Z/N)^{-1}`; a user-supplied
#'     `wmatrix` is echoed (it is the first-step weighting matrix under
#'     `method = "gmm2s"`).
#'
#'   Rows and columns are named by the R instrument columns — `(Intercept)`
#'   rather than Stata's `_cons`, and R column order (exogenous regressors
#'   before excluded instruments) rather than Stata's. Passing `fit$S` back
#'   via `smatrix` reproduces the corresponding efficient-GMM fit, mirroring
#'   the `e(S)` reuse workflow in Stata's `ivreg2` help file.
#'
#' @section Overidentification and endogeneity tests:
#'
#' When the model carries more instruments than regressors, `ivreg2()`
#' automatically reports an overidentification test: the Sargan (1958)
#' statistic under `vcov = "iid"`, and the Hansen J statistic under robust,
#' cluster, or HAC errors. The test evaluates the `L - K` overidentifying
#' restrictions jointly, where `L` is the number of instruments and `K` the
#' number of regressors. A rejection indicates that at least one instrument
#' violates its orthogonality condition, without identifying which one. A
#' non-rejection is necessary but not sufficient for instrument validity.
#' The test has no power against a bias shared by every instrument, so
#' orthogonality conditions that fail in the same direction can still appear
#' jointly satisfied. The test is also uninformative about the `K`
#' exactly-identifying moment conditions, whose validity must be defended on
#' substantive grounds. An exactly identified model has no overidentifying
#' restrictions, and therefore no overidentification test.
#'
#' The endogeneity test, also called the C-statistic or difference test,
#' asks whether regressors treated as endogenous could instead be treated
#' as exogenous. It is the overidentification statistic of the model that
#' adds the suspect regressors to the instrument set, treating them as
#' exogenous, minus that of the model that treats them as endogenous. Under
#' the null of exogeneity, this difference follows a chi-square distribution
#' with degrees of freedom equal to the number of regressors tested. Both
#' overidentification statistics are formed from the same estimated moment
#' covariance, which guarantees that the difference is non-negative
#' (Hayashi 2000).
#' `ivreg2()` computes this test automatically for every IV model, whereas
#' Stata computes it only when the `endog()` option is given explicitly.
#' The `endog` argument selects which endogenous regressors are tested. The
#' statistic takes the difference-in-Sargan form under iid errors and the
#' difference-in-Hansen-J form under robust, cluster, or HAC errors.
#'
#' @section LIML, Fuller, and k-class estimators:
#'
#' LIML, Fuller's modification, and the general k-class estimators are
#' alternatives to 2SLS motivated by weak instruments. LIML is approximately
#' median-unbiased, and its Stock-Yogo weak-identification thresholds are
#' more forgiving than those for 2SLS. LIML has no finite-sample moments of
#' any order, however. Its sampling distribution can therefore be dispersed,
#' and individual estimates can lie far from the truth. Fuller's (1977)
#' modification sets `k = lambda - fuller / (N - L)`, which yields an
#' estimator with finite moments. The choice `fuller = 1` is approximately
#' best-unbiased, and `fuller = 4` minimizes mean squared error, both within
#' the standard weak-instrument framework. Weak-identification pretesting for
#' LIML and Fuller uses estimator-specific critical values rather than the
#' 2SLS values. LIML coincides with 2SLS when the model is exactly
#' identified, because the LIML eigenvalue equals one in that case. Fuller
#' does not, since it subtracts `fuller / (N - L)` from that eigenvalue.
#'
#' @seealso [ivreg2r-conventions] for the statistical conventions (sigma
#'   normalization, `small`, weights, degrees of freedom) and how they differ
#'   from [lm()]; [ivreg2r-glossary] for the diagnostic-output acronyms;
#'   [summary.ivreg2()] and [print.ivreg2()] for console output;
#'   [tidy.ivreg2()], [glance.ivreg2()], and [augment.ivreg2()] for
#'   broom-style output; [first_stage()] for first-stage regressions.
#'   The vignettes `vignette("introduction", "ivreg2r")`,
#'   `vignette("advanced-iv", "ivreg2r")`, and
#'   `vignette("time-series-gmm", "ivreg2r")` work through applied examples.
#'
#' @family ivreg2 methods
#'
#' @examples
#' data(mroz)
#' mroz_work <- subset(mroz, inlf == 1)
#'
#' # --- Just-identified IV: return to schooling ---
#' # Wooldridge (2020), Example 15.1: father's education instruments
#' # education in a simple wage equation (replicates the published
#' # estimates: educ = 0.059, SE = 0.035). With one instrument for one
#' # endogenous variable, no overid test is possible; instrument validity
#' # must be defended on substantive grounds.
#' fit <- ivreg2(lwage ~ 1 | educ | fatheduc,
#'               data = mroz_work)
#' summary(fit)
#'
#' # --- Overidentified IV: testing instrument validity ---
#' # The Stata ivreg2 help file's illustrative baseline (line 1274):
#' # age, kidslt6, and kidsge6 give two overidentifying restrictions,
#' # so the Sargan test can detect misspecification. Note its weak
#' # first stage (Cragg-Donald F ~ 4.3) -- see the LIML example below.
#' fit_overid <- ivreg2(lwage ~ exper + expersq | educ |
#'                        age + kidslt6 + kidsge6, data = mroz_work)
#' summary(fit_overid)
#'
#' # --- Robust standard errors ---
#' # With heteroskedasticity, Kleibergen-Paap diagnostics replace
#' # Anderson/Cragg-Donald, and Hansen J replaces Sargan.
#' fit_robust <- ivreg2(lwage ~ exper + expersq | educ |
#'                        age + kidslt6 + kidsge6, data = mroz_work,
#'                        vcov = "robust")
#' summary(fit_robust)
#'
#' \donttest{
#' # --- LIML ---
#' # Same baseline as above (a documented variation on the help-file
#' # spec). Its weak first stage (Cragg-Donald F ~ 4.3) is the setting
#' # where 2SLS bias toward OLS grows with the overidentification degree
#' # relative to instrument strength, and LIML's approximate
#' # median-unbiasedness -- established within the iid Stock-Yogo
#' # framework -- is attractive. LIML equals 2SLS when just-identified.
#' fit_liml <- ivreg2(lwage ~ exper + expersq | educ |
#'                      age + kidslt6 + kidsge6, data = mroz_work,
#'                      method = "liml")
#' summary(fit_liml)
#'
#' # --- Endogeneity test ---
#' # Is education actually endogenous? The C-statistic tests
#' # H0: OLS is consistent (educ is exogenous).
#' # ivreg2 help file line 1283: the endog(educ) example on this equation.
#' fit_endog <- ivreg2(lwage ~ exper + expersq | educ |
#'                       age + kidslt6 + kidsge6, data = mroz_work,
#'                       endog = "educ")
#' if (requireNamespace("dplyr", quietly = TRUE)) {
#'   diagnostics(fit_endog) |> dplyr::filter(test == "endogeneity")
#' }
#'
#' # --- Clustering ---
#' # Griliches (1976) wage equation, cluster on year.
#' # Adapted from Stata `ivreg2` help-file example (line 1250): the help
#' # command also includes year dummies (xi i.year), omitted here, and no
#' # small option. Either way the fit warns -- 7 year clusters < the
#' # instrument count makes the moment covariance singular, which is the
#' # help file's own lead-in to partialling; see
#' # vignette("time-series-gmm").
#' data(griliches)
#' fit_cl <- ivreg2(lw ~ s + expr + tenure + rns + smsa |
#'                    iq | med + kww + age,
#'                  data = griliches, clusters = ~year, small = TRUE)
#' summary(fit_cl)
#'
#' # --- Two-step efficient GMM ---
#' # ivreg2 help file line 1154: efficient GMM with robust weighting.
#' # The optimal weighting matrix uses the first-step robust moment
#' # covariance; Hansen J replaces Sargan.
#' fit_gmm <- ivreg2(lw ~ s + expr + tenure + rns + smsa + factor(year) |
#'                     iq | med + kww + age + mrt,
#'                   data = griliches, method = "gmm2s", vcov = "robust")
#' summary(fit_gmm)
#' }
#'
#' @export
ivreg2 <- function(formula, data, weights, subset, na.action = stats::na.omit,
                   vcov = "iid", clusters = NULL, endog = NULL,
                   orthog = NULL,
                   redundant = NULL,
                   method = "2sls", kclass = NULL, fuller = 0,
                   coviv = FALSE,
                   small = FALSE,
                   dofminus = 0L, sdofminus = 0L,
                   weight_type = "aweight",
                   kernel = NULL, bw = NULL, tvar = NULL, ivar = NULL,
                   kiefer = FALSE, dkraay = NULL,
                   wmatrix = NULL, smatrix = NULL,
                   b0 = NULL,
                   noid = FALSE,
                   partial = NULL,
                   nopartialsmall = FALSE,
                   center = FALSE,
                   psd = NULL,
                   sw = FALSE,
                   reduced_form = "none",
                   first_stage = FALSE,
                   model = TRUE, x = FALSE, y = TRUE) {

  # --- 1. Capture call ---
  cl <- match.call()

  # --- 2. Validate and normalize arguments ---
  opts <- .validate_and_normalize_args(
    vcov = vcov, clusters = clusters, endog = endog, orthog = orthog,
    redundant = redundant, method = method, kclass = kclass, fuller = fuller,
    coviv = coviv, small = small, dofminus = dofminus, sdofminus = sdofminus,
    weight_type = weight_type, kernel = kernel, bw = bw, tvar = tvar,
    ivar = ivar, kiefer = kiefer, dkraay = dkraay,
    wmatrix = wmatrix, smatrix = smatrix, b0 = b0,
    partial = partial, nopartialsmall = nopartialsmall,
    center = center, psd = psd, reduced_form = reduced_form,
    first_stage = first_stage, noid = noid, sw = sw,
    model_flag = model, x_flag = x, y_flag = y
  )
  # Unpack options needed in ivreg2() scope (for .new_ivreg2() assembly).
  # Helpers read from `opts` directly.
  small       <- opts$small
  dofminus    <- opts$dofminus
  coviv       <- opts$coviv
  weight_type <- opts$weight_type
  kernel      <- opts$kernel
  kiefer      <- opts$kiefer
  dkraay      <- opts$dkraay
  tvar        <- opts$tvar
  ivar        <- opts$ivar
  psd         <- opts$psd
  smatrix     <- opts$smatrix

  # --- 3. Forward to parser ---
  # Build a call to .parse_formula() using the NSE arguments from ivreg2().
  # Evaluate in an environment where .parse_formula is visible (it's unexported)
  # with the user's frame as parent so data/subset/weights expressions resolve.
  pf_call <- cl[c(1L, match(c("formula", "data", "weights", "subset",
                               "na.action"), names(cl), 0L))]
  pf_call[[1L]] <- quote(.parse_formula)
  # Time/panel variable names for ts-operator resolution (literal strings,
  # not user expressions, so direct assignment into the call is safe).
  pf_call$tvar <- opts$tvar
  pf_call$ivar <- opts$ivar
  pf_env <- list2env(list(.parse_formula = .parse_formula),
                     parent = parent.frame())
  parsed <- eval(pf_call, pf_env)

  # --- 3a. Prepare model ---
  prep <- .prepare_model(parsed, opts, data, clusters)
  parsed       <- prep$parsed
  sdofminus    <- prep$sdofminus
  b0           <- prep$b0
  # NOTE: `smatrix` stays the RAW user matrix for the backward-compat
  # fit$smatrix echo field; the name-normalized copy (prep$smatrix) is what
  # estimation and fit$S consume.
  cluster_vec  <- prep$cluster_vec
  cluster_var_name <- prep$cluster_var_name
  M            <- prep$M
  M1           <- prep$M1
  M2           <- prep$M2
  time_index   <- prep$time_index
  unsort_order <- prep$unsort_order
  partial_ct   <- prep$partial_ct
  partial_names <- prep$partial_names
  partialcons  <- prep$partialcons
  n_physical   <- prep$n_physical
  w_raw        <- prep$w_raw
  bw           <- prep$bw

  # --- 4. Resolve estimation path and fit ---
  est <- .resolve_and_fit(parsed, prep, opts)
  fit      <- est$fit
  method   <- est$method
  bw       <- est$bw
  # Backward-compat fit$wmatrix echo keeps the RAW user matrix (NULL when
  # the wmatrix was ignored-with-warning); the name-normalized copy
  # (est$wmatrix) feeds estimation and fit$W.
  wmatrix  <- if (is.null(est$wmatrix)) NULL else opts$wmatrix
  omega_fn <- est$omega_fn

  # --- 5. VCV ---
  fit <- .compute_vcov(fit, parsed, prep, opts, method, bw)

  # --- 5b. Diagnostics, reduced-form, model F ---
  diag_result <- .compute_diagnostics(fit, parsed, prep, opts, method, bw)
  diagnostics         <- diag_result$diagnostics
  first_stage         <- diag_result$first_stage
  first_stage_models  <- diag_result$first_stage_models
  reduced_form_result <- diag_result$reduced_form_result
  model_f_result      <- diag_result$model_f_result
  effective_vcov_type <- diag_result$effective_vcov_type
  center              <- diag_result$center

  # --- 5c. Stored moment-covariance S and weighting matrix W ---
  # Must precede the unsort: HAC/AC moment covariances index the residuals
  # via the sorted prep$time_index. Uses the effective center from the
  # diagnostics step (reset for iid/AC, matching m_omega).
  stored_sw <- .compute_stored_sw(fit, parsed, prep, opts, est,
                                  bw = bw, center = center)

  # --- 5d. Unsort for user-facing output ---
  # If data was sorted for HAC/AC, restore original row order for user-facing
  # vectors and matrices (residuals, fitted values, y, X, Z).
  if (!is.null(unsort_order)) {
    fit$residuals <- fit$residuals[unsort_order]
    fit$fitted.values <- fit$fitted.values[unsort_order]
    parsed$y <- parsed$y[unsort_order]
    parsed$X <- parsed$X[unsort_order, , drop = FALSE]
    parsed$Z <- parsed$Z[unsort_order, , drop = FALSE]
    # Also unsort first-stage model residuals/fitted values
    if (!is.null(first_stage_models)) {
      for (nm in names(first_stage_models)) {
        first_stage_models[[nm]]$residuals <-
          first_stage_models[[nm]]$residuals[unsort_order]
        first_stage_models[[nm]]$fitted.values <-
          first_stage_models[[nm]]$fitted.values[unsort_order]
      }
    }
  }

  # --- 5e. Compute additional stored results (R2) ---
  # Condition numbers: Stata's cond(XX) = max(eigenvalue)/min(eigenvalue)
  # which equals kappa(XX, exact = TRUE) in R.
  if (is.null(parsed$weights)) {
    XX <- crossprod(parsed$X)
  } else {
    XX <- crossprod(sqrt(parsed$weights) * parsed$X)
  }
  condxx <- kappa(XX, exact = TRUE)

  if (parsed$is_iv) {
    if (is.null(parsed$weights)) {
      ZZ <- crossprod(parsed$Z)
    } else {
      ZZ <- crossprod(sqrt(parsed$weights) * parsed$Z)
    }
    condzz <- kappa(ZZ, exact = TRUE)
  } else {
    # OLS: Z = X, so condzz = condxx (matches Stata behavior)
    condzz <- condxx
  }

  # Gaussian log-likelihood (Stata line 2009)
  ll <- -0.5 * (parsed$N * log(2 * pi) + parsed$N * log(fit$rss / parsed$N) +
                 parsed$N)

  # --- 6. Assemble return object ---
  # Determine effective method for the return object
  est_method <- if (method == "gmmw") {
    "gmmw"
  } else if (method %in% c("liml", "kclass", "gmm2s", "cue")) {
    method
  } else if (parsed$is_iv && parsed$K1 > 0L) {
    "2sls"
  } else {
    # Includes the empty-endogenous (K1 = 0) form: Stata posts
    # e(model) = "ols" (ivreg2.ado:2101-2106)
    "ols"
  }

  .new_ivreg2(
    coefficients  = fit$coefficients,
    residuals     = fit$residuals,
    fitted.values = fit$fitted.values,
    vcov          = fit$vcov,
    sigma         = fit$sigma,
    df.residual   = fit$df.residual,
    rank          = fit$rank,
    r.squared     = fit$r.squared,
    adj.r.squared = fit$adj.r.squared,
    rss           = fit$rss,
    r2u           = fit$r2u,
    r2c           = fit$r2c,
    mss           = fit$mss,
    model_f       = model_f_result$model_f,
    model_f_p     = model_f_result$model_f_p,
    model_f_df1   = model_f_result$model_f_df1,
    model_f_df2   = model_f_result$model_f_df2,
    diagnostics   = diagnostics,
    first_stage   = first_stage,
    first_stage_models = first_stage_models,
    reduced_form  = reduced_form_result,
    call          = cl,
    formula       = parsed$formula,
    terms         = parsed$terms,
    nobs          = if (weight_type == "iweight") as.integer(floor(parsed$N)) else parsed$N,
    vcov_type     = effective_vcov_type,
    small         = small,
    dofminus      = dofminus,
    sdofminus     = sdofminus,
    cluster_var   = cluster_var_name,
    n_clusters    = M,
    n_clusters1   = M1,
    n_clusters2   = M2,
    na.action     = parsed$na.action,
    weights       = w_raw,
    weight_type   = weight_type,
    n_physical    = if (weight_type %in% c("fweight", "iweight")) n_physical else NULL,
    endogenous    = parsed$endo_names,
    endo_colnames = parsed$endo_colnames,
    instruments   = parsed$excluded_names,
    dropped_regressors      = parsed$dropped_regressors,
    dropped_instruments     = parsed$dropped_instruments,
    reclassified_endogenous = parsed$reclassified_endogenous,
    contrasts     = parsed$contrasts,
    xlevels       = parsed$xlevels,
    method            = est_method,
    lambda            = fit$lambda %||% NA_real_,
    kclass_value      = fit$kclass_value %||% NA_real_,
    fuller_parameter  = fit$fuller_param %||% 0,
    coviv             = coviv,
    S                 = stored_sw$S,
    W                 = stored_sw$W,
    wmatrix           = wmatrix,
    smatrix           = smatrix,
    b0                = b0,
    noid              = noid,
    cue_convergence   = fit$convergence,   # NULL for non-CUE
    cue_message       = fit$cue_message,   # NULL for non-CUE
    kernel            = kernel,
    bw                = bw,
    tvar              = tvar,
    kiefer            = kiefer,
    dkraay            = dkraay,
    ivar              = ivar,
    has_ts_ops        = isTRUE(parsed$has_ts_ops),
    center            = center,
    psd               = psd,
    sw                = isTRUE(opts$sw),
    partial_ct        = partial_ct,
    partial_names     = partial_names,
    partialcons       = partialcons,
    yy                = fit$tss_u,
    yyc               = fit$tss_c,
    rankzz            = if (parsed$is_iv) parsed$L else fit$rank,
    condxx            = condxx,
    condzz            = condzz,
    ll                = ll,
    model         = if (model) parsed$model_frame else NULL,
    x             = if (x) list(X = parsed$X, Z = parsed$Z) else NULL,
    y             = if (y) parsed$y else NULL
  )
}
