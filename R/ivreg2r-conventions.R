#' Statistical conventions in ivreg2r
#'
#' `ivreg2r` is a translation of Stata's `ivreg2`, and it follows Stata's
#' statistical conventions rather than those of R's [lm()]/[vcov()] ecosystem.
#' This topic collects those conventions in one place. The API itself is
#' idiomatic R (formula interface, tidyverse-friendly output); only the
#' *statistical* choices below lean on Stata. See [ivreg2()] for the
#' estimation function and [ivreg2r-glossary] for the diagnostic acronyms
#' printed by [summary()].
#'
#' @section Default VCE and small-sample corrections:
#'
#' `vcov` defaults to `"iid"` (classical) and `small` defaults to `FALSE`,
#' matching Stata's `ivreg2` defaults. With `small = FALSE`, inference uses
#' z and chi-squared statistics and the large-sample variance normalization.
#' Set `small = TRUE` for t and F statistics and the finite-sample `N - K`
#' correction, matching Stata's `small` option. `small` is the single, uniform
#' control for finite-sample corrections across every VCE type (`"iid"`,
#' `"robust"`, `"HAC"`, `"AC"`); it is not a VCE selector. There are no
#' `"HC0"`/`"HC1"` values --- use `vcov = "robust"` and toggle `small`.
#'
#' @section Sigma (RMSE) normalization:
#'
#' Under `small = FALSE`, \eqn{\sigma^2 = RSS / (N - dofminus)}. Under
#' `small = TRUE`, \eqn{\sigma^2 = RSS / (N - K - dofminus - sdofminus)}. This
#' differs from [lm()], which always divides the residual sum of squares by its
#' residual degrees of freedom. Because ivreg2r reports the large-sample sigma
#' by default, its RMSE and any quantity built from sigma differ from `lm()`
#' unless `small = TRUE`.
#'
#' @section Standard errors relative to lm():
#'
#' Coefficients match [lm()] for OLS. Standard errors and the full
#' variance-covariance matrix match `lm()` only when `small = TRUE`; the
#' default `small = FALSE` reports Stata-convention standard errors without the
#' `N - K` finite-sample correction. Under analytic (`"aweight"`) or
#' probability (`"pweight"`) weights, sigma additionally differs from
#' `lm(..., weights = w)` by a factor of `sqrt(N / sum(w))`, because `lm()`
#' uses raw (unnormalized) weights while ivreg2r normalizes them (see below).
#'
#' @section Weight normalization:
#'
#' Analytic and probability weights are normalized internally to sum to `N`,
#' following Stata. This makes sigma scale-invariant for those types:
#' multiplying every weight by a constant leaves all coefficients, standard
#' errors, and test statistics unchanged. Frequency (`"fweight"`) and
#' importance (`"iweight"`) weights are *not* normalized; they redefine `N` as
#' `sum(weights)`, so their scale is meaningful by construction. Probability
#' weights force a robust VCE. See the `weight_type` argument of [ivreg2()].
#'
#' @section Degrees-of-freedom adjustments:
#'
#' `dofminus` is a large-sample adjustment subtracted from `N` in large-sample
#' variance formulas (for example when fixed effects have been partialled out
#' upstream). `sdofminus` is a small-sample adjustment subtracted from the
#' residual degrees of freedom alongside `K` (for example the count of
#' regressors removed by `partial`). Both mirror Stata's `dofminus()` and
#' `sdofminus()` options. Stata's own source labels `dofminus` a large-sample
#' adjustment (e.g. number of fixed effects) and `sdofminus` a small-sample one
#' (e.g. number of partialled-out regressors).
#'
#' @section Confidence intervals and test statistics:
#'
#' [confint()] and the model Wald test use the standard normal and chi-squared
#' distributions when `small = FALSE`, and the t and F distributions when
#' `small = TRUE`. The reported model test is therefore a chi-squared Wald
#' statistic by default and an F statistic under `small = TRUE`.
#'
#' @name ivreg2r-conventions
#' @seealso [ivreg2()] for the estimation function and the full argument list;
#'   [ivreg2r-glossary] for the diagnostic-output acronyms.
#' @family ivreg2r reference
NULL

#' Glossary of diagnostic-output acronyms
#'
#' [summary()] and [print()] report a battery of identification, weak-instrument,
#' and overidentification diagnostics, many named by acronym. This topic glosses
#' every acronym that reaches the console. The underlying statistics are stored
#' on the fitted object under `fit$diagnostics`; the headline ones are also in
#' [glance()].
#'
#' @section Identification and weak instruments:
#'
#' \describe{
#'   \item{Underidentification test}{Tests whether the excluded instruments
#'     have rank high enough to identify the endogenous regressors. Reported as
#'     the Anderson (1951) canonical-correlations LM statistic under
#'     `vcov = "iid"`, and the Kleibergen-Paap rk LM statistic under robust,
#'     cluster, or HAC errors.}
#'   \item{Weak identification test}{Gauges instrument strength against the
#'     Stock-Yogo (2005) critical values. Reported as the Cragg-Donald (1993)
#'     Wald F under `vcov = "iid"`, and the Kleibergen-Paap rk Wald F under
#'     robust, cluster, or HAC errors.}
#'   \item{Cragg-Donald (CD)}{The iid weak-identification Wald F (and the
#'     underlying eigenvalue statistic); see "Weak identification test".}
#'   \item{Kleibergen-Paap (KP) rk}{Rank-based LM and Wald statistics that
#'     replace the Anderson and Cragg-Donald statistics when errors are not iid.}
#'   \item{Anderson-Rubin (AR)}{A weak-instrument-robust test of the joint
#'     significance of the endogenous regressors: valid regardless of
#'     instrument strength.}
#'   \item{Stock-Wright (S)}{A weak-instrument-robust LM statistic for the same
#'     hypothesis; equal to the Hansen J statistic under iid errors.}
#' }
#'
#' @section First-stage diagnostics:
#'
#' \describe{
#'   \item{Shea partial R-squared}{A partial R-squared for the first stage that
#'     accounts for correlation among instruments; the relevant measure with
#'     more than one endogenous regressor.}
#'   \item{SW F (Sanderson-Windmeijer)}{A conditional first-stage F statistic
#'     for the identification of an individual endogenous regressor, given the
#'     others.}
#'   \item{AP (Angrist-Pischke)}{A conditional first-stage F (or chi-squared,
#'     depending on the VCE) for an individual endogenous regressor.}
#' }
#'
#' @section Overidentification and endogeneity:
#'
#' \describe{
#'   \item{Sargan / Hansen J}{The overidentification test: Sargan (1958) under
#'     `vcov = "iid"`, and the Hansen J statistic under robust, cluster, or HAC
#'     errors. Tests the joint validity of the overidentifying restrictions.}
#'   \item{C-statistic (C-stat)}{A difference-of-J statistic used for the
#'     endogeneity test (`endog`), the orthogonality test (`orthog`), and the
#'     redundancy test (`redundant`).}
#' }
#'
#' @section VCE and estimator labels:
#'
#' \describe{
#'   \item{Kiefer}{The Kiefer (1980) VCE: autocorrelation-consistent with a
#'     Truncated kernel at bandwidth equal to the full time span.}
#'   \item{Driscoll-Kraay (DK)}{The Driscoll-Kraay (1998) panel VCE, robust to
#'     cross-sectional dependence.}
#'   \item{dofminus / sdofminus}{Large- and small-sample degrees-of-freedom
#'     adjustments; see [ivreg2r-conventions].}
#'   \item{psd0 / psda}{Positive-semidefinite corrections for the moment
#'     covariance matrix, zeroing (`psd0`) or taking the absolute value of
#'     (`psda`) any negative eigenvalues.}
#'   \item{COVIV}{The "covariance at the IV estimates" variance matrix for
#'     LIML and k-class estimators, robust to misspecification of the LIML
#'     model.}
#'   \item{HOLS}{Cragg's (1983) heteroskedastic OLS estimator, obtained from a
#'     no-endogenous-regressor model under `method = "gmm2s"` with a robust VCE.}
#'   \item{k-class}{The family of IV estimators indexed by a scalar k
#'     (2SLS at k = 1, OLS at k = 0, LIML at k equal to the LIML eigenvalue).}
#' }
#'
#' @name ivreg2r-glossary
#' @seealso [summary.ivreg2()] for the console output; [ivreg2r-conventions]
#'   for the statistical conventions; [ivreg2()] for the arguments that request
#'   each test.
#' @family ivreg2r reference
NULL
