# ============================================================================
# Broom methods: tidy(), glance(), augment() for ivreg2 objects
# ============================================================================

#' @importFrom generics tidy glance augment
#' @importFrom stats confint
NULL

#' @export
generics::tidy

#' @export
generics::glance

#' @export
generics::augment


# --------------------------------------------------------------------------
# Internal helper
# --------------------------------------------------------------------------

#' Safely extract a numeric diagnostic field
#' @keywords internal
#' @noRd
.safe_diag <- function(diag, test, field) {
  val <- diag[[test]][[field]]
  if (is.null(val) || length(val) == 0L) return(NA_real_)
  if (is.na(val)) return(NA_real_)
  as.double(val)
}


# --------------------------------------------------------------------------
# tidy.ivreg2
# --------------------------------------------------------------------------
#' Tidy an ivreg2 object
#'
#' Constructs a tibble summarizing coefficient estimates, standard errors,
#' test statistics, and p-values.
#'
#' @param x An object of class `"ivreg2"`.
#' @param conf.int Logical: include confidence intervals? Default `TRUE`.
#' @param conf.level Confidence level for intervals. Default `0.95`.
#' @param exponentiate Logical: exponentiate the coefficient estimates and
#'   confidence interval bounds? Default `FALSE`. Useful for log-linear models
#'   where `exp(estimate)` gives a multiplicative effect. Standard errors remain
#'   on the original (log) scale, following broom convention.
#' @param ... Additional arguments (ignored).
#' @return A [tibble::tibble()] with columns `term`, `estimate`, `std.error`,
#'   `statistic`, `p.value`, and optionally `conf.low`, `conf.high`.
#' @examples
#' data(mroz)
#' mroz_work <- subset(mroz, inlf == 1)
#' fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
#'               data = mroz_work, vcov = "robust")
#' tidy(fit)
#' tidy(fit, conf.int = FALSE)
#' tidy(fit, exponentiate = TRUE)
#'
#' \donttest{
#' # Compare 2SLS and LIML side-by-side on the help-file baseline spec
#' # (weak first stage; see the LIML example in ?ivreg2 for the framing)
#' fit_liml <- ivreg2(lwage ~ exper + expersq | educ |
#'                      age + kidslt6 + kidsge6,
#'                    data = mroz_work, method = "liml")
#' comparison <- rbind(
#'   cbind(method = "2SLS", tidy(fit)),
#'   cbind(method = "LIML", tidy(fit_liml))
#' )
#' comparison[comparison$term == "educ", c("method", "estimate", "std.error")]
#' }
#' @export
tidy.ivreg2 <- function(x, conf.int = TRUE, conf.level = 0.95,
                         exponentiate = FALSE, ...) {
  cf <- coef(x)
  se <- sqrt(diag(vcov(x)))
  stat <- cf / se
  if (x$small) {
    pval <- 2 * stats::pt(-abs(stat), df = x$df.residual)
  } else {
    pval <- 2 * stats::pnorm(-abs(stat))
  }
  out <- tibble::tibble(
    term      = names(cf),
    estimate  = unname(cf),
    std.error = unname(se),
    statistic = unname(stat),
    p.value   = unname(pval)
  )
  if (conf.int) {
    ci <- confint(x, level = conf.level)
    out$conf.low  <- unname(ci[, 1L])
    out$conf.high <- unname(ci[, 2L])
  }
  if (exponentiate) {
    out$estimate <- exp(out$estimate)
    if (conf.int) {
      out$conf.low  <- exp(out$conf.low)
      out$conf.high <- exp(out$conf.high)
    }
  }
  out
}


# --------------------------------------------------------------------------
# glance.ivreg2
# --------------------------------------------------------------------------
#' Glance at an ivreg2 object
#'
#' Returns a single-row tibble of model-level summary statistics and the
#' headline IV diagnostics. The column set is deliberately compact so that
#' table tools such as \pkg{modelsummary} render a sensible default.
#'
#' @param x An object of class `"ivreg2"`.
#' @param diagnostics Logical: include the headline IV diagnostic columns?
#'   Default `TRUE`. Set to `FALSE` for a goodness-of-fit summary without the
#'   test statistics. Follows the same convention as broom's
#'   \code{glance.ivreg()}.
#' @param ... Additional arguments (ignored).
#' @return A single-row [tibble::tibble()].
#'
#'   **Always present** (9 columns): `r.squared`, `adj.r.squared`, `sigma`,
#'   `statistic` (model F or Wald chi-squared), `p.value`, `df` (model
#'   numerator degrees of freedom), `df.residual`, `nobs`, `vcov_type`.
#'
#'   **When `diagnostics = TRUE`** (default, 6 additional columns): the
#'   headline IV specification tests — `weak_id_stat` (Cragg-Donald Wald F),
#'   `weak_id_robust_stat` (Kleibergen-Paap rk Wald F), `underid_stat` and
#'   `underid_p` (underidentification), `overid_stat` and `overid_p`
#'   (Sargan/Hansen J overidentification).
#'
#'   The remaining stored quantities and configuration flags are **not** in
#'   `glance()` — this keeps the goodness-of-fit block usable in rendered
#'   tables. They remain available as named elements on the fitted object: the
#'   estimation `method`, `lambda`/`kclass_value`/`fuller_parameter`, `coviv`,
#'   `center`, `psd`, `kernel`/`bw`, `kiefer`, `dkraay`, `sw`, cluster counts,
#'   `cue_convergence`, `partial_ct`, `small`, the cross-products `yy`/`yyc`,
#'   ranks and condition numbers (`rank`, `rankzz`, `condxx`, `condzz`), the
#'   log-likelihood `ll`, and the full diagnostic list `x$diagnostics`
#'   (endogeneity, orthogonality, redundancy, Anderson-Rubin, Stock-Wright,
#'   and the Cragg-Donald/Kleibergen-Paap eigenvalues).
#'
#' @details
#' \code{glance()} returns a fixed set of columns for a given value of
#' \code{diagnostics}, using \code{NA} for metrics that do not apply to the
#' fitted model. All diagnostic columns are \code{NA} for OLS models
#' (single-part formula). \code{overid_stat} and \code{overid_p} are also
#' \code{NA} when the model is exactly identified (the number of excluded
#' instruments equals the number of endogenous regressors), and
#' \code{weak_id_robust_stat} is \code{NA} under \code{vcov = "iid"} (the
#' Cragg-Donald F in \code{weak_id_stat} is reported instead of the
#' Kleibergen-Paap F).
#'
#' Set \code{diagnostics = FALSE} for a compact goodness-of-fit summary
#' without the IV test columns.
#'
#' @examples
#' data(mroz)
#' mroz_work <- subset(mroz, inlf == 1)
#' fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
#'               data = mroz_work, vcov = "robust")
#'
#' # Full output with diagnostics
#' glance(fit)
#'
#' # Compact output without diagnostics
#' glance(fit, diagnostics = FALSE)
#'
#' # Extract specific diagnostics
#' glance(fit)[, c("overid_stat", "overid_p")]
#' glance(fit)[, c("weak_id_stat", "weak_id_robust_stat")]
#'
#' # Diagnostics dropped from glance() remain on the fitted object
#' fit$diagnostics$endogeneity
#'
#' \donttest{
#' # Compare Sargan (IID) vs Hansen J (robust)
#' fit_iid <- ivreg2(lwage ~ exper + expersq | educ |
#'                     age + kidslt6 + kidsge6, data = mroz_work)
#' data.frame(
#'   vcov = c("iid", "robust"),
#'   overid = c(glance(fit_iid)$overid_stat, glance(fit)$overid_stat),
#'   overid_p = c(glance(fit_iid)$overid_p, glance(fit)$overid_p)
#' )
#'
#' # A compact modelsummary table built from the curated glance() columns
#' if (requireNamespace("modelsummary", quietly = TRUE)) {
#'   modelsummary::modelsummary(
#'     list("2SLS" = fit),
#'     statistic = "std.error",
#'     gof_map = c("nobs", "r.squared", "weak_id_stat", "overid_stat")
#'   )
#' }
#' }
#' @export
glance.ivreg2 <- function(x, diagnostics = TRUE, ...) {
  diag <- x$diagnostics

  out <- tibble::tibble(
    r.squared     = x$r.squared,
    adj.r.squared = x$adj.r.squared,
    sigma         = x$sigma,
    statistic     = x$model_f %||% NA_real_,
    p.value       = x$model_f_p %||% NA_real_,
    df            = x$model_f_df1 %||% NA_integer_,
    df.residual   = x$df.residual,
    nobs          = as.integer(x$nobs),
    vcov_type     = x$vcov_type
  )

  if (diagnostics) {
    # Overidentification: NA when exactly identified (df == 0)
    overid_stat <- NA_real_
    overid_p    <- NA_real_
    if (!is.null(diag$overid) && diag$overid$df > 0L) {
      overid_stat <- .safe_diag(diag, "overid", "stat")
      overid_p    <- .safe_diag(diag, "overid", "p")
    }
    out$weak_id_stat        <- .safe_diag(diag, "weak_id", "stat")
    out$weak_id_robust_stat <- .safe_diag(diag, "weak_id_robust", "stat")
    out$underid_stat        <- .safe_diag(diag, "underid", "stat")
    out$underid_p           <- .safe_diag(diag, "underid", "p")
    out$overid_stat         <- overid_stat
    out$overid_p            <- overid_p
  }

  out
}


# --------------------------------------------------------------------------
# augment.ivreg2
# --------------------------------------------------------------------------
#' Augment data with ivreg2 model predictions and residuals
#'
#' Adds `.fitted` and `.resid` columns to the model frame (or user-supplied
#' data).
#'
#' @param x An object of class `"ivreg2"`.
#' @param data A data frame to augment. If `NULL` (default), uses the stored
#'   model frame (`x$model`). An error is raised if `model = FALSE` was used
#'   at estimation time and `data` is not supplied.
#' @param ... Additional arguments (ignored).
#' @return A [tibble::tibble()] with all original data columns plus `.fitted`
#'   and `.resid`.
#' @examples
#' data(mroz)
#' mroz_work <- subset(mroz, inlf == 1)
#' fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
#'               data = mroz_work, vcov = "robust")
#' augment(fit) |> head()
#'
#' # Residuals vs fitted
#' aug <- augment(fit)
#' plot(aug$.fitted, aug$.resid, xlab = "Fitted", ylab = "Residuals")
#' abline(h = 0, lty = 2)
#' @export
augment.ivreg2 <- function(x, data = NULL, ...) {
  if (is.null(data)) {
    if (is.null(x$model)) {
      stop("Model frame not available (fitted with `model = FALSE`). ",
           "Supply `data` argument.", call. = FALSE)
    }
    out <- tibble::as_tibble(x$model)
    out$.fitted <- unname(x$fitted.values)
    out$.resid  <- unname(x$residuals)
  } else {
    out <- tibble::as_tibble(data)
    out$.fitted <- unname(stats::napredict(x$na.action, x$fitted.values))
    out$.resid  <- unname(stats::naresid(x$na.action, x$residuals))
  }
  out
}
