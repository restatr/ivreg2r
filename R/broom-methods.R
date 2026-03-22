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
#' data(card)
#' fit <- ivreg2(log(wage) ~ exper + expersq + black + south |
#'               educ | nearc4, data = card, vcov = "robust")
#' tidy(fit)
#' tidy(fit, conf.int = FALSE)
#' tidy(fit, exponentiate = TRUE)
#'
#' \donttest{
#' # Compare 2SLS and LIML side-by-side
#' fit_liml <- ivreg2(log(wage) ~ exper + expersq + black + south |
#'                     educ | nearc4 + nearc2, data = card, method = "liml")
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
#' Returns a single-row tibble of model-level summary statistics and
#' diagnostic test results.
#'
#' @param x An object of class `"ivreg2"`.
#' @param diagnostics Logical: include IV diagnostic test columns? Default
#'   `TRUE`. Set to `FALSE` for a compact summary without test statistics.
#'   Follows the same convention as broom's \code{glance.ivreg()}.
#' @param ... Additional arguments (ignored).
#' @return A single-row [tibble::tibble()].
#'
#'   **Always present** (33 columns):
#'   `r.squared`, `adj.r.squared`, `sigma`, `statistic`, `p.value`, `df`,
#'   `df.residual`, `nobs`, `vcov_type`, `small`, `weight_type`, `method`,
#'   `lambda`, `kclass_value`, `fuller_parameter`, `coviv`, `center`, `psd`,
#'   `kernel`, `bw`, `kiefer`, `dkraay`, `n_clusters1`, `n_clusters2`,
#'   `cue_convergence`, `partial_ct`,
#'   `yy`, `yyc`, `rankxx`, `rankzz`, `condxx`, `condzz`, `ll`.
#'
#'   **When `diagnostics = TRUE`** (default, 24 additional columns):
#'   `weak_id_stat`, `weak_id_robust_stat`,
#'   `underid_stat`, `underid_p`,
#'   `overid_stat`, `overid_p`,
#'   `ar_overid_lr_stat`, `ar_overid_lr_p`,
#'   `ar_overid_lin_stat`, `ar_overid_lin_p`, `ar_overid_df`,
#'   `endogeneity_stat`, `endogeneity_p`,
#'   `stock_wright_stat`, `stock_wright_p`, `stock_wright_df`,
#'   `orthog_stat`, `orthog_p`,
#'   `redundancy_stat`, `redundancy_p`,
#'   `rf_f_stat`, `rf_f_p`,
#'   `ccev_min`, `cdev_min`.
#'
#' @details
#' \code{glance()} always returns the same columns for a given value of the
#' \code{diagnostics} argument, using \code{NA} for metrics that do not apply
#' to the fitted model.
#'
#' The \code{small} column indicates whether finite-sample corrections were
#' applied. When \code{small = TRUE}, test statistics are F-distributed;
#' when \code{small = FALSE}, they are chi-squared.
#'
#' When \code{diagnostics = TRUE} (default), the output includes IV
#' specification tests. Columns that are conditionally \code{NA}:
#' \itemize{
#'   \item All diagnostic columns are \code{NA} for OLS models (1-part formula).
#'   \item \code{overid_stat}, \code{overid_p}: also \code{NA} when exactly
#'     identified (number of excluded instruments equals number of endogenous
#'     regressors).
#'   \item \code{weak_id_robust_stat}: \code{NA} when \code{vcov = "iid"}
#'     (Cragg-Donald F is used instead of Kleibergen-Paap).
#'   \item \code{ar_overid_*}: only non-\code{NA} for \code{method = "liml"}
#'     with \code{vcov = "iid"}.
#'   \item \code{orthog_*}, \code{redundancy_*}: \code{NA} unless \code{orthog}
#'     or \code{redundant} was specified.
#'   \item \code{rf_f_*}: \code{NA} unless \code{reduced_form = "rf"}.
#' }
#'
#' Set \code{diagnostics = FALSE} for a compact summary without test statistics.
#'
#' @examples
#' data(card)
#' fit <- ivreg2(log(wage) ~ exper + expersq + black + south |
#'               educ | nearc4 + nearc2, data = card, vcov = "robust")
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
#' \donttest{
#' # Compare Sargan (IID) vs Hansen J (robust)
#' fit_iid <- ivreg2(log(wage) ~ exper + expersq + black + south |
#'                    educ | nearc4 + nearc2, data = card)
#' data.frame(
#'   vcov = c("iid", "robust"),
#'   overid = c(glance(fit_iid)$overid_stat, glance(fit)$overid_stat),
#'   overid_p = c(glance(fit_iid)$overid_p, glance(fit)$overid_p)
#' )
#' }
#' @export
glance.ivreg2 <- function(x, diagnostics = TRUE, ...) {
  diag <- x$diagnostics

  # Overidentification: NA when exactly identified (df == 0)
  overid_stat <- NA_real_
  overid_p    <- NA_real_
  if (!is.null(diag$overid) && diag$overid$df > 0L) {
    overid_stat <- .safe_diag(diag, "overid", "stat")
    overid_p    <- .safe_diag(diag, "overid", "p")
  }

  out <- tibble::tibble(
    r.squared          = x$r.squared,
    adj.r.squared      = x$adj.r.squared,
    sigma              = x$sigma,
    statistic          = x$model_f %||% NA_real_,
    p.value            = x$model_f_p %||% NA_real_,
    df                 = x$model_f_df1 %||% NA_integer_,
    df.residual        = x$df.residual,
    nobs               = as.integer(x$nobs),
    vcov_type          = x$vcov_type,
    small              = x$small,
    weight_type        = x$weight_type %||% "aweight",
    method             = x$method %||% NA_character_,
    lambda             = x$lambda %||% NA_real_,
    kclass_value       = x$kclass_value %||% NA_real_,
    fuller_parameter   = x$fuller_parameter %||% NA_real_,
    coviv              = isTRUE(x$coviv),
    center             = isTRUE(x$center),
    psd                = x$psd %||% NA_character_,
    kernel             = x$kernel %||% NA_character_,
    bw                 = x$bw %||% NA_real_,
    kiefer             = isTRUE(x$kiefer),
    dkraay             = x$dkraay %||% NA_real_,
    n_clusters1        = x$n_clusters1 %||% NA_integer_,
    n_clusters2        = x$n_clusters2 %||% NA_integer_,
    cue_convergence    = x$cue_convergence %||% NA_integer_,
    partial_ct         = x$partial_ct %||% 0L
  )

  # Model-level stored results (always present, not test statistics)
  out$yy                 <- x$yy %||% NA_real_
  out$yyc                <- x$yyc %||% NA_real_
  out$rankxx             <- x$rank %||% NA_integer_
  out$rankzz             <- x$rankzz %||% NA_integer_
  out$condxx             <- x$condxx %||% NA_real_
  out$condzz             <- x$condzz %||% NA_real_
  out$ll                 <- x$ll %||% NA_real_

  if (diagnostics) {
    out$weak_id_stat       <- .safe_diag(diag, "weak_id", "stat")
    out$weak_id_robust_stat <- .safe_diag(diag, "weak_id_robust", "stat")
    out$underid_stat       <- .safe_diag(diag, "underid", "stat")
    out$underid_p          <- .safe_diag(diag, "underid", "p")
    out$overid_stat        <- overid_stat
    out$overid_p           <- overid_p
    out$ar_overid_lr_stat  <- .safe_diag(diag, "anderson_rubin_overid", "lr_stat")
    out$ar_overid_lr_p     <- .safe_diag(diag, "anderson_rubin_overid", "lr_p")
    out$ar_overid_lin_stat <- .safe_diag(diag, "anderson_rubin_overid", "lin_stat")
    out$ar_overid_lin_p    <- .safe_diag(diag, "anderson_rubin_overid", "lin_p")
    out$ar_overid_df       <- .safe_diag(diag, "anderson_rubin_overid", "df")
    out$endogeneity_stat   <- .safe_diag(diag, "endogeneity", "stat")
    out$endogeneity_p      <- .safe_diag(diag, "endogeneity", "p")
    out$stock_wright_stat  <- .safe_diag(diag, "stock_wright", "stat")
    out$stock_wright_p     <- .safe_diag(diag, "stock_wright", "p")
    out$stock_wright_df    <- .safe_diag(diag, "stock_wright", "df")
    out$orthog_stat        <- .safe_diag(diag, "orthog", "stat")
    out$orthog_p           <- .safe_diag(diag, "orthog", "p")
    out$redundancy_stat    <- .safe_diag(diag, "redundancy", "stat")
    out$redundancy_p       <- .safe_diag(diag, "redundancy", "p")
    out$rf_f_stat          <- if (!is.null(x$reduced_form) &&
                                    x$reduced_form$mode == "rf") {
                               x$reduced_form$f_stat %||% NA_real_
                             } else NA_real_
    out$rf_f_p             <- if (!is.null(x$reduced_form) &&
                                    x$reduced_form$mode == "rf") {
                               x$reduced_form$f_p %||% NA_real_
                             } else NA_real_
    out$ccev_min           <- if (!is.null(diag$ccev)) min(diag$ccev) else NA_real_
    out$cdev_min           <- if (!is.null(diag$cdev)) min(diag$cdev) else NA_real_
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
#' data(card)
#' fit <- ivreg2(log(wage) ~ exper + expersq + black + south |
#'               educ | nearc4, data = card)
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
