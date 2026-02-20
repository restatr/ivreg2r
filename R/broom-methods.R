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
#' @param ... Additional arguments (ignored).
#' @return A [tibble::tibble()] with columns `term`, `estimate`, `std.error`,
#'   `statistic`, `p.value`, and optionally `conf.low`, `conf.high`.
#' @export
tidy.ivreg2 <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
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
#' @param ... Additional arguments (ignored).
#' @return A single-row [tibble::tibble()] with columns:
#'   `r.squared`, `adj.r.squared`, `sigma`, `statistic`, `p.value`, `df`,
#'   `df.residual`, `nobs`, `vcov_type`,
#'   `weak_id_stat`, `weak_id_robust_stat`,
#'   `underid_stat`, `underid_p`,
#'   `overid_stat`, `overid_p`,
#'   `endogeneity_stat`, `endogeneity_p`,
#'   `stock_wright_stat`, `stock_wright_p`, `stock_wright_df`,
#'   `orthog_stat`, `orthog_p`.
#' @export
glance.ivreg2 <- function(x, ...) {
  diag <- x$diagnostics

  # Overidentification: NA when exactly identified (df == 0)
  overid_stat <- NA_real_
  overid_p    <- NA_real_
  if (!is.null(diag$overid) && diag$overid$df > 0L) {
    overid_stat <- .safe_diag(diag, "overid", "stat")
    overid_p    <- .safe_diag(diag, "overid", "p")
  }

  tibble::tibble(
    r.squared          = x$r.squared,
    adj.r.squared      = x$adj.r.squared,
    sigma              = x$sigma,
    statistic          = x$model_f %||% NA_real_,
    p.value            = x$model_f_p %||% NA_real_,
    df                 = x$model_f_df1 %||% NA_integer_,
    df.residual        = x$df.residual,
    nobs               = as.integer(x$nobs),
    vcov_type          = x$vcov_type,
    n_clusters1        = x$n_clusters1 %||% NA_integer_,
    n_clusters2        = x$n_clusters2 %||% NA_integer_,
    weak_id_stat       = .safe_diag(diag, "weak_id", "stat"),
    weak_id_robust_stat = .safe_diag(diag, "weak_id_robust", "stat"),
    underid_stat       = .safe_diag(diag, "underid", "stat"),
    underid_p          = .safe_diag(diag, "underid", "p"),
    overid_stat        = overid_stat,
    overid_p           = overid_p,
    endogeneity_stat   = .safe_diag(diag, "endogeneity", "stat"),
    endogeneity_p      = .safe_diag(diag, "endogeneity", "p"),
    stock_wright_stat  = .safe_diag(diag, "stock_wright", "stat"),
    stock_wright_p     = .safe_diag(diag, "stock_wright", "p"),
    stock_wright_df    = .safe_diag(diag, "stock_wright", "df"),
    orthog_stat        = .safe_diag(diag, "orthog", "stat"),
    orthog_p           = .safe_diag(diag, "orthog", "p"),
    rf_f_stat          = if (!is.null(x$reduced_form) &&
                              x$reduced_form$mode == "rf") {
                           x$reduced_form$f_stat %||% NA_real_
                         } else NA_real_,
    rf_f_p             = if (!is.null(x$reduced_form) &&
                              x$reduced_form$mode == "rf") {
                           x$reduced_form$f_p %||% NA_real_
                         } else NA_real_
  )
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
