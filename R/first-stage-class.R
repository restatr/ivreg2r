# ============================================================================
# First-stage model objects (Ticket S1)
# ============================================================================
# Lightweight S3 class wrapping each first-stage regression so users can
# extract coefficients, VCV, residuals, etc. for publication tables.
# Mirrors Stata's `savefirst` / `_ivreg2_<varname>` stored estimations.


# --------------------------------------------------------------------------
# first_stage() generic + method
# --------------------------------------------------------------------------

#' Extract first-stage regression objects
#'
#' Returns a named list of `ivreg2_first_stage` objects, one per endogenous
#' variable. Each `ivreg2_first_stage` object supports the following
#' methods:
#'
#' \describe{
#'   \item{[coef()]}{First-stage coefficient vector.}
#'   \item{[vcov()]}{First-stage variance-covariance matrix.}
#'   \item{[residuals()]}{First-stage residuals.}
#'   \item{[fitted()]}{First-stage fitted values.}
#'   \item{[nobs()]}{Number of observations used in the first stage.}
#'   \item{[confint()]}{Confidence intervals for the first-stage
#'     coefficients.}
#'   \item{[print()]}{Compact one-line-per-coefficient display.}
#'   \item{[summary()]}{Full first-stage regression table: coefficients,
#'     standard errors, t-values, p-values, the first-stage F-statistic, and
#'     the partial R-squared (its `print()` method renders the table).}
#'   \item{[tidy()]}{Broom-style tibble with one row per coefficient
#'     (`term`, `estimate`, `std.error`, `statistic`, `p.value`, and,
#'     optionally, confidence limits).}
#'   \item{[glance()]}{Broom-style one-row tibble of first-stage fit
#'     statistics (sigma, degrees of freedom, F-statistic, partial
#'     R-squared, and related diagnostics).}
#' }
#'
#' @param x A fitted model object.
#' @param ... Additional arguments (ignored).
#' @return A named list of `ivreg2_first_stage` objects.
#'
#' @section Small-sample correction: the first-stage F-statistic (surfaced
#'   via `glance()` and `summary()` on each `ivreg2_first_stage` object)
#'   always uses small-sample degrees of freedom, regardless of the main
#'   model's `small` argument. This matches Stata's `ivreg2`. See
#'   [ivreg2r-conventions] for the full statement of which corrections
#'   `small` controls.
#'
#' @examples
#' # Mirrors the Stata `ivreg2` help-file example at line 1325, which
#' # demonstrates the equivalence of the Kleibergen-Paap rk Wald F and the
#' # first-stage F with a single endogenous regressor (and, in passing,
#' # the `savefirst` option that `first_stage = TRUE` corresponds to).
#' data(mroz)
#' mroz_work <- subset(mroz, inlf == 1)
#' fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
#'               data = mroz_work, vcov = "robust", first_stage = TRUE)
#' fs <- first_stage(fit)
#' coef(fs$educ)
#' summary(fs$educ)
#' # KP rk Wald F = the robust first-stage F (single endogenous regressor)
#' c(kp_rk_wald_F = fit$diagnostics$weak_id_robust$stat,
#'   first_stage_F = glance(fs$educ)$f_stat)
#'
#' @family ivreg2 methods
#' @seealso [ivreg2()]
#' @export
first_stage <- function(x, ...) UseMethod("first_stage")

#' @rdname first_stage
#' @export
first_stage.ivreg2 <- function(x, ...) {
  if (is.null(x$first_stage_models)) {
    stop("First-stage models not stored. Re-fit with `first_stage = TRUE`.",
         call. = FALSE)
  }
  x$first_stage_models
}


# --------------------------------------------------------------------------
# .new_ivreg2_first_stage
# --------------------------------------------------------------------------
#' Construct an ivreg2_first_stage object
#' @keywords internal
#' @noRd
.new_ivreg2_first_stage <- function(coefficients, vcov, residuals,
                                     fitted.values, sigma, df.residual,
                                     nobs, rank, depvar, vcov_type,
                                     cluster_var = NULL, n_clusters = NULL,
                                     kernel = NULL, bw = NULL,
                                     f_stat, f_p, f_df1, f_df2,
                                     partial_r2, shea_partial_r2,
                                     sw_f = NULL, sw_f_p = NULL,
                                     sw_f_df1 = NULL, sw_f_df2 = NULL,
                                     ap_f = NULL, ap_f_p = NULL,
                                     ap_f_df1 = NULL, ap_f_df2 = NULL) {
  structure(
    list(
      coefficients    = coefficients,
      vcov            = vcov,
      residuals       = residuals,
      fitted.values   = fitted.values,
      sigma           = sigma,
      df.residual     = df.residual,
      nobs            = nobs,
      rank            = rank,
      depvar          = depvar,
      vcov_type       = vcov_type,
      cluster_var     = cluster_var,
      n_clusters      = n_clusters,
      kernel          = kernel,
      bw              = bw,
      f_stat          = f_stat,
      f_p             = f_p,
      f_df1           = f_df1,
      f_df2           = f_df2,
      partial_r2      = partial_r2,
      shea_partial_r2 = shea_partial_r2,
      sw_f            = sw_f,
      sw_f_p          = sw_f_p,
      sw_f_df1        = sw_f_df1,
      sw_f_df2        = sw_f_df2,
      ap_f            = ap_f,
      ap_f_p          = ap_f_p,
      ap_f_df1        = ap_f_df1,
      ap_f_df2        = ap_f_df2
    ),
    class = "ivreg2_first_stage"
  )
}


# --------------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------------

#' @export
coef.ivreg2_first_stage <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.ivreg2_first_stage <- function(object, ...) {
  object$vcov
}

#' @export
residuals.ivreg2_first_stage <- function(object, ...) {
  object$residuals
}

#' @export
fitted.ivreg2_first_stage <- function(object, ...) {
  object$fitted.values
}

#' @export
nobs.ivreg2_first_stage <- function(object, ...) {
  object$nobs
}

#' @export
print.ivreg2_first_stage <- function(x,
                                      digits = max(3L, getOption("digits") - 3L),
                                      ...) {
  cat("\nFirst-stage regression: ", x$depvar, " ~ instruments\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @export
summary.ivreg2_first_stage <- function(object, ...) {
  cf <- coef(object)
  se <- sqrt(diag(vcov(object)))
  stat <- cf / se
  # First-stage always uses small-sample (t-test)
  pval <- 2 * stats::pt(-abs(stat), df = object$df.residual)
  coef_table <- cbind(Estimate = cf, `Std. Error` = se,
                      `t value` = stat, `Pr(>|t|)` = pval)
  structure(
    c(object, list(coef_table = coef_table)),
    class = "summary.ivreg2_first_stage"
  )
}

#' @export
print.summary.ivreg2_first_stage <- function(x,
                                              digits = max(3L, getOption("digits") - 3L),
                                              signif.stars = getOption("show.signif.stars", TRUE),
                                              ...) {
  cat("\nFirst-stage regression: ", x$depvar, " ~ instruments\n", sep = "")
  cat("Observations: ", format(x$nobs, big.mark = ","), "\n", sep = "")
  cat("VCV type:     ", x$vcov_type, "\n", sep = "")

  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coef_table, digits = digits,
                      signif.stars = signif.stars,
                      na.print = "NA", ...)

  cat("---\n")
  cat("Root MSE:      ", formatC(x$sigma, digits = 4, format = "f"), "\n")
  if (!is.na(x$f_stat)) {
    cat("F(", x$f_df1, ", ", x$f_df2, ") = ",
        formatC(x$f_stat, digits = 2, format = "f"),
        " (p = ", formatC(x$f_p, digits = 4, format = "f"), ")\n", sep = "")
  }
  if (!is.na(x$partial_r2)) {
    cat("Partial R2:    ", formatC(x$partial_r2, digits = 4, format = "f"),
        "\n")
  }
  if (!is.null(x$shea_partial_r2) && !is.na(x$shea_partial_r2)) {
    cat("Shea PR2:      ", formatC(x$shea_partial_r2, digits = 4, format = "f"),
        "\n")
  }

  cat("\n")
  invisible(x)
}

#' @export
confint.ivreg2_first_stage <- function(object, parm, level = 0.95, ...) {
  cf <- coef(object)
  ses <- sqrt(diag(vcov(object)))
  a <- (1 - level) / 2
  crit <- stats::qt(a, df = object$df.residual)
  ci <- cbind(cf + crit * ses, cf - crit * ses)
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE,
                       digits = 3), "%")
  colnames(ci) <- pct
  if (!missing(parm)) {
    if (is.character(parm)) parm <- match(parm, names(cf))
    ci <- ci[parm, , drop = FALSE]
  }
  ci
}


# --------------------------------------------------------------------------
# Broom methods for ivreg2_first_stage
# --------------------------------------------------------------------------

#' @export
tidy.ivreg2_first_stage <- function(x, conf.int = TRUE, conf.level = 0.95,
                                     exponentiate = FALSE, ...) {
  cf <- coef(x)
  se <- sqrt(diag(vcov(x)))
  stat <- cf / se
  pval <- 2 * stats::pt(-abs(stat), df = x$df.residual)
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

#' @export
glance.ivreg2_first_stage <- function(x, ...) {
  tibble::tibble(
    sigma           = x$sigma,
    df.residual     = x$df.residual,
    nobs            = as.integer(x$nobs),
    rank            = as.integer(x$rank),
    vcov_type       = x$vcov_type,
    f_stat          = x$f_stat,
    f_p             = x$f_p,
    f_df1           = x$f_df1,
    f_df2           = x$f_df2,
    partial_r2      = x$partial_r2,
    shea_partial_r2 = x$shea_partial_r2 %||% NA_real_,
    sw_f            = x$sw_f %||% NA_real_,
    ap_f            = x$ap_f %||% NA_real_
  )
}


# --------------------------------------------------------------------------
# Builder: assemble first-stage model objects from existing diagnostics
# --------------------------------------------------------------------------
#' Build first-stage model objects from .compute_first_stage() output
#'
#' @param fs_results Output from `.compute_first_stage()` — named list of
#'   per-endogenous-variable diagnostic results.
#' @param Z N x L instrument matrix.
#' @param weights Normalized weights or NULL.
#' @param cluster_vec Cluster membership vector or NULL.
#' @param vcov_type Effective VCE type.
#' @param N,L,M Integer dimensions.
#' @param dofminus,sdofminus DoF adjustments.
#' @param weight_type Weight type string.
#' @param kernel,bw,time_index HAC parameters.
#' @param center Logical: centering flag.
#' @param cluster_var Character: cluster variable name(s).
#' @return Named list of `ivreg2_first_stage` objects.
#' @keywords internal
#' @noRd
.build_first_stage_models <- function(fs_results, Z, weights, cluster_vec,
                                       vcov_type, N, L, M,
                                       dofminus, sdofminus,
                                       weight_type, kernel, bw,
                                       time_index, center,
                                       cluster_var = NULL) {
  sqrt_w <- if (!is.null(weights)) sqrt(weights) else NULL
  endo_names <- names(fs_results)
  K1 <- length(endo_names)

  # Bread: (Z'WZ)^{-1}
  if (is.null(sqrt_w)) {
    ZtWZ <- crossprod(Z)
  } else {
    ZtWZ <- crossprod(sqrt_w * Z)
  }
  ZtWZ_inv <- chol2inv(chol(ZtWZ))

  # Clustered models use M-1 as residual df (matching main model and Stata's
  # ereturn display behavior). Non-clustered uses N - L - dofminus - sdofminus.
  df_resid <- if (!is.null(cluster_vec)) {
    as.integer(M - 1L)
  } else {
    as.integer(N - L - dofminus - sdofminus)
  }

  # Determine if cluster+kernel (Thompson two-way) — needs special meat path
  has_cluster_kernel <- !is.null(cluster_vec) && !is.null(kernel)

  models <- vector("list", K1)
  names(models) <- endo_names

  for (j in seq_len(K1)) {
    nm <- endo_names[j]
    fs <- fs_results[[nm]]
    resid_j <- fs$residuals

    # --- Compute VCV with PostFirstRF scaling ---
    if (vcov_type == "iid") {
      rss <- if (is.null(weights)) {
        drop(crossprod(resid_j))
      } else {
        sum(weights * resid_j^2)
      }
      sigma2 <- rss / (N - dofminus)
      vcov_raw <- sigma2 * ZtWZ_inv
    } else if (has_cluster_kernel) {
      # Thompson two-way: cluster + kernel — same logic as .compute_first_stage
      if (is.list(cluster_vec)) {
        scores <- .cl_scores(Z, resid_j, weights,
                             center = center, weight_type = weight_type)
        shat1 <- crossprod(rowsum(scores, cluster_vec[[1L]], reorder = FALSE))
        shat1 <- (shat1 + t(shat1)) / 2
        shat2 <- .cluster_kernel_meat(Z, resid_j, time_index, kernel, bw,
                                       weights, weight_type, center = center)
        shat3 <- .hac_scores_meat(scores, time_index, kernel, bw)
        meat <- shat1 + shat2 - shat3
      } else {
        meat <- .cluster_kernel_meat(Z, resid_j, time_index, kernel, bw,
                                     weights, weight_type, center = center)
      }
      vcov_raw <- ZtWZ_inv %*% meat %*% ZtWZ_inv
    } else if (!is.null(cluster_vec)) {
      scores <- .cl_scores(Z, resid_j, weights,
                           center = center, weight_type = weight_type)
      meat <- .cluster_meat(scores, cluster_vec)
      vcov_raw <- ZtWZ_inv %*% meat %*% ZtWZ_inv
    } else if (!is.null(kernel)) {
      meat <- .hac_meat(Z, resid_j, time_index, kernel, bw,
                        weights, weight_type, center = center)
      vcov_raw <- ZtWZ_inv %*% meat %*% ZtWZ_inv
    } else {
      meat <- .hc_meat(Z, resid_j, weights, weight_type, center = center)
      vcov_raw <- ZtWZ_inv %*% meat %*% ZtWZ_inv
    }

    # PostFirstRF small-sample correction (always applied)
    if (!is.null(cluster_vec)) {
      vcov_raw <- vcov_raw * (N - 1) / (N - L - sdofminus) *
        M / (M - 1)
    } else {
      vcov_raw <- vcov_raw * (N - dofminus) / (N - L - dofminus - sdofminus)
    }

    vcov_mat <- (vcov_raw + t(vcov_raw)) / 2
    colnames(vcov_mat) <- rownames(vcov_mat) <- colnames(Z)

    models[[nm]] <- .new_ivreg2_first_stage(
      coefficients    = fs$coefficients,
      vcov            = vcov_mat,
      residuals       = fs$residuals,
      fitted.values   = fs$fitted.values,
      sigma           = fs$rmse,
      df.residual     = df_resid,
      nobs            = N,
      rank            = as.integer(L),
      depvar          = nm,
      vcov_type       = vcov_type,
      cluster_var     = cluster_var,
      n_clusters      = M,
      kernel          = kernel,
      bw              = bw,
      f_stat          = fs$f_stat,
      f_p             = fs$f_p,
      f_df1           = fs$f_df1,
      f_df2           = fs$f_df2,
      partial_r2      = fs$partial_r2,
      shea_partial_r2 = fs$shea_partial_r2,
      sw_f            = fs$sw_f,
      sw_f_p          = fs$sw_f_p,
      sw_f_df1        = fs$sw_f_df1,
      sw_f_df2        = fs$sw_f_df2,
      ap_f            = fs$ap_f,
      ap_f_p          = fs$ap_f_p,
      ap_f_df1        = fs$ap_f_df1,
      ap_f_df2        = fs$ap_f_df2
    )
  }

  models
}
