# --------------------------------------------------------------------------
# .chol_solve
# --------------------------------------------------------------------------
#' Cholesky-first linear solver
#'
#' Tries Cholesky factorization first (fast, stable for PD matrices),
#' falls back to QR if Cholesky fails. Matches Stata's `cholqrsolve()`.
#'
#' @param A Symmetric positive-definite matrix.
#' @param b Right-hand side vector or matrix.
#' @return Solution x such that A x = b.
#' @keywords internal
.chol_solve <- function(A, b) {
  R <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(R)) {
    backsolve(R, forwardsolve(t(R), b))
  } else {
    qr.solve(A, b)
  }
}


# --------------------------------------------------------------------------
# .compute_omega
# --------------------------------------------------------------------------
#' Compute L x L moment covariance matrix Omega
#'
#' Computes the instrument-space score covariance for the Hansen J test.
#' This is a *new* computation in Z-space, not the K x K meat from vcov-robust.R.
#'
#' @param Z N x L instrument matrix.
#' @param residuals N x 1 residual vector.
#' @param weights Normalized weights (sum to N), or NULL.
#' @param cluster_vec Cluster membership vector, or NULL.
#' @param N Number of observations.
#' @return L x L symmetric matrix Omega.
#' @keywords internal
.compute_omega <- function(Z, residuals, weights, cluster_vec, N) {
  if (!is.null(cluster_vec)) {
    # Cluster path
    if (!is.null(weights)) {
      sqrt_w <- sqrt(weights)
      scores <- (sqrt_w * Z) * (sqrt_w * residuals)
    } else {
      scores <- Z * residuals
    }
    cluster_scores <- rowsum(scores, cluster_vec, reorder = FALSE)
    Omega <- crossprod(cluster_scores) / N
  } else {
    # HC path (no clusters)
    # Scores are s_i = w_i * z_i * e_i, so the outer product has w_i^2.
    # Stata (livreg2.do line 248): wv = (e .* wvar * wf):^2
    if (!is.null(weights)) {
      wv <- weights^2 * residuals^2
    } else {
      wv <- residuals^2
    }
    Omega <- crossprod(Z, wv * Z) / N
  }
  # Force symmetry
  (Omega + t(Omega)) / 2
}


# --------------------------------------------------------------------------
# .sargan_test
# --------------------------------------------------------------------------
#' Sargan overidentification test (IID path)
#'
#' Computes the Sargan statistic as a quadratic form in the instrument-space
#' score, normalized by the large-sample sigma^2 (rss/N, always).
#'
#' @param Z N x L instrument matrix.
#' @param residuals N x 1 residual vector.
#' @param rss Residual sum of squares (already weighted if applicable).
#' @param N Number of observations.
#' @param overid_df Degrees of freedom (L - K).
#' @param weights Normalized weights or NULL.
#' @return Named list with `stat`, `p`, `df`, `test_name`.
#' @keywords internal
.sargan_test <- function(Z, residuals, rss, N, overid_df, weights) {
  sigma_sq <- rss / N  # always large-sample, regardless of small

  if (!is.null(weights)) {
    Ze <- crossprod(Z, weights * residuals)
    ZtZ <- crossprod(Z, weights * Z)
  } else {
    Ze <- crossprod(Z, residuals)
    ZtZ <- crossprod(Z)
  }

  ZtZ_inv_Ze <- .chol_solve(ZtZ, Ze)
  stat <- drop(crossprod(Ze, ZtZ_inv_Ze)) / sigma_sq
  p <- stats::pchisq(stat, df = overid_df, lower.tail = FALSE)

  list(stat = stat, p = p, df = overid_df, test_name = "Sargan")
}


# --------------------------------------------------------------------------
# .compute_j_with_omega
# --------------------------------------------------------------------------
#' Compute J statistic given a pre-computed Omega matrix
#'
#' Performs closed-form 2-step efficient GMM re-estimation using a provided
#' Omega (instrument-space score covariance) and returns the J statistic.
#' This is the inner loop shared by `.hansen_j_test()` and
#' `.compute_endogeneity_test()`.
#'
#' @param Z N x L instrument matrix.
#' @param X N x K regressor matrix.
#' @param y N x 1 response vector.
#' @param Omega L x L symmetric score covariance matrix.
#' @param weights Normalized weights (sum to N), or NULL.
#' @param N Number of observations.
#' @return Scalar J statistic, or `NA_real_` if Omega is rank-deficient.
#' @keywords internal
.compute_j_with_omega <- function(Z, X, y, Omega, weights, N) {
  L <- ncol(Z)

  # Check Omega rank via Cholesky
  R_chol <- tryCatch(chol(Omega), error = function(e) NULL)
  if (is.null(R_chol)) {
    if (qr(Omega)$rank < L) {
      return(NA_real_)
    }
    R_chol <- NULL
  }

  # Solve Omega^{-1} %*% b
  omega_inv <- function(b) {
    if (!is.null(R_chol)) {
      backsolve(R_chol, forwardsolve(t(R_chol), b))
    } else {
      qr.solve(Omega, b)
    }
  }

  # 2-step GMM re-estimation (closed-form)
  if (!is.null(weights)) {
    QXZ <- crossprod(X, weights * Z) / N   # K x L
    QZy <- crossprod(Z, weights * y) / N   # L x 1
  } else {
    QXZ <- crossprod(X, Z) / N             # K x L
    QZy <- crossprod(Z, y) / N             # L x 1
  }

  aux1 <- omega_inv(t(QXZ))   # L x K: Omega^{-1} QXZ'
  aux2 <- omega_inv(QZy)      # L x 1: Omega^{-1} QZy

  M <- QXZ %*% aux1           # K x K
  M <- (M + t(M)) / 2         # force symmetry

  beta_2s <- .chol_solve(M, QXZ %*% aux2)  # K x 1

  # J statistic from 2-step residuals
  e_2s <- y - X %*% beta_2s
  if (!is.null(weights)) {
    gbar <- crossprod(Z, weights * e_2s) / N  # L x 1
  } else {
    gbar <- crossprod(Z, e_2s) / N            # L x 1
  }

  Omega_inv_gbar <- omega_inv(gbar)
  drop(N * crossprod(gbar, Omega_inv_gbar))
}


# --------------------------------------------------------------------------
# .hansen_j_test
# --------------------------------------------------------------------------
#' Hansen J overidentification test (robust/cluster path)
#'
#' Computes Omega from first-step residuals and delegates to
#' `.compute_j_with_omega()` for 2-step GMM re-estimation and J statistic.
#'
#' @param Z N x L instrument matrix.
#' @param X N x K regressor matrix.
#' @param y N x 1 response vector.
#' @param residuals N x 1 first-step residual vector.
#' @param weights Normalized weights or NULL.
#' @param cluster_vec Cluster membership vector, or NULL.
#' @param N Number of observations.
#' @param K Number of regressors.
#' @param L Number of instruments.
#' @param overid_df Degrees of freedom (L - K).
#' @return Named list with `stat`, `p`, `df`, `test_name`.
#' @keywords internal
.hansen_j_test <- function(Z, X, y, residuals, weights, cluster_vec,
                           N, K, L, overid_df) {
  Omega <- .compute_omega(Z, residuals, weights, cluster_vec, N)
  J <- .compute_j_with_omega(Z, X, y, Omega, weights, N)

  if (is.na(J)) {
    warning("Omega is rank-deficient; Hansen J statistic not computed.",
            call. = FALSE)
    return(list(stat = NA_real_, p = NA_real_, df = overid_df,
                test_name = "Hansen J"))
  }

  p <- stats::pchisq(J, df = overid_df, lower.tail = FALSE)
  list(stat = J, p = p, df = overid_df, test_name = "Hansen J")
}


# --------------------------------------------------------------------------
# .compute_overid_test
# --------------------------------------------------------------------------
#' Dispatcher for overidentification tests
#'
#' Returns Sargan (IID) or Hansen J (robust/cluster), or a zero-stat
#' placeholder for exactly identified models. Returns NULL for OLS.
#'
#' @param Z N x L instrument matrix.
#' @param X N x K regressor matrix.
#' @param y N x 1 response vector.
#' @param residuals N x 1 residual vector.
#' @param rss Residual sum of squares.
#' @param weights Normalized weights or NULL.
#' @param cluster_vec Cluster membership vector, or NULL.
#' @param vcov_type Character: "iid", "HC0", "HC1", or "CL".
#' @param is_iv Logical: TRUE if this is an IV model.
#' @param N Number of observations.
#' @param K Number of regressors.
#' @param L Number of instruments.
#' @param overid_df Degrees of freedom (L - K).
#' @return Named list with `stat`, `p`, `df`, `test_name`, or NULL.
#' @note The Sargan statistic is normalized by the large-sample sigma-squared
#'   (e'e/N). No small-sample correction is applied even when
#'   \code{small = TRUE}. This matches Stata's \code{ivreg2}.
#' @keywords internal
.compute_overid_test <- function(Z, X, y, residuals, rss, weights,
                                 cluster_vec, vcov_type, is_iv,
                                 N, K, L, overid_df) {
  if (!is_iv) return(NULL)

  if (overid_df == 0L) {
    test_name <- if (vcov_type == "iid") "Sargan" else "Hansen J"
    return(list(stat = 0, p = NA_real_, df = 0L, test_name = test_name))
  }

  if (vcov_type == "iid") {
    .sargan_test(Z, residuals, rss, N, overid_df, weights)
  } else {
    .hansen_j_test(Z, X, y, residuals, weights, cluster_vec,
                   N, K, L, overid_df)
  }
}
