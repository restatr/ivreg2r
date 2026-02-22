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
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#'   HC path divides by `N - dofminus` (Stata livreg2.do line 326);
#'   cluster path divides by `N` (line 545, no dofminus adjustment).
#' @param kernel Canonical kernel name, or NULL for non-HAC.
#' @param bw Numeric bandwidth, or NULL.
#' @param time_index List from `.build_time_index()`, or NULL.
#' @return L x L symmetric matrix Omega.
#' @keywords internal
.compute_omega <- function(Z, residuals, weights, cluster_vec, N,
                            dofminus = 0L, weight_type = "aweight",
                            kernel = NULL, bw = NULL, time_index = NULL) {
  if (!is.null(cluster_vec) && !is.null(kernel)) {
    # Cluster + kernel (DK or Thompson)
    if (is.list(cluster_vec)) {
      # Thompson: CGM decomposition with kernel-smoothed time dimension
      # Third term is HAC meat (not just HC), matching Stata livreg2.do line 336.
      scores <- .cl_scores(Z, residuals, weights)
      shat1 <- crossprod(rowsum(scores, cluster_vec[[1L]], reorder = FALSE))
      shat1 <- (shat1 + t(shat1)) / 2
      shat2 <- .cluster_kernel_meat(Z, residuals, time_index, kernel, bw,
                                     weights, weight_type)
      shat3 <- .hac_scores_meat(scores, time_index, kernel, bw)
      Omega <- (shat1 + shat2 - shat3) / N
    } else {
      # DK: one-way cluster+kernel on tvar
      meat <- .cluster_kernel_meat(Z, residuals, time_index, kernel, bw,
                                   weights, weight_type)
      Omega <- meat / N
    }
  } else if (!is.null(cluster_vec)) {
    # Cluster path — divisor is N (no dofminus)
    scores <- .cl_scores(Z, residuals, weights)
    Omega <- .cluster_meat(scores, cluster_vec) / N
  } else if (!is.null(kernel)) {
    # HAC path — same structure as HC but with autocovariance lags
    meat <- .hac_meat(Z, residuals, time_index, kernel, bw,
                      weights, weight_type)
    Omega <- meat / (N - dofminus)
  } else {
    # HC path — divisor is N - dofminus (Stata livreg2.do line 326)
    Omega <- .hc_meat(Z, residuals, weights, weight_type) / (N - dofminus)
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
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @return Named list with `stat`, `p`, `df`, `test_name`.
#' @keywords internal
.sargan_test <- function(Z, residuals, rss, N, overid_df, weights,
                          dofminus = 0L) {
  sigma_sq <- rss / (N - dofminus)  # large-sample with dofminus adjustment

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
#' @return Scalar J statistic, or `NA_real_` if Omega is rank-deficient
#'   or the GMM Hessian is singular.
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

  # M (GMM Hessian) can be singular even when Omega is full-rank
  # (e.g., cluster+kernel with few time periods). Return NA gracefully.
  R_M <- tryCatch(chol(M), error = function(e) NULL)
  if (is.null(R_M)) {
    if (qr(M)$rank < ncol(M)) return(NA_real_)
    beta_2s <- qr.solve(M, QXZ %*% aux2)
  } else {
    beta_2s <- backsolve(R_M, forwardsolve(t(R_M), QXZ %*% aux2))
  }

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
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @return Named list with `stat`, `p`, `df`, `test_name`.
#' @keywords internal
.hansen_j_test <- function(Z, X, y, residuals, weights, cluster_vec,
                           N, K, L, overid_df, dofminus = 0L,
                           weight_type = "aweight",
                           kernel = NULL, bw = NULL, time_index = NULL) {
  Omega <- .compute_omega(Z, residuals, weights, cluster_vec, N,
                           dofminus = dofminus, weight_type = weight_type,
                           kernel = kernel, bw = bw, time_index = time_index)
  J <- .compute_j_with_omega(Z, X, y, Omega, weights, N)

  if (is.na(J)) {
    warning("Hansen J statistic not computed; ",
            "singular moment covariance or Hessian matrix.",
            call. = FALSE)
    return(list(stat = NA_real_, p = NA_real_, df = overid_df,
                test_name = "Hansen J"))
  }

  p <- stats::pchisq(J, df = overid_df, lower.tail = FALSE)
  list(stat = J, p = p, df = overid_df, test_name = "Hansen J")
}


# --------------------------------------------------------------------------
# .compute_stock_wright
# --------------------------------------------------------------------------
#' Stock-Wright LM S statistic
#'
#' Weak-instrument-robust LM test of H0: all endogenous coefficients = 0
#' and orthogonality conditions are valid. This is the LM counterpart of
#' the Anderson-Rubin Wald test.
#'
#' Unlike Sargan/Hansen J, the S statistic uses residuals from regressing
#' `y` on the included exogenous regressors only (constraining endogenous
#' coefficients to zero). For IID, a homoskedastic omega is used
#' (`sigma^2 * Z'WZ / N`); for robust/cluster, the HC/CL omega is used.
#'
#' @param Z N x L instrument matrix.
#' @param X N x K regressor matrix.
#' @param y N x 1 response vector.
#' @param weights Normalized weights (sum to N), or NULL.
#' @param cluster_vec Cluster membership vector, or NULL.
#' @param vcov_type Character: `"iid"`, `"HC0"`, `"HC1"`, or `"CL"`.
#' @param N Number of observations.
#' @param K1 Number of endogenous regressors.
#' @param L1 Number of excluded instruments.
#' @param endo_names Character vector of endogenous variable names.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @return Named list with `stat`, `p`, `df`.
#' @keywords internal
.compute_stock_wright <- function(Z, X, y, weights, cluster_vec,
                                   vcov_type, N, K1, L1,
                                   endo_names, dofminus = 0L,
                                   weight_type = "aweight",
                                   kernel = NULL, bw = NULL,
                                   time_index = NULL) {
  # 1. Extract X2 (included exogenous regressors)
  endo_idx <- match(endo_names, colnames(X))
  exog_idx <- setdiff(seq_len(ncol(X)), endo_idx)
  K2 <- length(exog_idx)

  # 2. Constrained residuals: e_0 = y - X2 * b_ols
  #    (endogenous coefficients constrained to zero)
  if (K2 > 0L) {
    X2 <- X[, exog_idx, drop = FALSE]
    if (is.null(weights)) {
      e <- stats::lm.fit(X2, y)$residuals
    } else {
      e <- stats::lm.wfit(X2, y, weights)$residuals
    }
  } else {
    e <- y
  }

  # 3. Compute gbar = Z'We / N
  if (!is.null(weights)) {
    gbar <- crossprod(Z, weights * e) / N
  } else {
    gbar <- crossprod(Z, e) / N
  }

  # 4. Compute omega (VCE-dependent, matching Stata m_omega)
  if (vcov_type %in% c("iid", "AC")) {
    # Homoskedastic / AC: omega = sigma^2_0 * Z'WZ / N
    # (Stata livreg2.do lines 194-235)
    if (!is.null(weights)) {
      rss_0 <- sum(weights * e^2)
      ZWZ <- crossprod(Z, weights * Z)
    } else {
      rss_0 <- sum(e^2)
      ZWZ <- crossprod(Z)
    }
    sigma2_0 <- rss_0 / (N - dofminus)
    if (vcov_type == "AC" && !is.null(kernel)) {
      # AC: use AC meat (Kronecker structure with autocovariances)
      Omega <- .ac_meat(Z, e, time_index, kernel, bw,
                         N, dofminus, weights, weight_type, ZWZ)
    } else {
      Omega <- sigma2_0 * ZWZ / N
    }
  } else {
    # HC/CL/HAC: heteroskedasticity-robust omega
    Omega <- .compute_omega(Z, e, weights, cluster_vec, N,
                             dofminus = dofminus, weight_type = weight_type,
                             kernel = kernel, bw = bw,
                             time_index = time_index)
  }

  # 5. S = N * gbar' * Omega^{-1} * gbar
  R_chol <- tryCatch(chol(Omega), error = function(e) NULL)
  if (is.null(R_chol)) {
    if (qr(Omega)$rank < ncol(Z)) {
      warning("Stock-Wright: Omega is rank-deficient; S statistic not computed.",
              call. = FALSE)
      return(list(stat = NA_real_, p = NA_real_, df = as.integer(L1)))
    }
    Omega_inv_gbar <- qr.solve(Omega, gbar)
  } else {
    Omega_inv_gbar <- backsolve(R_chol, forwardsolve(t(R_chol), gbar))
  }

  stat <- drop(N * crossprod(gbar, Omega_inv_gbar))
  df <- as.integer(L1)
  p <- stats::pchisq(stat, df = df, lower.tail = FALSE)

  list(stat = stat, p = p, df = df)
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
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @return Named list with `stat`, `p`, `df`, `test_name`, or NULL.
#' @note The Sargan statistic is normalized by the large-sample sigma-squared
#'   `e'e/(N-dofminus)`. No small-sample correction is applied even when
#'   \code{small = TRUE}. This matches Stata's \code{ivreg2}.
#' @keywords internal
.compute_overid_test <- function(Z, X, y, residuals, rss, weights,
                                 cluster_vec, vcov_type, is_iv,
                                 N, K, L, overid_df, dofminus = 0L,
                                 weight_type = "aweight",
                                 kernel = NULL, bw = NULL,
                                 time_index = NULL) {
  if (!is_iv) return(NULL)

  if (overid_df == 0L) {
    test_name <- if (vcov_type %in% c("iid", "AC")) "Sargan" else "Hansen J"
    return(list(stat = 0, p = NA_real_, df = 0L, test_name = test_name))
  }

  if (vcov_type == "iid") {
    # Plain IID: standard Sargan
    .sargan_test(Z, residuals, rss, N, overid_df, weights,
                  dofminus = dofminus)
  } else if (vcov_type == "AC") {
    # AC: J-test with AC omega (Stata uses kernel-weighted omega even for
    # the "Sargan" test when a kernel is specified)
    if (!is.null(weights)) {
      ZWZ <- crossprod(Z, weights * Z)
    } else {
      ZWZ <- crossprod(Z)
    }
    Omega <- .ac_meat(Z, residuals, time_index, kernel, bw,
                       N, dofminus, weights, weight_type, ZWZ)
    J <- .compute_j_with_omega(Z, X, y, Omega, weights, N)
    if (is.na(J)) {
      warning("Sargan statistic not computed; ",
              "singular moment covariance or Hessian matrix.",
              call. = FALSE)
      return(list(stat = NA_real_, p = NA_real_, df = overid_df,
                  test_name = "Sargan"))
    }
    p <- stats::pchisq(J, df = overid_df, lower.tail = FALSE)
    list(stat = J, p = p, df = overid_df, test_name = "Sargan")
  } else {
    # HC/CL/HAC: Hansen J with robust omega
    .hansen_j_test(Z, X, y, residuals, weights, cluster_vec,
                   N, K, L, overid_df, dofminus = dofminus,
                   weight_type = weight_type,
                   kernel = kernel, bw = bw, time_index = time_index)
  }
}


# --------------------------------------------------------------------------
# .compute_ar_liml_overid
# --------------------------------------------------------------------------
#' Anderson-Rubin LIML overidentification statistics
#'
#' Computes the LR and linearized forms of the Anderson-Rubin
#' overidentification test for LIML estimation under IID errors.
#'
#' Both statistics are chi-squared(L - K) under the null that all
#' overidentifying restrictions are valid. The LR form uses
#' \eqn{(N - \text{dofminus}) \ln(\lambda)}{(N - dofminus) * ln(lambda)}
#' and the linearized form uses
#' \eqn{(N - \text{dofminus})(\lambda - 1)}{(N - dofminus) * (lambda - 1)}.
#'
#' @param lambda Numeric: LIML eigenvalue from `.fit_kclass()`.
#' @param N Integer: number of observations.
#' @param overid_df Integer: degree of overidentification (L - K).
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @return Named list with `lr_stat`, `lr_p`, `lin_stat`, `lin_p`, `df`,
#'   or a zero-stat placeholder when exactly identified (df == 0).
#' @keywords internal
.compute_ar_liml_overid <- function(lambda, N, overid_df, dofminus = 0L) {
  if (overid_df == 0L) {
    return(list(lr_stat = 0, lr_p = NA_real_,
                lin_stat = 0, lin_p = NA_real_,
                df = 0L))
  }

  # Guard: lambda must be >= 1 for a well-posed problem. Values < 1 indicate

  # numerical issues (rank-deficient concentration matrices, many-instrument
  # pathology). Stata produces missing silently; we warn and return NA.
  if (!is.finite(lambda) || lambda <= 0) {
    warning("LIML eigenvalue is non-positive (lambda = ",
            format(lambda, digits = 4),
            "); AR overidentification statistics not computed.",
            call. = FALSE)
    return(list(lr_stat = NA_real_, lr_p = NA_real_,
                lin_stat = NA_real_, lin_p = NA_real_,
                df = as.integer(overid_df)))
  }

  scale <- N - dofminus
  lr_stat  <- scale * log(lambda)
  lin_stat <- scale * (lambda - 1)

  lr_p  <- stats::pchisq(lr_stat,  df = overid_df, lower.tail = FALSE)
  lin_p <- stats::pchisq(lin_stat, df = overid_df, lower.tail = FALSE)

  list(lr_stat = lr_stat, lr_p = lr_p,
       lin_stat = lin_stat, lin_p = lin_p,
       df = as.integer(overid_df))
}
