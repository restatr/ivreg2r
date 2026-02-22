# --------------------------------------------------------------------------
# .fit_gmm2s
# --------------------------------------------------------------------------
#' Fit two-step efficient GMM
#'
#' Computes two-step efficient GMM estimates using the inverse of the
#' moment covariance matrix (Omega) as the optimal weighting matrix.
#' Mirrors Stata's `s_egmm()` (ivreg2.ado:5384-5490).
#'
#' Algorithm:
#' 1. First step: call [.fit_2sls()] to get IV coefficients and residuals.
#' 2. Compute Omega from step-1 residuals via `omega_fn`.
#' 3. Rank-check Omega; error if singular.
#' 4. Efficient GMM re-estimation using optimal weighting matrix W = Omega inverse.
#' 5. Compute J statistic, R-squared, sigma, etc.
#'
#' The efficient GMM VCV is V = N * (X'Z Omega_inv Z'X)_inv — the
#' sandwich collapses, so no separate meat computation is needed.
#'
#' @param parsed A `parsed_formula` object from `.parse_formula()`.
#' @param small Logical: if `TRUE`, use `N-K` denominator for sigma.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @param omega_fn Function: closure that takes a residual vector and returns
#'   the L x L moment covariance matrix Omega. Built in `ivreg2()`.
#' @return A named list with: `coefficients`, `residuals`, `fitted.values`,
#'   `vcov`, `sigma`, `df.residual`, `rank`, `r.squared`, `adj.r.squared`,
#'   `rss`, `r2u`, `r2c`, `mss`, `bread`, `bread_gmm`, `X_hat`,
#'   `j_stat`, `j_df`, `j_p`, `omega`, `method`.
#' @keywords internal
.fit_gmm2s <- function(parsed, small = FALSE, dofminus = 0L, sdofminus = 0L,
                       omega_fn) {
  y <- parsed$y
  X <- parsed$X
  Z <- parsed$Z
  N <- parsed$N
  K <- parsed$K
  L <- parsed$L
  has_intercept <- parsed$has_intercept
  w <- parsed$weights
  overid_df <- parsed$overid_df

  # --- Step 1: first-step 2SLS ---
  fit_1 <- .fit_2sls(parsed, small = FALSE, dofminus = dofminus,
                     sdofminus = sdofminus)

  # --- Step 2: compute Omega from step-1 residuals ---
  Omega <- omega_fn(fit_1$residuals)

  # --- Step 3: rank-check Omega ---
  R_chol <- tryCatch(chol(Omega), error = function(e) NULL)
  if (is.null(R_chol)) {
    if (qr(Omega)$rank < L) {
      stop("estimated covariance matrix of moment conditions not of full rank;\n",
           "optimal GMM weighting matrix not unique.", call. = FALSE)
    }
    # Omega is full rank but not PD — use QR path
    R_chol <- NULL
  }

  # Omega^{-1} solver
  omega_inv <- function(b) {
    if (!is.null(R_chol)) {
      backsolve(R_chol, forwardsolve(t(R_chol), b))
    } else {
      qr.solve(Omega, b)
    }
  }

  # --- Step 4: efficient GMM re-estimation ---
  # QXZ = X'WZ / N  (K x L),  QZy = Z'Wy / N  (L x 1)
  if (!is.null(w)) {
    QXZ <- crossprod(X, w * Z) / N
    QZy <- crossprod(Z, w * y) / N
  } else {
    QXZ <- crossprod(X, Z) / N
    QZy <- crossprod(Z, y) / N
  }

  aux1 <- omega_inv(t(QXZ))   # L x K: Omega^{-1} QXZ'
  aux2 <- omega_inv(QZy)      # L x 1: Omega^{-1} QZy

  M_hess <- QXZ %*% aux1      # K x K: GMM Hessian
  M_hess <- (M_hess + t(M_hess)) / 2  # force symmetry

  # Solve for beta: M_hess * beta = QXZ * Omega^{-1} * QZy
  rhs <- QXZ %*% aux2
  R_M <- tryCatch(chol(M_hess), error = function(e) NULL)
  if (is.null(R_M)) {
    if (qr(M_hess)$rank < K) {
      stop("GMM Hessian matrix is singular; model not identified.",
           call. = FALSE)
    }
    beta <- drop(qr.solve(M_hess, rhs))
  } else {
    beta <- drop(backsolve(R_M, forwardsolve(t(R_M), rhs)))
  }
  names(beta) <- colnames(X)

  # --- Step 5: GMM VCV (base, pre-small-sample) ---
  # V_base = (1/N) * M_hess^{-1} = N * (X'Z Omega^{-1} Z'X)^{-1}
  if (!is.null(R_M)) {
    M_hess_inv <- chol2inv(R_M)
  } else {
    M_hess_inv <- qr.solve(M_hess)
  }
  V <- M_hess_inv / N
  colnames(V) <- rownames(V) <- colnames(X)

  # bread_gmm for diagnostics
  bread_gmm <- M_hess_inv / N
  colnames(bread_gmm) <- rownames(bread_gmm) <- colnames(X)

  # --- Step 6: residuals and RSS ---
  fitted <- drop(X %*% beta)
  names(fitted) <- names(y)
  resid <- y - fitted
  rss <- if (is.null(w)) drop(crossprod(resid)) else sum(w * resid^2)

  # --- Step 7: J statistic ---
  if (overid_df > 0L) {
    if (!is.null(w)) {
      gbar <- crossprod(Z, w * resid) / N
    } else {
      gbar <- crossprod(Z, resid) / N
    }
    Omega_inv_gbar <- omega_inv(gbar)
    j_stat <- drop(N * crossprod(gbar, Omega_inv_gbar))
    j_p <- stats::pchisq(j_stat, df = overid_df, lower.tail = FALSE)
  } else {
    j_stat <- 0
    j_p <- NA_real_
  }

  # --- Step 8: sigma, R-squared, etc. ---
  denom <- if (small) N - K - dofminus - sdofminus else N - dofminus
  sigma <- sqrt(rss / denom)

  if (is.null(w)) {
    tss_c <- sum((y - mean(y))^2)
    tss_u <- sum(y^2)
  } else {
    wmean <- sum(w * y) / sum(w)
    tss_c <- sum(w * (y - wmean)^2)
    tss_u <- sum(w * y^2)
  }
  r2c <- 1 - rss / tss_c
  r2u <- 1 - rss / tss_u
  r2  <- if (has_intercept) r2c else r2u
  tss <- if (has_intercept) tss_c else tss_u
  mss <- tss - rss

  adj_r2 <- if (has_intercept) {
    1 - (1 - r2) * (N - 1) / (N - K - dofminus - sdofminus)
  } else {
    1 - (1 - r2) * N / (N - K - dofminus - sdofminus)
  }

  list(
    coefficients  = beta,
    residuals     = resid,
    fitted.values = fitted,
    vcov          = V,
    sigma         = sigma,
    df.residual   = as.integer(N - K - dofminus - sdofminus),
    rank          = as.integer(K),
    r.squared     = r2,
    adj.r.squared = adj_r2,
    rss           = rss,
    r2u           = r2u,
    r2c           = r2c,
    mss           = mss,
    bread         = fit_1$bread,     # 2SLS bread for first-stage diagnostics
    bread_gmm     = bread_gmm,      # GMM bread for VCV / diagnostics
    X_hat         = fit_1$X_hat,     # from step 1
    j_stat        = j_stat,
    j_df          = as.integer(overid_df),
    j_p           = j_p,
    omega         = Omega,
    method        = "gmm2s"
  )
}
