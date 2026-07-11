# --------------------------------------------------------------------------
# Shared error helpers
# --------------------------------------------------------------------------
#' Stop: GMM Hessian singular (shared by all GMM estimation paths)
#' @noRd
.stop_gmm_hessian_singular <- function() {
  stop("The GMM Hessian matrix is singular, so the model is not identified. ",
       "Check for collinearity among the regressors and instruments.",
       call. = FALSE)
}

#' Stop: moment covariance rank-deficient (shared by all GMM estimation paths)
#' @noRd
.stop_gmm_omega_singular <- function() {
  stop("The estimated covariance matrix of the moment conditions is not of ",
       "full rank, so the optimal GMM weighting matrix is not unique. Too ",
       "few clusters or too many instruments is the usual cause; reduce the ",
       "instrument set or change the VCE.", call. = FALSE)
}

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
#' @noRd
.fit_gmm2s <- function(parsed, small = FALSE, dofminus = 0L, sdofminus = 0L,
                       omega_fn, omega_rank_bound = Inf) {
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
  # Structural gate first (platform-stable; rationale and Stata anchors in
  # .cluster_rank_bound). Stata's efficient-GMM path exits r(506) here
  # (s_egmm, ivreg2.ado:5436-5440).
  if (omega_rank_bound < L) {
    .stop_gmm_omega_singular()
  }
  R_chol <- tryCatch(chol(Omega), error = function(e) NULL)
  if (is.null(R_chol)) {
    if (qr(Omega)$rank < L) {
      .stop_gmm_omega_singular()
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
      .stop_gmm_hessian_singular()
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
    method        = "gmm2s",
    tss_u         = tss_u,
    tss_c         = tss_c
  )
}


# --------------------------------------------------------------------------
# .wmatrix_first_step_resid
# --------------------------------------------------------------------------
#' Compute residuals from W-weighted GMM first step
#'
#' Uses a user-supplied weighting matrix W to estimate beta via GMM, then
#' returns the residuals. Used in Path 3 (wmatrix + gmm2s) to feed the
#' omega_fn closure with W-step residuals for the efficient second step.
#'
#' @param parsed A `parsed_formula` object from `.parse_formula()`.
#' @param W L x L user-supplied weighting matrix.
#' @return Numeric vector of residuals from the W-weighted first step.
#' @noRd
.wmatrix_first_step_resid <- function(parsed, W) {
  X <- parsed$X
  Z <- parsed$Z
  y <- parsed$y
  N <- parsed$N
  K <- parsed$K
  w <- parsed$weights

  if (!is.null(w)) {
    QXZ <- crossprod(X, w * Z) / N
    QZy <- crossprod(Z, w * y) / N
  } else {
    QXZ <- crossprod(X, Z) / N
    QZy <- crossprod(Z, y) / N
  }

  H <- QXZ %*% W %*% t(QXZ)
  H <- (H + t(H)) / 2
  rhs <- QXZ %*% W %*% QZy

  R_H <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(R_H)) {
    if (qr(H)$rank < K)
      .stop_gmm_hessian_singular()
    beta <- drop(qr.solve(H, rhs))
  } else {
    beta <- drop(backsolve(R_H, forwardsolve(t(R_H), rhs)))
  }

  y - drop(X %*% beta)
}


# --------------------------------------------------------------------------
# .fit_gmm_wmatrix
# --------------------------------------------------------------------------
#' Fit GMM with user-supplied weighting matrix (inefficient GMM)
#'
#' Computes GMM estimates using user W as the weighting matrix. The VCV is
#' the full sandwich form (not the collapsed efficient form), and the J
#' statistic is computed via efficient 2-step re-estimation using the
#' estimated Omega (not the wmatrix-weighted residuals).
#'
#' Mirrors Stata's `s_iegmm()` (ivreg2.ado:5495-5563).
#'
#' @param parsed A `parsed_formula` object from `.parse_formula()`.
#' @param small Logical: if `TRUE`, use `N-K` denominator for sigma.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @param W L x L user-supplied weighting matrix.
#' @param omega_fn Function: closure that takes a residual vector and returns
#'   the L x L moment covariance matrix Omega.
#' @return A named list with the same fields as `.fit_gmm2s()` plus `method =
#'   "gmmw"`.
#' @noRd
.fit_gmm_wmatrix <- function(parsed, small = FALSE, dofminus = 0L,
                              sdofminus = 0L, W, omega_fn,
                              omega_rank_bound = Inf) {
  y <- parsed$y
  X <- parsed$X
  Z <- parsed$Z
  N <- parsed$N
  K <- parsed$K
  L <- parsed$L
  w <- parsed$weights
  overid_df <- parsed$overid_df
  has_intercept <- parsed$has_intercept

  # Get 2SLS bread and X_hat for diagnostics
  fit_2sls <- .fit_2sls(parsed, small = FALSE, dofminus = dofminus,
                         sdofminus = sdofminus)

  # --- Step 1: GMM estimation with user W ---
  if (!is.null(w)) {
    QXZ <- crossprod(X, w * Z) / N
    QZy <- crossprod(Z, w * y) / N
  } else {
    QXZ <- crossprod(X, Z) / N
    QZy <- crossprod(Z, y) / N
  }

  H <- QXZ %*% W %*% t(QXZ)
  H <- (H + t(H)) / 2   # force symmetry
  rhs <- QXZ %*% W %*% QZy

  # Rank check with Chol/QR fallback
  R_H <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(R_H)) {
    if (qr(H)$rank < K)
      .stop_gmm_hessian_singular()
    beta <- drop(qr.solve(H, rhs))
    H_inv <- qr.solve(H)
  } else {
    beta <- drop(backsolve(R_H, forwardsolve(t(R_H), rhs)))
    H_inv <- chol2inv(R_H)
  }
  names(beta) <- colnames(X)

  # --- Step 2: Residuals ---
  fitted <- drop(X %*% beta)
  names(fitted) <- names(y)
  resid <- y - fitted
  rss <- if (is.null(w)) drop(crossprod(resid)) else sum(w * resid^2)

  # --- Step 3: Omega from residuals ---
  Omega <- omega_fn(resid)

  # --- Step 4: Inefficient GMM VCV (Stata s_iegmm lines 5556-5557) ---
  # V = (1/N) * H^{-1} * QXZ * W * Omega * W * QXZ' * H^{-1}
  aux5 <- H_inv %*% QXZ %*% W    # K x L
  V <- (1 / N) * aux5 %*% Omega %*% t(aux5)
  V <- (V + t(V)) / 2
  colnames(V) <- rownames(V) <- colnames(X)

  # bread_gmm for VCV corrections (analogous to M_hess_inv / N in gmm2s)
  bread_gmm <- H_inv / N
  colnames(bread_gmm) <- rownames(bread_gmm) <- colnames(X)

  # --- Step 5: J statistic via efficient 2-step re-estimation ---
  if (overid_df > 0L) {
    # rank_bound: with user W the estimation itself needs no Omega inverse,
    # but the J does — Stata sets j to missing when rankS < iv1_ct
    # (ivreg2.ado:1793-1808; rationale in .cluster_rank_bound).
    j_stat <- .compute_j_with_omega(Z, X, y, Omega, w, N,
                                    rank_bound = omega_rank_bound)
    if (is.na(j_stat)) {
      j_p <- NA_real_
    } else {
      j_p <- stats::pchisq(j_stat, df = overid_df, lower.tail = FALSE)
    }
  } else {
    j_stat <- 0
    j_p <- NA_real_
  }

  # --- Step 6: sigma, R-squared, etc. (same as .fit_gmm2s) ---
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
    bread         = fit_2sls$bread,      # 2SLS bread for first-stage diagnostics
    bread_gmm     = bread_gmm,           # inefficient GMM bread
    X_hat         = fit_2sls$X_hat,
    j_stat        = j_stat,
    j_df          = as.integer(overid_df),
    j_p           = j_p,
    omega         = Omega,
    method        = "gmmw",
    tss_u         = tss_u,
    tss_c         = tss_c
  )
}


# --------------------------------------------------------------------------
# .cue_convergence
# --------------------------------------------------------------------------
#' Decide the CUE convergence code from the two optimizer results
#'
#' Nelder-Mead exits with code 10 (simplex degeneracy) when started at an
#' already-converged BFGS optimum on low-dimensional problems, so code 10 is
#' forgiven only when BFGS converged AND the restart provably moved the
#' OBJECTIVE by no more than 1e-6 relative (plus a 1e-10 absolute floor for
#' near-zero J; J >= 0 bounds the forgivable movement by J itself there).
#' The agreement tested is on the objective value, deliberately not on the
#' parameter vectors: in a flat CUE valley two well-separated betas can share
#' the same J, but that non-uniqueness is weak identification — surfaced by
#' the KP/AR identification diagnostics — not an optimization failure, and
#' the reported coefficients are always the Nelder-Mead point regardless of
#' this flag. The 1e-6 bound is derived, not fitted: BFGS's finite-difference
#' gradients leave it short of the optimum by O(sqrt(machine eps)) ~ 1.5e-8
#' relative on the objective, so genuine refinement noise sits two orders
#' below the bound, while a genuine move to a different basin shifts the
#' objective by many orders more (on the Griliches help-file CUE example,
#' H22, the gap between the degenerate basin and a mid-valley stall point is
#' order 1e0 in J). The additional absolute ceiling of 1 closes the huge-J
#' scale gap: at J ~ 1e9 a purely relative bound would forgive
#' chi-square-scale movement, and movement that could visibly change a J
#' p-value is never forgiven. isTRUE() keeps a non-finite objective from
#' being forgiven. Code 1 (maxit) is never forgiven: an optimizer still
#' improving at its iteration cap is genuinely unconverged. No bundled
#' fixture reaches the non-convergence warning (H22 converges and exercises
#' the forgiveness path), so this decision logic is pinned by direct unit
#' tests in test-cue.R; only the one-line warning emission in .fit_gmm_cue
#' is fixture-uncovered.
#'
#' @param opt_bfgs `stats::optim()` result of the BFGS stage (uses
#'   `$convergence`, `$value`).
#' @param opt_nm `stats::optim()` result of the Nelder-Mead restart (uses
#'   `$convergence`, `$value`).
#' @return Integer convergence code: `0L` for convergence (including forgiven
#'   code-10 exits), otherwise the Nelder-Mead code.
#' @noRd
.cue_convergence <- function(opt_bfgs, opt_nm) {
  improvement <- opt_bfgs$value - opt_nm$value
  improvement_small <- isTRUE(
    improvement <= abs(opt_bfgs$value) * 1e-6 + 1e-10 && improvement <= 1
  )
  if (opt_nm$convergence == 0L ||
      (opt_nm$convergence == 10L &&
       opt_bfgs$convergence == 0L && improvement_small)) {
    0L
  } else {
    opt_nm$convergence
  }
}


# --------------------------------------------------------------------------
# .fit_cue
# --------------------------------------------------------------------------
#' Fit Continuously Updated GMM Estimator (CUE)
#'
#' Minimizes the CUE objective `J(b) = N * gbar(b)' * Omega(b)^{-1} * gbar(b)`
#' where both moment conditions `gbar` and moment covariance `Omega` are
#' re-evaluated at each candidate beta. Uses BFGS + Nelder-Mead for
#' optimization. Starting values are 2SLS (IID) or two-step efficient GMM
#' (non-IID), matching Stata ivreg2.ado lines 5918-5948.
#'
#' When `b0` is supplied, evaluates the CUE objective at the fixed parameter
#' vector without optimization (Stata's `b0()` option).
#'
#' @param parsed A `parsed_formula` object from `.parse_formula()`.
#' @param small Logical: if `TRUE`, use `N-K` denominator for sigma.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @param omega_fn Function: closure that takes a residual vector and returns
#'   the L x L moment covariance matrix Omega. Built in `ivreg2()`.
#' @param b0 Optional numeric vector: fixed parameter vector for J evaluation
#'   (no optimization). NULL for standard CUE optimization.
#' @param iid Logical: if `TRUE` (default), use 2SLS starting values. If
#'   `FALSE` (non-IID VCE: robust, cluster, HAC, AC), use two-step efficient
#'   GMM starting values, matching Stata's behavior (ivreg2.ado changelog
#'   line 6596).
#' @return A named list with the same fields as `.fit_gmm2s()` plus
#'   `convergence` (integer) and `cue_message` (string). `method = "cue"`.
#' @noRd
.fit_cue <- function(parsed, small = FALSE, dofminus = 0L, sdofminus = 0L,
                     omega_fn, b0 = NULL, iid = TRUE,
                     omega_rank_bound = Inf) {
  y <- parsed$y
  X <- parsed$X
  Z <- parsed$Z
  N <- parsed$N
  K <- parsed$K
  L <- parsed$L
  has_intercept <- parsed$has_intercept
  w <- parsed$weights
  overid_df <- parsed$overid_df

  # --- Step 1: Starting values or fixed beta ---
  if (is.null(b0)) {
    # When VCE is non-IID (robust, cluster, HAC, AC), use two-step efficient
    # GMM starting values instead of plain 2SLS. The efficient GMM estimates
    # are already in the right neighborhood for the robust CUE objective.
    # Matches Stata's documented behavior (ivreg2.ado changelog line 6596:
    # "Fixed cue bug - was not starting with robust 2-step gmm if
    # robust/cluster").
    if (iid) {
      fit_init <- .fit_2sls(parsed, small = FALSE, dofminus = dofminus,
                            sdofminus = sdofminus)
    } else {
      fit_init <- .fit_gmm2s(parsed, small = FALSE, dofminus = dofminus,
                             sdofminus = sdofminus, omega_fn = omega_fn,
                             omega_rank_bound = omega_rank_bound)
    }
    beta_init <- unname(fit_init$coefficients)

    if (overid_df == 0L) {
      # Exactly-identified: CUE = 2SLS, no optimization needed.
      # The CUE objective is flat near the minimum (J=0 at any consistent
      # estimator), so numerical optimizers report false non-convergence.
      beta <- beta_init
      convergence <- 0L
      cue_message <- "exactly identified (CUE = 2SLS)"
      j_from_opt <- NULL
    } else {
      # --- Step 2: CUE objective function ---
      # Precompute weighted cross-products that don't depend on beta
      if (!is.null(w)) {
        wZ <- w * Z
      } else {
        wZ <- Z
      }

      cue_obj <- function(beta) {
        e <- y - drop(X %*% beta)
        Omega <- omega_fn(e)

        # gbar = Z'(w*e) / N
        gbar <- crossprod(Z, if (!is.null(w)) w * e else e) / N

        # Solve Omega^{-1} * gbar via Cholesky with QR fallback
        R_chol <- tryCatch(chol(Omega), error = function(e) NULL)
        if (is.null(R_chol)) {
          # Check rank — if singular, return penalty
          if (qr(Omega)$rank < L) return(.Machine$double.xmax)
          Omega_inv_gbar <- tryCatch(qr.solve(Omega, gbar),
                                     error = function(e) NULL)
          if (is.null(Omega_inv_gbar)) return(.Machine$double.xmax)
        } else {
          Omega_inv_gbar <- backsolve(R_chol,
                                      forwardsolve(t(R_chol), gbar))
        }

        # J = N * gbar' * Omega^{-1} * gbar
        drop(N * crossprod(gbar, Omega_inv_gbar))
      }

      # --- Step 3: Optimize ---
      # Two-step strategy: BFGS for fast initial convergence, then Nelder-Mead
      # for robust refinement. nlminb struggles with the CUE ratio objective
      # due to noisy finite-difference gradients; this combination matches
      # Stata's Mata optimize() results reliably.
      # Per-coordinate scaling: optim's Nelder-Mead builds its initial
      # simplex with a SINGLE scalar step — 0.1 * max|scaled coefficient|
      # applied to every coordinate (nmmin C source) — and BFGS takes its
      # finite-difference gradient steps (ndeps) in scaled units, so under
      # the default parscale = 1 a coefficient vector of mixed magnitudes
      # gets a distorted simplex and gradient geometry. Scaling by the
      # starting values makes every scaled coordinate O(1), so both
      # optimizers step proportionally to each coefficient's own magnitude.
      # The floor at 1e-4 of the largest |coefficient| covers start
      # coefficients that carry no usable scale of their own — exact zeros,
      # and coordinates accidentally started many orders below the rest,
      # which an unfloored parscale would near-freeze (BFGS reads a ~zero
      # gradient through a ~zero scale). The floor is inert on every
      # benchmarked model (observed min/max coefficient ratios are all
      # >= 3.2e-4) and caps the scale spread the optimizers must handle at
      # 1e4. An all-zero start (no scale information at all) falls back to
      # optim's default scale of 1. Non-finite starts are NOT handled here:
      # optim rejects them identically with or without parscale
      # ("non-finite value supplied by optim", verified), so a guard would
      # be dead code.
      parscale <- abs(beta_init)
      scale_ref <- max(parscale)
      if (scale_ref == 0) {
        parscale <- rep(1, length(parscale))
      } else {
        parscale <- pmax(parscale, 1e-4 * scale_ref)
      }
      opt_bfgs <- stats::optim(par = beta_init, fn = cue_obj, method = "BFGS",
                               control = list(maxit = 500L, reltol = 1e-12,
                                              parscale = parscale))
      opt_nm <- stats::optim(par = opt_bfgs$par, fn = cue_obj,
                             method = "Nelder-Mead",
                             control = list(maxit = 5000L * K, reltol = 1e-14,
                                            parscale = parscale))
      # Nelder-Mead starts at the BFGS point, which is a vertex of its initial simplex, and its tracked best vertex is monotonically non-increasing, so opt_nm can never be worse than opt_bfgs — take it unconditionally.
      beta <- opt_nm$par
      j_from_opt <- opt_nm$value

      convergence <- .cue_convergence(opt_bfgs, opt_nm)
      cue_message <- if (convergence == 0L) {
        "relative convergence (4)"
      } else {
        "Nelder-Mead did not converge"
      }

      if (convergence != 0L) {
        warning("CUE optimization did not converge, so the reported CUE ",
                "estimates may be unreliable. A flat or ill-behaved GMM ",
                "objective is the usual cause, often from weak instruments; ",
                "check the weak-identification diagnostics and compare against ",
                'two-step GMM (`method = "gmm2s"`).', call. = FALSE)
      }
    }

  } else {
    # b0 path: fixed beta, no optimization
    beta <- unname(b0)
    convergence <- 0L
    cue_message <- "b0 evaluation (no optimization)"
    j_from_opt <- NULL

    # Still need 2SLS for bread and X_hat
    fit_init <- .fit_2sls(parsed, small = FALSE, dofminus = dofminus,
                          sdofminus = sdofminus)
  }

  names(beta) <- colnames(X)

  # --- Step 4: Recompute Omega at final beta (Stata line 5946-5948) ---
  fitted <- drop(X %*% beta)
  names(fitted) <- names(y)
  resid <- y - fitted
  Omega <- omega_fn(resid)

  # --- Step 5: Omega rank check (Stata line 5972-5979) ---
  # Structural gate first, as in .fit_gmm2s step 3 (rationale in
  # .cluster_rank_bound). For non-iid CUE the GMM2S init above already
  # errors on this condition; this belt covers the b0 path, which skips the
  # init's Omega inversion.
  if (omega_rank_bound < L) {
    .stop_gmm_omega_singular()
  }
  R_chol <- tryCatch(chol(Omega), error = function(e) NULL)
  if (is.null(R_chol)) {
    if (qr(Omega)$rank < L) {
      .stop_gmm_omega_singular()
    }
    R_chol <- NULL
  }

  omega_inv <- function(b) {
    if (!is.null(R_chol)) {
      backsolve(R_chol, forwardsolve(t(R_chol), b))
    } else {
      qr.solve(Omega, b)
    }
  }

  # --- Step 6: VCV (same formula as .fit_gmm2s) ---
  if (!is.null(w)) {
    QXZ <- crossprod(X, w * Z) / N
  } else {
    QXZ <- crossprod(X, Z) / N
  }

  aux1 <- omega_inv(t(QXZ))   # L x K: Omega^{-1} QXZ'
  M_hess <- QXZ %*% aux1      # K x K: GMM Hessian
  M_hess <- (M_hess + t(M_hess)) / 2  # force symmetry

  R_M <- tryCatch(chol(M_hess), error = function(e) NULL)
  if (is.null(R_M)) {
    if (qr(M_hess)$rank < K) {
      .stop_gmm_hessian_singular()
    }
    M_hess_inv <- qr.solve(M_hess)
  } else {
    M_hess_inv <- chol2inv(R_M)
  }

  # VCV rank check (Stata line 5983-5988)
  diag_V <- diag(M_hess_inv) / N
  if (any(diag_V <= 0)) {
    stop("The estimated variance-covariance matrix has zero or negative ",
         "diagonal elements, so the model may not be identified. Check for ",
         "collinearity among the regressors and instruments.", call. = FALSE)
  }

  V <- M_hess_inv / N
  colnames(V) <- rownames(V) <- colnames(X)

  bread_gmm <- M_hess_inv / N
  colnames(bread_gmm) <- rownames(bread_gmm) <- colnames(X)

  # --- Step 7: RSS ---
  rss <- if (is.null(w)) drop(crossprod(resid)) else sum(w * resid^2)

  # --- Step 8: J statistic ---
  if (overid_df > 0L || !is.null(b0)) {
    # For overidentified CUE without b0: use optimizer J directly
    # For b0: compute J at b0 (always, even if exactly identified)
    if (!is.null(b0)) {
      # Compute J at b0
      gbar <- crossprod(Z, if (!is.null(w)) w * resid else resid) / N
      Omega_inv_gbar <- omega_inv(gbar)
      j_stat <- drop(N * crossprod(gbar, Omega_inv_gbar))
      j_p <- if (overid_df > 0L) {
        stats::pchisq(j_stat, df = overid_df, lower.tail = FALSE)
      } else {
        NA_real_
      }
    } else if (overid_df > 0L) {
      # Overidentified CUE: use J from optimizer
      j_stat <- j_from_opt
      j_p <- stats::pchisq(j_stat, df = overid_df, lower.tail = FALSE)
    } else {
      # Exactly-identified CUE without b0: force J = 0
      j_stat <- 0
      j_p <- NA_real_
    }
  } else {
    # Exactly-identified without b0
    j_stat <- 0
    j_p <- NA_real_
  }

  # --- Step 9: sigma, R-squared, etc. (same as .fit_gmm2s) ---
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
    bread         = fit_init$bread,   # 2SLS bread for first-stage diagnostics
    bread_gmm     = bread_gmm,       # CUE bread for VCV / diagnostics
    X_hat         = fit_init$X_hat,   # from 2SLS init
    j_stat        = j_stat,
    j_df          = as.integer(overid_df),
    j_p           = j_p,
    omega         = Omega,
    convergence   = as.integer(convergence),
    cue_message   = cue_message,
    method        = "cue",
    tss_u         = tss_u,
    tss_c         = tss_c
  )
}
