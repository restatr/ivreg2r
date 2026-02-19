# --------------------------------------------------------------------------
# .fit_2sls
# --------------------------------------------------------------------------
#' Fit 2SLS via two calls to lm.fit / lm.wfit
#'
#' Stage 1 projects regressors X onto the instrument space Z.
#' Stage 2 regresses y on the projected regressors X_hat.
#' Residuals and fitted values use the *original* X, not X_hat.
#'
#' When weights are present, uses [lm.wfit()] for both stages. `lm.wfit`
#' returns fitted values on the original (unweighted) scale and a QR of
#' the weighted design matrix, so the bread is `(X_hat'WX_hat)^{-1}`.
#'
#' @importFrom stats lm.fit lm.wfit
#' @param parsed A `parsed_formula` object from `.parse_formula()`.
#' @param small Logical: if `TRUE`, use `N-K` denominator for sigma;
#'   if `FALSE`, use `N`.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @return A named list with: `coefficients`, `residuals`, `fitted.values`,
#'   `vcov`, `sigma`, `df.residual`, `rank`, `r.squared`, `adj.r.squared`,
#'   `rss`, `bread`, `X_hat`.
#' @keywords internal
.fit_2sls <- function(parsed, small = FALSE, dofminus = 0L, sdofminus = 0L) {
  y <- parsed$y
  X <- parsed$X
  Z <- parsed$Z
  N <- parsed$N
  K <- parsed$K
  has_intercept <- parsed$has_intercept
  w <- parsed$weights

  # --- Stage 1: project X onto column space of Z ---
  if (is.null(w)) {
    lm_Z <- lm.fit(Z, X)
  } else {
    lm_Z <- lm.wfit(Z, X, w)
  }
  X_hat <- as.matrix(lm_Z$fitted.values)
  colnames(X_hat) <- colnames(X)

  # --- Stage 2: regress y on projected regressors ---
  if (is.null(w)) {
    lm_XZ <- lm.fit(X_hat, y)
  } else {
    lm_XZ <- lm.wfit(X_hat, y, w)
  }

  if (lm_XZ$rank < K) {
    stop("Projected regressor matrix X_hat is rank-deficient (rank ",
         lm_XZ$rank, " < ", K, " regressors). ",
         "Instruments may not identify the endogenous regressors.",
         call. = FALSE)
  }

  coef <- lm_XZ$coefficients

  # --- Residuals use original X, not X_hat ---
  fitted <- drop(X %*% coef)
  names(fitted) <- names(y)
  resid <- y - fitted
  # lm.wfit returns unweighted residuals; RSS needs weighting
  rss <- if (is.null(w)) drop(crossprod(resid)) else sum(w * resid^2)

  # --- Bread: (X_hat'WX_hat)^{-1} via stored QR from stage 2 ---
  # lm.wfit stores QR of sqrt(w)*X_hat, so chol2inv gives (X_hat'WX_hat)^{-1}
  p <- lm_XZ$rank
  R <- lm_XZ$qr$qr[seq_len(p), seq_len(p), drop = FALSE]
  XtPX_inv <- chol2inv(R)
  piv <- lm_XZ$qr$pivot
  piv_order <- order(piv[seq_len(p)])
  XtPX_inv <- XtPX_inv[piv_order, piv_order, drop = FALSE]
  colnames(XtPX_inv) <- rownames(XtPX_inv) <- colnames(X)[sort(piv[seq_len(p)])]

  # --- Sigma ---
  denom <- if (small) N - K - dofminus - sdofminus else N - dofminus
  sigma <- sqrt(rss / denom)

  # --- VCV ---
  V <- sigma^2 * XtPX_inv
  colnames(V) <- rownames(V) <- colnames(XtPX_inv)

  # --- R-squared (centered and uncentered) ---
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

  # --- Adjusted R-squared ---
  adj_r2 <- if (has_intercept) {
    1 - (1 - r2) * (N - 1) / (N - K - dofminus - sdofminus)
  } else {
    1 - (1 - r2) * N / (N - K - dofminus - sdofminus)
  }

  list(
    coefficients  = coef,
    residuals     = resid,
    fitted.values = fitted,
    vcov          = V,
    sigma         = sigma,
    df.residual   = as.integer(N - K - dofminus - sdofminus),
    rank          = as.integer(lm_XZ$rank),
    r.squared     = r2,
    adj.r.squared = adj_r2,
    rss           = rss,
    r2u           = r2u,
    r2c           = r2c,
    mss           = mss,
    bread         = XtPX_inv,
    X_hat         = X_hat
  )
}
