# --------------------------------------------------------------------------
# .fit_ols
# --------------------------------------------------------------------------
#' Fit OLS via lm.fit / lm.wfit (QR-based least squares)
#'
#' When weights are present, uses [lm.wfit()] which internally transforms by
#' `sqrt(w)` and returns unweighted residuals/fitted values with a QR of the
#' weighted design matrix. The bread is `(X'WX)^{-1}` automatically.
#'
#' @importFrom stats lm.fit lm.wfit
#' @param parsed A `parsed_formula` object from `.parse_formula()`.
#' @param small Logical: if `TRUE`, use `N-K` denominator for sigma;
#'   if `FALSE`, use `N`.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @return A named list with: `coefficients`, `residuals`, `fitted.values`,
#'   `vcov`, `sigma`, `df.residual`, `rank`, `r.squared`, `adj.r.squared`,
#'   `rss`, `bread`.
#' @keywords internal
.fit_ols <- function(parsed, small = FALSE, dofminus = 0L, sdofminus = 0L) {
  y <- parsed$y
  X <- parsed$X
  N <- parsed$N
  K <- parsed$K
  has_intercept <- parsed$has_intercept
  w <- parsed$weights

  # --- Fit via lm.fit / lm.wfit (QR-based least squares) ---
  if (is.null(w)) {
    lm_out <- lm.fit(X, y)
  } else {
    lm_out <- lm.wfit(X, y, w)
  }

  coef <- lm_out$coefficients
  fitted <- lm_out$fitted.values
  resid <- lm_out$residuals
  # lm.wfit returns unweighted residuals; RSS needs weighting
  rss <- if (is.null(w)) drop(crossprod(resid)) else sum(w * resid^2)

  # --- Sigma ---
  denom <- if (small) N - K - dofminus - sdofminus else N - dofminus
  sigma <- sqrt(rss / denom)

  # --- Bread: (X'X)^{-1} or (X'WX)^{-1} via stored QR ---
  # lm.wfit stores QR of sqrt(w)*X, so chol2inv gives (X'WX)^{-1}
  p <- lm_out$rank
  R <- lm_out$qr$qr[seq_len(p), seq_len(p), drop = FALSE]
  XtX_inv <- chol2inv(R)
  piv <- lm_out$qr$pivot
  piv_order <- order(piv[seq_len(p)])
  # Reorder rows/cols from pivot order back to original column order.
  # Names must use sort(piv[1:p]) to match: after reordering, row/col i
  # corresponds to original column sort(piv[1:p])[i], not piv[i].
  XtX_inv <- XtX_inv[piv_order, piv_order, drop = FALSE]
  colnames(XtX_inv) <- rownames(XtX_inv) <- colnames(X)[sort(piv[seq_len(p)])]

  # --- VCV ---
  V <- sigma^2 * XtX_inv
  colnames(V) <- rownames(V) <- colnames(XtX_inv)

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
    rank          = as.integer(lm_out$rank),
    r.squared     = r2,
    adj.r.squared = adj_r2,
    rss           = rss,
    r2u           = r2u,
    r2c           = r2c,
    mss           = mss,
    bread         = XtX_inv
  )
}
