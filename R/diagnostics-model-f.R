# --------------------------------------------------------------------------
# Model F-test (Ticket F1)
# --------------------------------------------------------------------------
# Wald test of H0: all slope coefficients = 0, reported as F.
# Works for both OLS and IV; uses the VCV already computed by the main
# ivreg2() pipeline (IID, HC, or cluster-robust).


# --------------------------------------------------------------------------
# .compute_model_f
# --------------------------------------------------------------------------
#' Model F-test (all slopes = 0)
#'
#' Computes the Wald F-statistic for the joint null that all slope
#' coefficients equal zero. The intercept (if present) is excluded from
#' the test.
#'
#' @param coefficients Named numeric vector of coefficient estimates.
#' @param vcov Variance-covariance matrix of coefficients.
#' @param N Integer: number of observations.
#' @param K Integer: number of regressors (including intercept if present).
#' @param has_intercept Logical: does the model include an intercept?
#' @param vcov_type Character: `"iid"`, `"HC0"`, `"HC1"`, or `"CL"`.
#' @param small Logical: whether small-sample corrections were applied.
#' @param M Integer or NULL: number of clusters (only for `vcov_type = "CL"`).
#' @return Named list with `model_f`, `model_f_p`, `model_f_df1`,
#'   `model_f_df2` (all `NA` when df1 = 0 or VCV is not positive definite).
#' @note The model F-statistic is always reported as an F-distributed
#'   statistic. For non-IID VCE types, the underlying Wald chi-squared is
#'   converted to F using VCE-branch-specific degrees of freedom adjustments.
#'   This matches Stata's \code{ivreg2} behavior.
#' @keywords internal
.compute_model_f <- function(coefficients, vcov, N, K,
                              has_intercept, vcov_type, small,
                              M = NULL) {

  # --- Identify slope coefficients (exclude intercept) ---
  coef_names <- names(coefficients)
  intercept_idx <- if (has_intercept) {
    match("(Intercept)", coef_names, nomatch = 0L)
  } else {
    0L
  }
  slope_idx <- setdiff(seq_along(coefficients), intercept_idx)
  df1 <- length(slope_idx)

  # Edge case: intercept-only model
  if (df1 == 0L) {
    return(list(model_f = NA_real_, model_f_p = NA_real_,
                model_f_df1 = NA_integer_, model_f_df2 = NA_integer_))
  }

  # --- Extract slope subvector and submatrix ---
  beta_s <- unname(coefficients[slope_idx])
  V_s <- vcov[slope_idx, slope_idx, drop = FALSE]

  # --- Wald chi-squared ---
  # Try Cholesky for full-rank VCV; fall back to Gauss-Jordan sweep with
  # partial pivoting for rank-deficient cases (e.g., M < K clusters).
  # The fallback matches Stata's `syminv()` algorithm used by `test`.
  chi2 <- tryCatch({
    R_chol <- chol(V_s)
    z <- forwardsolve(t(R_chol), beta_s)
    drop(crossprod(z))
  }, error = function(e) {
    G <- .syminv_sweep(V_s)
    if (is.null(G)) return(NA_real_)
    drop(crossprod(beta_s, G %*% beta_s))
  })

  if (is.na(chi2) || !is.finite(chi2) || chi2 < 0) {
    df2 <- if (!is.null(M)) as.integer(M - 1L) else as.integer(N - K)
    return(list(model_f = NA_real_, model_f_p = NA_real_,
                model_f_df1 = as.integer(df1), model_f_df2 = df2))
  }

  # --- Convert chi-squared to F ---
  # The Wald chi2 is computed from the VCV that may or may not include
  # finite-sample corrections. The F conversion must undo those corrections
  # so the F-statistic is invariant to whether corrections were applied.
  #
  # Three branches:
  #   1. Corrected VCV (iid+small, HC1, CL+small): F = chi2 / df1
  #   2. Uncorrected, no cluster (iid+!small, HC0): F = chi2/df1 * (N-K)/N
  #   3. Uncorrected, cluster (CL+!small): F = chi2/df1 * (M-1)/M * (N-K)/(N-1)

  is_corrected <- (vcov_type == "iid" && small) ||
                  (vcov_type == "HC1") ||
                  (vcov_type == "CL" && small)

  if (is_corrected) {
    f_stat <- chi2 / df1
  } else if (is.null(M)) {
    # Uncorrected, no cluster
    f_stat <- chi2 / df1 * (N - K) / N
  } else {
    # Uncorrected, cluster
    f_stat <- chi2 / df1 * (M - 1) / M * (N - K) / (N - 1)
  }

  # --- Degrees of freedom ---
  df2 <- if (!is.null(M)) as.integer(M - 1L) else as.integer(N - K)

  # --- p-value ---
  f_p <- stats::pf(f_stat, df1 = df1, df2 = df2, lower.tail = FALSE)

  list(
    model_f     = f_stat,
    model_f_p   = f_p,
    model_f_df1 = as.integer(df1),
    model_f_df2 = df2
  )
}


# --------------------------------------------------------------------------
# .syminv_sweep
# --------------------------------------------------------------------------
#' Generalized inverse via Gauss-Jordan elimination with partial pivoting
#'
#' Matches the behavior of Stata's `syminv()`: at each step, the largest
#' remaining diagonal element is selected as the pivot (partial pivoting).
#' Rows/columns whose pivot falls below a tolerance are zeroed out.
#' This differs from the Moore-Penrose pseudoinverse for rank-deficient
#' matrices and matches how Stata's `test` command computes Wald statistics.
#'
#' @param A Symmetric n x n matrix (typically a VCV submatrix).
#' @return The generalized inverse matrix, or `NULL` if no pivots succeed.
#' @keywords internal
.syminv_sweep <- function(A) {
  n <- nrow(A)
  # Augmented matrix [A | I] — Gauss-Jordan on left half, inverse in right
  AI <- cbind(A, diag(n))
  tol <- max(abs(diag(A))) * sqrt(.Machine$double.eps)
  pivoted <- logical(n)

  for (step in seq_len(n)) {
    # Partial pivoting: largest *current* diagonal among remaining indices
    remaining <- which(!pivoted)
    if (length(remaining) == 0L) break
    diag_vals <- abs(AI[remaining, remaining, drop = FALSE])
    diag_vals <- diag(diag_vals)
    best <- remaining[which.max(diag_vals)]

    if (abs(AI[best, best]) <= tol) break

    # Scale pivot row
    AI[best, ] <- AI[best, ] / AI[best, best]
    # Eliminate from all other rows
    for (j in setdiff(seq_len(n), best)) {
      AI[j, ] <- AI[j, ] - AI[j, best] * AI[best, ]
    }
    pivoted[best] <- TRUE
  }

  if (!any(pivoted)) return(NULL)

  # Extract inverse from right half
  G <- AI[, (n + 1L):(2L * n), drop = FALSE]

  # Zero out un-pivoted rows/columns
  not_pivoted <- which(!pivoted)
  if (length(not_pivoted) > 0L) {
    G[not_pivoted, ] <- 0
    G[, not_pivoted] <- 0
  }

  # Enforce symmetry
  (G + t(G)) / 2
}
