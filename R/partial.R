# --------------------------------------------------------------------------
# FWL Projection / Partial Option (Ticket O1)
# --------------------------------------------------------------------------
# Partials out specified exogenous regressors from y, X, and Z via the
# Frisch-Waugh-Lovell theorem.  Mirroring Stata's s_partial() algorithm
# (ivreg2.ado lines 4993-5058).

#' Partial out exogenous regressors via FWL projection
#'
#' Projects y, X, and Z onto the orthogonal complement of the specified
#' exogenous regressor columns (P).  Uses three branches matching Stata:
#'
#' 1. **Variables + constant**: demean everything, then regress on P, take
#'    residuals.
#' 2. **Variables only (no constant)**: regress on P directly, take residuals.
#' 3. **Constant only**: demean (subtract weighted/unweighted column means).
#'
#' @param parsed A `parsed_formula` object from `.parse_formula()`.
#' @param partial_colnames Character vector of X column names to partial out.
#' @param partialcons Logical: whether the constant is being partialled.
#' @return A modified copy of `parsed` with:
#'   - `y`, `X`, `Z` projected onto M_P
#'   - Partialled columns removed from X and Z
#'   - `has_intercept` set to FALSE (intercept absorbed)
#' @keywords internal
.partial_out <- function(parsed, partial_colnames, partialcons) {
  y <- parsed$y
  X <- parsed$X
  Z <- parsed$Z
  w <- parsed$weights
  N <- length(y)

  # Identify columns to partial from X (and from Z, since exogenous
  # regressors appear in both)
  p_idx_x <- match(partial_colnames, colnames(X))
  p_idx_z <- if (!is.null(Z)) match(partial_colnames, colnames(Z)) else NULL

  # Build P matrix from X columns
  P <- X[, p_idx_x, drop = FALSE]

  # Remove partialled columns from X and Z
  keep_x <- setdiff(seq_len(ncol(X)), p_idx_x)
  X_remain <- X[, keep_x, drop = FALSE]

  Z_remain <- NULL
  if (!is.null(Z)) {
    keep_z <- setdiff(seq_len(ncol(Z)), p_idx_z)
    Z_remain <- Z[, keep_z, drop = FALSE]
  }

  # Remove intercept column (it's absorbed by partialling)
  icept_x <- which(colnames(X_remain) == "(Intercept)")
  if (length(icept_x) > 0L) {
    X_remain <- X_remain[, -icept_x, drop = FALSE]
  }
  if (!is.null(Z_remain)) {
    icept_z <- which(colnames(Z_remain) == "(Intercept)")
    if (length(icept_z) > 0L) {
      Z_remain <- Z_remain[, -icept_z, drop = FALSE]
    }
  }

  # --- Project via FWL ---
  # Combine all columns to project: y, remaining X, remaining Z (excluding
  # columns that are already in X)
  # We need to project y, every column of X_remain, and every column of Z_remain
  n_x <- ncol(X_remain)
  n_z <- if (!is.null(Z_remain)) ncol(Z_remain) else 0L

  # Columns unique to Z (excluded instruments + anything not in X_remain)
  if (n_z > 0L) {
    z_only_idx <- which(!colnames(Z_remain) %in% colnames(X_remain))
    Z_only <- Z_remain[, z_only_idx, drop = FALSE]
    n_z_only <- ncol(Z_only)
  } else {
    Z_only <- NULL
    n_z_only <- 0L
  }

  # Build Y matrix: [y, X_remain columns, Z_only columns]
  Y_all <- cbind(y, X_remain)
  if (n_z_only > 0L) {
    Y_all <- cbind(Y_all, Z_only)
  }

  n_p <- ncol(P)

  if (n_p == 0L && partialcons) {
    # Branch 3: constant only — just demean
    if (!is.null(w)) {
      col_means <- colSums(w * Y_all) / sum(w)
    } else {
      col_means <- colMeans(Y_all)
    }
    Y_resid <- sweep(Y_all, 2L, col_means)
  } else if (partialcons && n_p > 0L) {
    # Branch 1: variables + constant — demean first, then regress on demeaned P
    if (!is.null(w)) {
      # Weighted demeaning
      P_dm_means <- colSums(w * P) / sum(w)
      P_dm <- sweep(P, 2L, P_dm_means)
      Y_dm_means <- colSums(w * Y_all) / sum(w)
      Y_dm <- sweep(Y_all, 2L, Y_dm_means)
      # Weighted projection: use lm.wfit
      fit_proj <- stats::lm.wfit(P_dm, Y_dm, w = w)
    } else {
      P_dm <- scale(P, center = TRUE, scale = FALSE)
      Y_dm <- scale(Y_all, center = TRUE, scale = FALSE)
      fit_proj <- stats::lm.fit(P_dm, Y_dm)
    }
    Y_resid <- fit_proj$residuals
  } else {
    # Branch 2: variables only (no constant in model) — regress directly
    if (!is.null(w)) {
      fit_proj <- stats::lm.wfit(P, Y_all, w = w)
    } else {
      fit_proj <- stats::lm.fit(P, Y_all)
    }
    Y_resid <- fit_proj$residuals
  }

  # Ensure Y_resid is a matrix
  if (is.null(dim(Y_resid))) {
    Y_resid <- matrix(Y_resid, ncol = 1L)
  }

  # Extract projected components
  parsed$y <- Y_resid[, 1L]
  names(parsed$y) <- names(y)

  if (n_x > 0L) {
    parsed$X <- Y_resid[, 1L + seq_len(n_x), drop = FALSE]
    colnames(parsed$X) <- colnames(X_remain)
  } else {
    parsed$X <- matrix(nrow = N, ncol = 0L)
  }

  if (!is.null(Z_remain)) {
    # Reconstruct Z: shared columns from X_remain + Z-only projected columns
    shared_cols <- colnames(Z_remain)[colnames(Z_remain) %in% colnames(parsed$X)]
    if (n_z_only > 0L) {
      Z_only_proj <- Y_resid[, (1L + n_x) + seq_len(n_z_only), drop = FALSE]
      colnames(Z_only_proj) <- colnames(Z_only)
    } else {
      Z_only_proj <- NULL
    }
    # Rebuild Z: columns from projected X (shared) + projected Z-only
    if (length(shared_cols) > 0L && n_z_only > 0L) {
      parsed$Z <- cbind(parsed$X[, shared_cols, drop = FALSE], Z_only_proj)
    } else if (length(shared_cols) > 0L) {
      parsed$Z <- parsed$X[, shared_cols, drop = FALSE]
    } else if (n_z_only > 0L) {
      parsed$Z <- Z_only_proj
    } else {
      parsed$Z <- matrix(nrow = N, ncol = 0L)
    }
    # Ensure Z column order matches original Z_remain
    parsed$Z <- parsed$Z[, colnames(Z_remain), drop = FALSE]
  }

  # Intercept is absorbed
  parsed$has_intercept <- FALSE

  # Update dimensions
  parsed$K <- ncol(parsed$X)
  if (!is.null(parsed$Z)) {
    parsed$L <- ncol(parsed$Z)
  } else {
    parsed$L <- parsed$K
  }
  # Recompute K2 (exogenous regressors including intercept — but intercept gone)
  parsed$K2 <- parsed$K - parsed$K1
  # L1 = excluded instruments
  parsed$L1 <- parsed$L - parsed$K2

  parsed
}
