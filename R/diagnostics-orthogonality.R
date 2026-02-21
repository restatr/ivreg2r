# --------------------------------------------------------------------------
# .compute_orthog_test
# --------------------------------------------------------------------------
#' Instrument orthogonality test (C-statistic / difference-in-J)
#'
#' Tests H0: specified instruments satisfy orthogonality conditions.
#' Computes a difference-of-J (robust/cluster) or difference-of-Sargan (IID)
#' statistic by comparing the full model to a restricted model where the
#' tested instruments are removed.
#'
#' This is the mirror image of the endogeneity test: the **full** model
#' (more instruments) provides the S (Omega) matrix, and the restricted
#' model (fewer instruments) is the constrained model. The same-S-matrix
#' constraint guarantees C >= 0 (Hayashi 2000, p. 220).
#'
#' @param Z N x L instrument matrix (full model).
#' @param X N x K regressor matrix.
#' @param y N x 1 response vector.
#' @param residuals N x 1 residual vector from the full model.
#' @param rss Residual sum of squares from the full model.
#' @param weights Normalized weights (sum to N), or NULL.
#' @param cluster_vec Cluster membership vector, or NULL.
#' @param vcov_type Character: "iid", "HC0", "HC1", or "CL".
#' @param N Number of observations.
#' @param K Number of regressors.
#' @param L Number of instruments.
#' @param orthog_vars Character vector of instrument names to test.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @return Named list with `stat`, `p`, `df`, `test_name`, `tested_vars`,
#'   or NULL if this is not an IV model or orthog_vars is NULL.
#' @keywords internal
.compute_orthog_test <- function(Z, X, y, residuals, rss, weights,
                                  cluster_vec, vcov_type, N, K, L,
                                  orthog_vars, dofminus = 0L,
                                  weight_type = "aweight") {
  q <- length(orthog_vars)

  # --- Build restricted instrument matrix (remove tested columns) ---
  tested_idx <- match(orthog_vars, colnames(Z))
  keep_idx <- setdiff(seq_len(L), tested_idx)
  Z_r <- Z[, keep_idx, drop = FALSE]
  L_r <- ncol(Z_r)

  # Restricted model underidentified: L_r < K
  if (L_r < K) {
    return(list(stat = 0, p = NA_real_, df = 0L,
                test_name = "C (orthog)",
                tested_vars = orthog_vars))
  }

  # --- Compute Omega_full from full model's residuals ---
  if (vcov_type == "iid") {
    sigma_sq <- rss / (N - dofminus)
    if (is.null(weights)) {
      ZwZ <- crossprod(Z)
    } else {
      ZwZ <- crossprod(Z, weights * Z)
    }
    Omega_full <- sigma_sq * ZwZ / N
  } else {
    Omega_full <- .compute_omega(Z, residuals, weights, cluster_vec, N,
                                  dofminus = dofminus, weight_type = weight_type)
  }

  # --- J_full: J statistic of full model ---
  J_full <- .compute_j_with_omega(Z, X, y, Omega_full, weights, N)
  if (is.na(J_full)) {
    return(list(stat = NA_real_, p = NA_real_, df = as.integer(q),
                test_name = "C (orthog)",
                tested_vars = orthog_vars))
  }

  # --- Extract Omega submatrix for restricted model ---
  Omega_sub <- Omega_full[keep_idx, keep_idx, drop = FALSE]

  # --- J_r: J statistic of restricted model (using full model's S) ---
  J_r <- .compute_j_with_omega(Z_r, X, y, Omega_sub, weights, N)
  if (is.na(J_r)) {
    return(list(stat = NA_real_, p = NA_real_, df = as.integer(q),
                test_name = "C (orthog)",
                tested_vars = orthog_vars))
  }

  # --- C-statistic ---
  C <- J_full - J_r

  # Collinearity df check: overid_df - overid_df_r should equal q
  overid_df <- L - K
  overid_df_r <- L_r - K
  if (overid_df - overid_df_r != q) {
    return(list(stat = 0, p = NA_real_, df = 0L,
                test_name = "C (orthog)",
                tested_vars = orthog_vars))
  }

  p <- stats::pchisq(C, df = q, lower.tail = FALSE)

  list(
    stat = C,
    p = p,
    df = as.integer(q),
    test_name = "C (orthog)",
    tested_vars = orthog_vars
  )
}
