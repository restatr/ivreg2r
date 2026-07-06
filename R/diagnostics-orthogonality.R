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
#' @param vcov_type Character: "iid", "robust", "HAC", "AC", or "CL".
#' @param N Number of observations.
#' @param K Number of regressors.
#' @param L Number of instruments.
#' @param orthog_vars Character vector of instrument names to test.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param omega Pre-computed Omega (S) matrix from the full model, or NULL.
#'   When non-NULL (GMM2S/CUE paths), used directly instead of recomputing
#'   from residuals. This ensures the C-statistic uses the same first-step
#'   Omega as the main model's J statistic (Hayashi 2000, p. 220).
#' @param j_full The full model's own overidentification statistic (Stata's `e(j)`/`e(sargan)`), or NULL. Stata computes the C-statistic as `cstat = j - cj` (ivreg2.ado:1547) where `j` is the main model's reported statistic, not a re-minimization of the fixed-S quadratic; the two coincide for 2SLS/GMM2S but differ for CUE (optimum of the self-consistent objective) and LIML (Sargan at the LIML estimate), since the recursive orthog call forwards `gmm2s` but not `cue`/`liml` (ado:1521-1538). When non-NULL and finite, used as J_full directly.
#' @return Named list with `stat`, `p`, `df`, `test_name`, `tested_vars`,
#'   or NULL if this is not an IV model or orthog_vars is NULL.
#' @keywords internal
.compute_orthog_test <- function(Z, X, y, residuals, rss, weights,
                                  cluster_vec, vcov_type, N, K, L,
                                  orthog_vars, dofminus = 0L,
                                  weight_type = "aweight",
                                  kernel = NULL, bw = NULL,
                                  time_index = NULL,
                                  center = FALSE, psd = NULL,
                                  omega = NULL, j_full = NULL,
                                  sw = FALSE, ivar_vec = NULL) {
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

  # --- Omega_full: use pre-computed omega if provided (GMM2S/CUE paths),
  #     otherwise compute from full model's residuals ---
  if (!is.null(omega)) {
    Omega_full <- omega
  } else if (vcov_type %in% c("iid", "AC")) {
    sigma_sq <- rss / (N - dofminus)
    if (is.null(weights)) {
      ZwZ <- crossprod(Z)
    } else {
      ZwZ <- crossprod(Z, weights * Z)
    }
    if (vcov_type == "AC" && !is.null(kernel)) {
      Omega_full <- .ac_meat(Z, residuals, time_index, kernel, bw,
                              N, dofminus, weights, weight_type, ZwZ)
      Omega_full <- .psd_correct(Omega_full, psd)
    } else {
      Omega_full <- sigma_sq * ZwZ / N
    }
  } else {
    Omega_full <- .compute_omega(Z, residuals, weights, cluster_vec, N,
                                  dofminus = dofminus, weight_type = weight_type,
                                  kernel = kernel, bw = bw,
                                  time_index = time_index,
                                  center = center, psd = psd,
                                  sw = sw, ivar_vec = ivar_vec)
  }

  # Structural rank bound of the cluster meat (also bounds every submatrix).
  # Stata sets cstat to missing when rankS < iv1_ct (ivreg2.ado:1793-1808);
  # the structural gate reproduces that suppression platform-stably.
  cl_rank_bound <- .cluster_rank_bound(cluster_vec, kernel, center, weights)

  # --- J_full: the full model's own J (Stata: cstat = j - cj, ado:1547). Recompute via the fixed-Omega minimization only when the caller did not supply the model's reported statistic (see @param j_full). ---
  J_full <- if (!is.null(j_full) && is.finite(j_full)) {
    j_full
  } else {
    .compute_j_with_omega(Z, X, y, Omega_full, weights, N,
                          rank_bound = cl_rank_bound)
  }
  if (is.na(J_full)) {
    return(list(stat = NA_real_, p = NA_real_, df = as.integer(q),
                test_name = "C (orthog)",
                tested_vars = orthog_vars))
  }

  # --- Extract Omega submatrix for restricted model ---
  Omega_sub <- Omega_full[keep_idx, keep_idx, drop = FALSE]

  # --- J_r: J statistic of restricted model (using full model's S) ---
  J_r <- .compute_j_with_omega(Z_r, X, y, Omega_sub, weights, N,
                               rank_bound = cl_rank_bound)
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
