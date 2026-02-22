# --------------------------------------------------------------------------
# Anderson-Rubin test (Ticket E3)
# --------------------------------------------------------------------------
# Weak-instrument-robust test of H0: all endogenous coefficients = 0.
# Valid regardless of instrument strength.
#
# Ported from Stata's ivreg2.ado (lines ~6141-6184).
# Uses the y-equation-only sandwich simplification: the full multi-equation
# Kronecker VCV reduces to a standard sandwich on y-equation residuals.


# --------------------------------------------------------------------------
# .compute_anderson_rubin
# --------------------------------------------------------------------------
#' Anderson-Rubin test for weak-instrument-robust inference
#'
#' Computes the Anderson-Rubin F and chi-squared statistics for testing
#' H0: all endogenous variable coefficients are zero in the structural
#' equation. The test is valid regardless of instrument strength.
#'
#' @param Z N x L instrument matrix.
#' @param X N x K regressor matrix.
#' @param y N-vector of responses.
#' @param weights Normalized weights (sum to N), or NULL.
#' @param cluster_vec Cluster membership vector, or NULL.
#' @param vcov_type Character: `"iid"`, `"HC0"`, `"HC1"`, or `"CL"`.
#' @param N,K,L,K1,L1 Integer dimensions.
#' @param M Number of clusters (or NULL).
#' @param endo_names Character vector of endogenous variable names.
#' @param excluded_names Character vector of excluded instrument names.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @return Named list: `f_stat`, `f_p`, `f_df1`, `f_df2`, `chi2_stat`,
#'   `chi2_p`, `chi2_df`.
#' @keywords internal
.compute_anderson_rubin <- function(Z, X, y, weights, cluster_vec,
                                     vcov_type, N, K, L, K1, L1, M,
                                     endo_names, excluded_names,
                                     dofminus = 0L, sdofminus = 0L,
                                     weight_type = "aweight",
                                     kernel = NULL, bw = NULL,
                                     time_index = NULL) {

  # --- A. Index vectors ---
  endo_idx <- match(endo_names, colnames(X))
  excl_idx <- match(excluded_names, colnames(Z))

  # --- B. Reduced-form fit: regress y on Z ---
  if (is.null(weights)) {
    rf_fit <- stats::lm.fit(Z, y)
  } else {
    rf_fit <- stats::lm.wfit(Z, y, weights)
  }
  rf_coefs <- rf_fit$coefficients   # L-vector
  rf_resid <- rf_fit$residuals      # N-vector

  # --- C. Bread: (Z'WZ)^{-1} from QR ---
  p <- rf_fit$rank
  R_qr <- rf_fit$qr$qr[seq_len(p), seq_len(p), drop = FALSE]
  ZtWZ_inv <- chol2inv(R_qr)
  piv <- rf_fit$qr$pivot
  piv_order <- order(piv[seq_len(p)])
  ZtWZ_inv <- ZtWZ_inv[piv_order, piv_order, drop = FALSE]

  # --- D. RVR computation ---
  Rb <- rf_coefs[excl_idx]

  if (vcov_type == "iid") {
    # Classical VCV: sigma2_y * (Z'WZ)^{-1}
    rss_y <- if (is.null(weights)) {
      drop(crossprod(rf_resid))
    } else {
      sum(weights * rf_resid^2)
    }
    sigma2_y <- rss_y / (N - dofminus)
    RVR <- sigma2_y * ZtWZ_inv[excl_idx, excl_idx, drop = FALSE]
  } else {
    # Robust sandwich (HC0/HC1/CL/HAC/cluster+kernel)
    if (!is.null(cluster_vec) && !is.null(kernel)) {
      # Cluster + kernel (DK or Thompson)
      if (is.list(cluster_vec)) {
        scores <- .cl_scores(Z, rf_resid, weights)
        shat1 <- crossprod(rowsum(scores, cluster_vec[[1L]], reorder = FALSE))
        shat1 <- (shat1 + t(shat1)) / 2
        shat2 <- .cluster_kernel_meat(Z, rf_resid, time_index, kernel, bw,
                                       weights, weight_type)
        shat3 <- .hac_scores_meat(scores, time_index, kernel, bw)
        meat <- shat1 + shat2 - shat3
      } else {
        meat <- .cluster_kernel_meat(Z, rf_resid, time_index, kernel, bw,
                                     weights, weight_type)
      }
    } else if (!is.null(cluster_vec)) {
      scores <- .cl_scores(Z, rf_resid, weights)
      meat <- .cluster_meat(scores, cluster_vec)
    } else if (!is.null(kernel)) {
      meat <- .hac_meat(Z, rf_resid, time_index, kernel, bw,
                        weights, weight_type)
    } else {
      meat <- .hc_meat(Z, rf_resid, weights, weight_type)
    }
    sandwich_full <- ZtWZ_inv %*% meat %*% ZtWZ_inv

    # Extract excluded-IV block and force symmetry
    RVR <- sandwich_full[excl_idx, excl_idx, drop = FALSE]
    RVR <- (RVR + t(RVR)) / 2
  }

  # --- E. Wald statistic ---
  wald <- tryCatch({
    # Guard against numerically singular RVR (e.g., very few clusters make
    # the cluster meat rank-deficient; chol may succeed but give garbage).
    if (ncol(RVR) > 1L && rcond(RVR) < sqrt(.Machine$double.eps))
      stop("nearly singular")
    R_chol <- chol(RVR)
    z <- forwardsolve(t(R_chol), Rb)
    drop(crossprod(z))
  }, error = function(e) {
    warning("Anderson-Rubin: RVR matrix is singular; returning NA.",
            call. = FALSE)
    NA_real_
  })

  # --- E2. Omega normalization adjustment ---
  # Stata's HC omega divides by (N-dofminus) while our raw sandwich is unscaled.
  # For IID, sigma2_y already absorbs the dofminus factor.
  # For cluster, Stata's omega divides by N (no dofminus), so no adjustment.
  # For HC (non-cluster): scale Wald by (N-dofminus)/N to match Stata.
  if (!is.na(wald) && vcov_type != "iid" && is.null(cluster_vec)) {
    wald <- wald * (N - dofminus) / N
  }

  # --- F. F conversion ---
  # Non-cluster: F = Wald / (N-dofminus) * (N-L-dofminus-sdofminus) / L1
  # Cluster: F = Wald / (N-1) * (N-L-sdofminus) * (M-1)/M / L1
  if (is.na(wald)) {
    f_stat <- NA_real_
    f_p    <- NA_real_
  } else if (!is.null(cluster_vec)) {
    f_stat <- wald / (N - 1) * (N - L - sdofminus) * (M - 1) / M / L1
    f_p    <- stats::pf(f_stat, df1 = L1, df2 = M - 1L, lower.tail = FALSE)
  } else {
    f_stat <- wald / (N - dofminus) * (N - L - dofminus - sdofminus) / L1
    f_p    <- stats::pf(f_stat, df1 = L1,
                         df2 = N - L - dofminus - sdofminus,
                         lower.tail = FALSE)
  }

  f_df1 <- as.integer(L1)
  f_df2 <- if (!is.null(cluster_vec)) {
    as.integer(M - 1L)
  } else {
    as.integer(N - L - dofminus - sdofminus)
  }

  # --- G. Chi-sq ---
  chi2_stat <- if (is.na(wald)) NA_real_ else wald
  chi2_p    <- if (is.na(wald)) {
    NA_real_
  } else {
    stats::pchisq(wald, df = L1, lower.tail = FALSE)
  }
  chi2_df <- as.integer(L1)

  list(
    f_stat    = f_stat,
    f_p       = f_p,
    f_df1     = f_df1,
    f_df2     = f_df2,
    chi2_stat = chi2_stat,
    chi2_p    = chi2_p,
    chi2_df   = chi2_df
  )
}
