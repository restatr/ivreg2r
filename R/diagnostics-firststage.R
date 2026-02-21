# --------------------------------------------------------------------------
# First-stage diagnostics (Ticket E1)
# --------------------------------------------------------------------------
# For each endogenous variable: first-stage F-test of excluded instruments,
# partial R², Shea partial R², and RMSE.
#
# Ported from Stata's ivreg2.ado first-stage block (lines ~6130-6350).


# --------------------------------------------------------------------------
# .compute_first_stage
# --------------------------------------------------------------------------
#' First-stage regression diagnostics
#'
#' For each endogenous variable, fits the first-stage OLS regression on the
#' full instrument set Z, then computes: F-test of excluded instruments,
#' partial R², Shea partial R², and RMSE.
#'
#' @param X N x K regressor matrix.
#' @param Z N x L instrument matrix.
#' @param weights Normalized weights (sum to N), or NULL.
#' @param cluster_vec Cluster membership vector, or NULL.
#' @param vcov_type Character: `"iid"`, `"HC0"`, `"HC1"`, or `"CL"`.
#' @param endo_names Character vector of endogenous variable names.
#' @param excluded_names Character vector of excluded instrument names.
#' @param N,K,L,K1,L1 Integer dimensions.
#' @param M Number of clusters (or NULL).
#' @param bread_2sls K x K bread from 2SLS: \eqn{(X_{hat}'WX_{hat})^{-1}}.
#' @details
#' First-stage F-statistics always use small-sample degrees of freedom
#' (F(L1, N-L) or F(L1, M-1) for clustered models), regardless of the main
#' model's \code{small} option. This matches Stata's \code{ivreg2} behavior.
#'
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @return Named list keyed by endogenous variable names, each element
#'   containing: `f_stat`, `f_p`, `f_df1`, `f_df2`, `partial_r2`,
#'   `shea_partial_r2`, `rmse`, `coefficients`, `residuals`, `fitted.values`,
#'   plus SW/AP diagnostics: `sw_f`, `sw_f_p`, `sw_f_df1`, `sw_f_df2`,
#'   `sw_chi2`, `sw_chi2_p`, `sw_partial_r2`, `ap_f`, `ap_f_p`, `ap_f_df1`,
#'   `ap_f_df2`, `ap_chi2`, `ap_chi2_p`, `ap_partial_r2`.
#' @keywords internal
.compute_first_stage <- function(X, Z, weights, cluster_vec,
                                  vcov_type, endo_names, excluded_names,
                                  N, K, L, K1, L1, M,
                                  bread_2sls,
                                  dofminus = 0L, sdofminus = 0L,
                                  weight_type = "aweight",
                                  kernel = NULL, bw = NULL,
                                  time_index = NULL) {

  # --- A4: Index vectors ---
  endo_idx <- match(endo_names, colnames(X))
  excl_idx <- match(excluded_names, colnames(Z))

  # --- A1: Fit all first-stage regressions at once ---
  X1 <- X[, endo_idx, drop = FALSE]
  if (is.null(weights)) {
    fs_fit <- stats::lm.fit(Z, X1)
  } else {
    fs_fit <- stats::lm.wfit(Z, X1, weights)
  }
  fs_coefs <- as.matrix(fs_fit$coefficients)      # L x K1
  fs_resid <- as.matrix(fs_fit$residuals)          # N x K1
  fs_fitted <- as.matrix(fs_fit$fitted.values)     # N x K1

  # --- A2: Bread (Z'WZ)^{-1} from shared QR ---
  p <- fs_fit$rank
  R_qr <- fs_fit$qr$qr[seq_len(p), seq_len(p), drop = FALSE]
  ZtWZ_inv <- chol2inv(R_qr)
  piv <- fs_fit$qr$pivot
  piv_order <- order(piv[seq_len(p)])
  ZtWZ_inv <- ZtWZ_inv[piv_order, piv_order, drop = FALSE]

  # --- A3: Precompute sqrt(weights); (X'WX)^{-1} for Shea partial R² ---
  sqrt_w <- if (!is.null(weights)) sqrt(weights) else NULL

  XXinv <- tryCatch({
    if (is.null(sqrt_w)) {
      chol2inv(chol(crossprod(X)))
    } else {
      chol2inv(chol(crossprod(sqrt_w * X)))
    }
  }, error = function(e) NULL)

  # --- B: Per-endogenous-variable loop ---
  results <- vector("list", K1)
  names(results) <- endo_names
  wald_classical_vec <- numeric(K1)
  wald_vec <- numeric(K1)

  for (j in seq_len(K1)) {
    coef_j <- fs_coefs[, j]
    resid_j <- fs_resid[, j]
    fitted_j <- fs_fitted[, j]

    # B1: RSS
    rss_j <- if (is.null(weights)) {
      drop(crossprod(resid_j))
    } else {
      sum(weights * resid_j^2)
    }

    # B2: RMSE (always small-sample denominator)
    rmse_j <- sqrt(rss_j / (N - L - dofminus - sdofminus))

    # B3: Classical Wald (for partial R²; VCE-independent)
    sigma2_j <- rss_j / (N - dofminus)
    Rb_j <- coef_j[excl_idx]
    RVR_classical <- sigma2_j * ZtWZ_inv[excl_idx, excl_idx, drop = FALSE]

    wald_classical_j <- tryCatch({
      R_chol <- chol(RVR_classical)
      z <- forwardsolve(t(R_chol), Rb_j)
      drop(crossprod(z))
    }, error = function(e) NA_real_)

    # B4: Partial R² (always from classical Wald)
    # Stata: pr2 = (Wald/(N-dofminus)) / (1 + (Wald/(N-dofminus)))
    partial_r2_j <- if (is.na(wald_classical_j)) {
      NA_real_
    } else {
      wald_scaled <- wald_classical_j / (N - dofminus)
      wald_scaled / (1 + wald_scaled)
    }

    # B5: Robust Wald (only when VCE is not IID/AC)
    if (vcov_type %in% c("iid", "AC")) {
      wald_j <- wald_classical_j
    } else {
      # Raw robust sandwich (no finite-sample corrections)
      if (!is.null(cluster_vec)) {
        scores <- .cl_scores(Z, resid_j, weights)
        meat <- .cluster_meat(scores, cluster_vec)
      } else if (!is.null(kernel)) {
        meat <- .hac_meat(Z, resid_j, time_index, kernel, bw,
                          weights, weight_type)
      } else {
        meat <- .hc_meat(Z, resid_j, weights, weight_type)
      }
      robust_vcov <- ZtWZ_inv %*% meat %*% ZtWZ_inv
      RVR_robust <- robust_vcov[excl_idx, excl_idx, drop = FALSE]

      wald_j <- tryCatch({
        if (ncol(RVR_robust) > 1L && rcond(RVR_robust) < sqrt(.Machine$double.eps))
          stop("nearly singular")
        R_chol <- chol(RVR_robust)
        z <- forwardsolve(t(R_chol), Rb_j)
        drop(crossprod(z))
      }, error = function(e) NA_real_)

      # Omega normalization: HC omega divides by (N-dofminus) in Stata,
      # our raw sandwich is unscaled. Cluster omega uses N (no adjustment).
      if (!is.na(wald_j) && is.null(cluster_vec)) {
        wald_j <- wald_j * (N - dofminus) / N
      }
    }

    # B6: F conversion
    # Non-cluster: F = Wald / (N-dofminus) * (N-L-dofminus-sdofminus) / L1
    # Cluster: F = Wald / (N-1) * (N-L-sdofminus) * (M-1)/M / L1
    if (is.na(wald_j)) {
      f_stat_j <- NA_real_
      f_p_j <- NA_real_
    } else if (!is.null(cluster_vec)) {
      f_stat_j <- wald_j / (N - 1) * (N - L - sdofminus) *
        (M - 1) / M / L1
      f_p_j <- stats::pf(f_stat_j, df1 = L1, df2 = M - 1L,
                          lower.tail = FALSE)
    } else {
      f_stat_j <- wald_j / (N - dofminus) *
        (N - L - dofminus - sdofminus) / L1
      f_p_j <- stats::pf(f_stat_j, df1 = L1,
                          df2 = N - L - dofminus - sdofminus,
                          lower.tail = FALSE)
    }

    f_df1_j <- as.integer(L1)
    f_df2_j <- if (!is.null(cluster_vec)) {
      as.integer(M - 1L)
    } else {
      as.integer(N - L - dofminus - sdofminus)
    }

    wald_classical_vec[j] <- wald_classical_j
    wald_vec[j] <- wald_j

    results[[j]] <- list(
      f_stat          = f_stat_j,
      f_p             = f_p_j,
      f_df1           = f_df1_j,
      f_df2           = f_df2_j,
      partial_r2      = partial_r2_j,
      shea_partial_r2 = NA_real_,   # filled in Phase C
      rmse            = rmse_j,
      coefficients    = coef_j,
      residuals       = resid_j,
      fitted.values   = fitted_j
    )
  }

  # --- C: Shea partial R² (after loop) ---
  if (!is.null(XXinv)) {
    diag_XX <- diag(XXinv)
    diag_bread <- diag(bread_2sls)
    for (j in seq_len(K1)) {
      idx <- endo_idx[j]
      denom <- diag_bread[idx]
      if (abs(denom) >= .Machine$double.eps) {
        results[[j]]$shea_partial_r2 <- unname(diag_XX[idx] / denom)
      }
    }
  }

  # --- D: Sanderson-Windmeijer / Angrist-Pischke diagnostics ---
  Fdf1_sw <- as.integer(L1 - K1 + 1L)
  df2_sw <- if (!is.null(cluster_vec)) {
    as.integer(M - 1L)
  } else {
    as.integer(N - L - dofminus - sdofminus)
  }

  if (K1 == 1L) {
    # K1=1 shortcut: SW/AP collapse to standard first-stage values
    for (j in seq_len(K1)) {
      res <- results[[j]]
      wald_j <- wald_vec[j]
      chi2_p_j <- if (is.na(wald_j)) {
        NA_real_
      } else {
        stats::pchisq(wald_j, df = Fdf1_sw, lower.tail = FALSE)
      }
      results[[j]]$sw_f          <- res$f_stat
      results[[j]]$sw_f_p        <- res$f_p
      results[[j]]$sw_f_df1      <- Fdf1_sw
      results[[j]]$sw_f_df2      <- res$f_df2
      results[[j]]$sw_chi2       <- wald_j
      results[[j]]$sw_chi2_p     <- chi2_p_j
      results[[j]]$sw_partial_r2 <- res$partial_r2
      results[[j]]$ap_f          <- res$f_stat
      results[[j]]$ap_f_p        <- res$f_p
      results[[j]]$ap_f_df1      <- Fdf1_sw
      results[[j]]$ap_f_df2      <- res$f_df2
      results[[j]]$ap_chi2       <- wald_j
      results[[j]]$ap_chi2_p     <- chi2_p_j
      results[[j]]$ap_partial_r2 <- res$partial_r2
    }
  } else {
    # K1 > 1: full SW/AP computation
    Xhat <- X
    Xhat[, endo_idx] <- fs_fitted

    for (j in seq_len(K1)) {
      X1_j <- X[, endo_idx[j]]
      Xhat_mj <- Xhat[, -endo_idx[j], drop = FALSE]
      X_mj <- X[, -endo_idx[j], drop = FALSE]

      # D1: Projection coefficients (OLS of X1_j on Xhat_mj)
      b1 <- tryCatch({
        if (is.null(sqrt_w)) {
          solve(crossprod(Xhat_mj), crossprod(Xhat_mj, X1_j))
        } else {
          W_Xhat_mj <- sqrt_w * Xhat_mj
          W_X1_j <- sqrt_w * X1_j
          solve(crossprod(W_Xhat_mj), crossprod(W_Xhat_mj, W_X1_j))
        }
      }, error = function(e) NULL)

      if (is.null(b1)) {
        # Singular projection — all SW/AP are NA
        for (prefix in c("sw", "ap")) {
          results[[j]][[paste0(prefix, "_f")]]          <- NA_real_
          results[[j]][[paste0(prefix, "_f_p")]]        <- NA_real_
          results[[j]][[paste0(prefix, "_f_df1")]]      <- Fdf1_sw
          results[[j]][[paste0(prefix, "_f_df2")]]      <- df2_sw
          results[[j]][[paste0(prefix, "_chi2")]]       <- NA_real_
          results[[j]][[paste0(prefix, "_chi2_p")]]     <- NA_real_
          results[[j]][[paste0(prefix, "_partial_r2")]] <- NA_real_
        }
        next
      }

      # D2: AP and SW residuals
      e_ap <- drop(X1_j - Xhat_mj %*% b1)
      e_sw <- drop(X1_j - X_mj %*% b1)

      # D3: Compute stats for both types
      for (type in c("sw", "ap")) {
        e_j <- if (type == "ap") e_ap else e_sw

        # Regress e on Z
        if (is.null(sqrt_w)) {
          b2 <- drop(ZtWZ_inv %*% crossprod(Z, e_j))
        } else {
          b2 <- drop(ZtWZ_inv %*% crossprod(sqrt_w * Z, sqrt_w * e_j))
        }
        resid_aux <- e_j - drop(Z %*% b2)

        # Classical Wald
        rss_aux <- if (is.null(sqrt_w)) {
          drop(crossprod(resid_aux))
        } else {
          sum(weights * resid_aux^2)
        }
        sigma2_aux <- rss_aux / (N - dofminus)
        Rb <- b2[excl_idx]
        RVR_cl <- sigma2_aux * ZtWZ_inv[excl_idx, excl_idx, drop = FALSE]

        wald_cl <- tryCatch({
          R_chol <- chol(RVR_cl)
          z <- forwardsolve(t(R_chol), Rb)
          drop(crossprod(z))
        }, error = function(e) NA_real_)

        # Partial R² (always classical)
        # Stata: pr2 = (Wald/(N-dofminus)) / (1 + (Wald/(N-dofminus)))
        pr2 <- if (is.na(wald_cl)) {
          NA_real_
        } else {
          wald_cl_scaled <- wald_cl / (N - dofminus)
          wald_cl_scaled / (1 + wald_cl_scaled)
        }

        # Robust Wald (if VCE is not IID/AC)
        if (vcov_type %in% c("iid", "AC")) {
          wald_sw <- wald_cl
        } else {
          if (!is.null(cluster_vec)) {
            scores <- .cl_scores(Z, resid_aux, weights)
            meat <- .cluster_meat(scores, cluster_vec)
          } else if (!is.null(kernel)) {
            meat <- .hac_meat(Z, resid_aux, time_index, kernel, bw,
                              weights, weight_type)
          } else {
            meat <- .hc_meat(Z, resid_aux, weights, weight_type)
          }
          robust_vcov <- ZtWZ_inv %*% meat %*% ZtWZ_inv
          robust_vcov <- (robust_vcov + t(robust_vcov)) / 2
          RVR_robust <- robust_vcov[excl_idx, excl_idx, drop = FALSE]

          wald_sw <- tryCatch({
            if (ncol(RVR_robust) > 1L && rcond(RVR_robust) < sqrt(.Machine$double.eps))
              stop("nearly singular")
            R_chol <- chol(RVR_robust)
            z <- forwardsolve(t(R_chol), Rb)
            drop(crossprod(z))
          }, error = function(e) NA_real_)

          # Omega normalization (HC: Stata divides by N-dofminus, not N)
          if (!is.na(wald_sw) && is.null(cluster_vec)) {
            wald_sw <- wald_sw * (N - dofminus) / N
          }
        }

        # F conversion
        if (is.na(wald_sw)) {
          f_val <- NA_real_
          f_p <- NA_real_
        } else if (!is.null(cluster_vec)) {
          f_val <- wald_sw / (N - 1) * (N - L - sdofminus) *
            (M - 1) / M / Fdf1_sw
          f_p <- stats::pf(f_val, df1 = Fdf1_sw, df2 = M - 1L,
                           lower.tail = FALSE)
        } else {
          f_val <- wald_sw / (N - dofminus) *
            (N - L - dofminus - sdofminus) / Fdf1_sw
          f_p <- stats::pf(f_val, df1 = Fdf1_sw, df2 = df2_sw,
                           lower.tail = FALSE)
        }

        # Chi-sq
        chi2 <- if (is.na(wald_sw)) NA_real_ else wald_sw
        chi2_p <- if (is.na(wald_sw)) {
          NA_real_
        } else {
          stats::pchisq(wald_sw, df = Fdf1_sw, lower.tail = FALSE)
        }

        # Store
        results[[j]][[paste0(type, "_f")]]          <- f_val
        results[[j]][[paste0(type, "_f_p")]]        <- f_p
        results[[j]][[paste0(type, "_f_df1")]]      <- Fdf1_sw
        results[[j]][[paste0(type, "_f_df2")]]      <- df2_sw
        results[[j]][[paste0(type, "_chi2")]]       <- chi2
        results[[j]][[paste0(type, "_chi2_p")]]     <- chi2_p
        results[[j]][[paste0(type, "_partial_r2")]] <- pr2
      }
    }
  }

  results
}
