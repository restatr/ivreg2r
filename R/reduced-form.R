# --------------------------------------------------------------------------
# Reduced-form regression (Ticket J3)
# --------------------------------------------------------------------------
# Computes reduced-form (y ~ Z) and optionally system (y + all X_j ~ Z)
# regressions with Stata-parity small-sample corrections.
#
# Ported from Stata's ivreg2.ado PostFirstRF (lines ~3033-3244).
# The VCV scaling always applies OLS-style small-sample df corrections
# (lines 3091-3097), regardless of the main model's `small` option.


# --------------------------------------------------------------------------
# .compute_reduced_form
# --------------------------------------------------------------------------
#' Compute reduced-form regression(s) for IV models
#'
#' In `"rf"` mode, fits y ~ Z and returns coefficients, VCV, residuals,
#' fitted values, F-test of excluded instruments, and RMSE.
#' In `"system"` mode, additionally fits each endogenous X_j ~ Z and
#' constructs the cross-equation VCV.
#'
#' @param mode Character: `"rf"` or `"system"`.
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
#' @param depvar_name Character: name of the dependent variable.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @return A list with `mode` and per-equation (rf) or multi-equation (system)
#'   regression results.
#' @keywords internal
.compute_reduced_form <- function(mode, Z, X, y, weights, cluster_vec,
                                   vcov_type, N, K, L, K1, L1, M,
                                   endo_names, excluded_names,
                                   depvar_name,
                                   dofminus = 0L, sdofminus = 0L) {

  # --- A. Common setup ---
  excl_idx <- match(excluded_names, colnames(Z))
  endo_idx <- match(endo_names, colnames(X))
  sqrt_w <- if (!is.null(weights)) sqrt(weights) else NULL

  # Number of equations: 1 for rf, 1+K1 for system
  n_eq <- if (mode == "rf") 1L else 1L + K1

  # --- B. Fit all RF/first-stage regressions ---
  # Combine y and endogenous columns into a single response matrix
  if (mode == "rf") {
    Y <- matrix(y, ncol = 1L)
    colnames(Y) <- depvar_name
  } else {
    X1 <- X[, endo_idx, drop = FALSE]
    Y <- cbind(y, X1)
    colnames(Y) <- c(depvar_name, endo_names)
  }

  if (is.null(weights)) {
    rf_fit <- stats::lm.fit(Z, Y)
  } else {
    rf_fit <- stats::lm.wfit(Z, Y, weights)
  }

  rf_coefs  <- as.matrix(rf_fit$coefficients)    # L x n_eq
  rf_resid  <- as.matrix(rf_fit$residuals)        # N x n_eq
  rf_fitted <- as.matrix(rf_fit$fitted.values)    # N x n_eq

  # --- C. Bread: (Z'WZ)^{-1} from QR ---
  p <- rf_fit$rank
  R_qr <- rf_fit$qr$qr[seq_len(p), seq_len(p), drop = FALSE]
  ZtWZ_inv <- chol2inv(R_qr)
  piv <- rf_fit$qr$pivot
  piv_order <- order(piv[seq_len(p)])
  ZtWZ_inv <- ZtWZ_inv[piv_order, piv_order, drop = FALSE]

  # --- D. VCV computation ---
  # PostFirstRF always applies OLS-style small-sample correction
  # (lines 3091-3097), regardless of the main model's `small` option:
  #   Non-cluster: V * (N - dofminus) / (N - L - dofminus - sdofminus)
  #   Cluster: V * (N - 1) / (N - L - sdofminus) * M / (M - 1)
  ss_df <- N - L - dofminus - sdofminus  # small-sample residual df

  if (mode == "rf") {
    # --- Single-equation VCV ---
    vcov_rf <- .rf_equation_vcov(
      Z = Z, resid = rf_resid[, 1L], coefs = rf_coefs[, 1L],
      ZtWZ_inv = ZtWZ_inv, weights = weights, sqrt_w = sqrt_w,
      cluster_vec = cluster_vec, vcov_type = vcov_type,
      N = N, L = L, M = M, dofminus = dofminus, sdofminus = sdofminus,
      excl_idx = excl_idx, L1 = L1
    )
    colnames(vcov_rf$vcov) <- rownames(vcov_rf$vcov) <- colnames(Z)

    # RMSE (always small-sample denominator)
    rss_y <- if (is.null(weights)) {
      drop(crossprod(rf_resid[, 1L]))
    } else {
      sum(weights * rf_resid[, 1L]^2)
    }
    sigma_y <- sqrt(rss_y / ss_df)

    # Partial R² of excluded instruments
    partial_r2 <- .rf_partial_r2(vcov_rf$wald_classical, N, dofminus)

    result <- list(
      mode          = "rf",
      depvar        = depvar_name,
      coefficients  = rf_coefs[, 1L],
      vcov          = vcov_rf$vcov,
      residuals     = rf_resid[, 1L],
      fitted.values = rf_fitted[, 1L],
      sigma         = sigma_y,
      f_stat        = vcov_rf$f_stat,
      f_p           = vcov_rf$f_p,
      f_df1         = vcov_rf$f_df1,
      f_df2         = vcov_rf$f_df2,
      chi2_stat     = vcov_rf$chi2_stat,
      chi2_p        = vcov_rf$chi2_p,
      chi2_df       = vcov_rf$chi2_df,
      partial_r2    = partial_r2
    )

  } else {
    # --- System (multi-equation) VCV ---
    eq_names <- colnames(Y)

    if (vcov_type == "iid") {
      # IID: Kronecker product V = sigma_matrix %x% (Z'WZ)^{-1}
      # Sigma matrix: E[e_i * e_j'] estimated as sum(w * e_i * e_j) / (N - dofminus)
      sigma_mat <- matrix(0, n_eq, n_eq)
      for (i in seq_len(n_eq)) {
        for (j in seq(i, n_eq)) {
          if (is.null(weights)) {
            sigma_mat[i, j] <- drop(crossprod(rf_resid[, i], rf_resid[, j]))
          } else {
            sigma_mat[i, j] <- sum(weights * rf_resid[, i] * rf_resid[, j])
          }
          sigma_mat[j, i] <- sigma_mat[i, j]
        }
      }
      sigma_mat <- sigma_mat / (N - dofminus)
      system_vcov <- kronecker(sigma_mat, ZtWZ_inv)

    } else {
      # Robust/cluster: stack scores across equations
      # scores = [Z*e_1, Z*e_2, ..., Z*e_{n_eq}]  →  N x (n_eq * L)
      scores_list <- vector("list", n_eq)
      for (j in seq_len(n_eq)) {
        if (is.null(weights)) {
          scores_list[[j]] <- Z * rf_resid[, j]
        } else {
          scores_list[[j]] <- (sqrt_w * Z) * (sqrt_w * rf_resid[, j])
        }
      }
      scores_stacked <- do.call(cbind, scores_list)  # N x (n_eq * L)

      if (!is.null(cluster_vec)) {
        meat <- .cluster_meat(scores_stacked, cluster_vec)
      } else {
        meat <- crossprod(scores_stacked)
      }

      # Bread: I_{n_eq} %x% (Z'WZ)^{-1}
      bread_kron <- kronecker(diag(n_eq), ZtWZ_inv)
      system_vcov <- bread_kron %*% meat %*% bread_kron
    }

    # Apply Stata small-sample correction (PostFirstRF lines 3091-3097)
    if (!is.null(cluster_vec)) {
      system_vcov <- system_vcov * (N - 1) / (N - L - sdofminus) *
        M / (M - 1)
    } else {
      system_vcov <- system_vcov * (N - dofminus) /
        (N - L - dofminus - sdofminus)
    }

    system_vcov <- (system_vcov + t(system_vcov)) / 2

    # Label system VCV
    sys_names <- as.vector(outer(colnames(Z), eq_names,
                                  function(z, e) paste0(e, ":", z)))
    colnames(system_vcov) <- rownames(system_vcov) <- sys_names

    # Per-equation results
    sigma_vec <- numeric(n_eq)
    equations <- vector("list", n_eq)
    names(equations) <- eq_names

    for (j in seq_len(n_eq)) {
      rss_j <- if (is.null(weights)) {
        drop(crossprod(rf_resid[, j]))
      } else {
        sum(weights * rf_resid[, j]^2)
      }
      sigma_vec[j] <- sqrt(rss_j / ss_df)

      # Per-equation VCV block
      idx <- ((j - 1L) * L + 1L):(j * L)
      vcov_j <- system_vcov[idx, idx, drop = FALSE]

      # F-test and chi-sq of excluded instruments from equation VCV
      eq_stats <- .rf_excl_test_from_vcov(
        coefs = rf_coefs[, j], vcov = vcov_j, excl_idx = excl_idx,
        L1 = L1, N = N, L = L, M = M, ss_df = ss_df,
        cluster_vec = cluster_vec
      )

      # Partial R² (always classical, VCE-independent)
      wald_classical_j <- .rf_classical_wald(
        Z = Z, resid = rf_resid[, j], coefs = rf_coefs[, j],
        ZtWZ_inv = ZtWZ_inv, weights = weights, N = N,
        dofminus = dofminus, excl_idx = excl_idx
      )

      equations[[j]] <- list(
        f_stat     = eq_stats$f_stat,
        f_p        = eq_stats$f_p,
        f_df1      = eq_stats$f_df1,
        f_df2      = eq_stats$f_df2,
        chi2_stat  = eq_stats$chi2_stat,
        chi2_p     = eq_stats$chi2_p,
        chi2_df    = eq_stats$chi2_df,
        partial_r2 = .rf_partial_r2(wald_classical_j, N, dofminus)
      )
    }
    names(sigma_vec) <- eq_names

    result <- list(
      mode          = "system",
      depvar        = eq_names,
      coefficients  = rf_coefs,
      vcov          = system_vcov,
      residuals     = rf_resid,
      fitted.values = rf_fitted,
      sigma         = sigma_vec,
      equations     = equations
    )
  }

  result
}


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

#' Compute single-equation VCV with Stata PostFirstRF scaling
#' @keywords internal
#' @noRd
.rf_equation_vcov <- function(Z, resid, coefs, ZtWZ_inv, weights, sqrt_w,
                                cluster_vec, vcov_type,
                                N, L, M, dofminus, sdofminus,
                                excl_idx, L1) {
  if (vcov_type == "iid") {
    # Classical: sigma2 * (Z'WZ)^{-1}
    rss <- if (is.null(weights)) {
      drop(crossprod(resid))
    } else {
      sum(weights * resid^2)
    }
    sigma2 <- rss / (N - dofminus)
    vcov_raw <- sigma2 * ZtWZ_inv
  } else {
    # Robust sandwich
    if (is.null(weights)) {
      scores <- Z * resid
    } else {
      scores <- (sqrt_w * Z) * (sqrt_w * resid)
    }
    if (!is.null(cluster_vec)) {
      meat <- .cluster_meat(scores, cluster_vec)
    } else {
      meat <- crossprod(scores)
    }
    vcov_raw <- ZtWZ_inv %*% meat %*% ZtWZ_inv
  }

  # Apply PostFirstRF small-sample correction (always applied)
  if (!is.null(cluster_vec)) {
    vcov_raw <- vcov_raw * (N - 1) / (N - L - sdofminus) *
      M / (M - 1)
  } else {
    vcov_raw <- vcov_raw * (N - dofminus) / (N - L - dofminus - sdofminus)
  }

  vcov_raw <- (vcov_raw + t(vcov_raw)) / 2

  # F-test of excluded instruments from this VCV
  Rb <- coefs[excl_idx]
  RVR <- vcov_raw[excl_idx, excl_idx, drop = FALSE]

  wald <- tryCatch({
    R_chol <- chol(RVR)
    z <- forwardsolve(t(R_chol), Rb)
    drop(crossprod(z))
  }, error = function(e) NA_real_)

  # F conversion: Wald / L1
  if (is.na(wald)) {
    f_stat <- NA_real_
    f_p    <- NA_real_
  } else {
    f_stat <- wald / L1
    f_df2 <- if (!is.null(cluster_vec)) {
      as.integer(M - 1L)
    } else {
      as.integer(N - L - dofminus - sdofminus)
    }
    f_p <- stats::pf(f_stat, df1 = L1, df2 = f_df2, lower.tail = FALSE)
  }

  f_df1 <- as.integer(L1)
  f_df2 <- if (!is.null(cluster_vec)) {
    as.integer(M - 1L)
  } else {
    as.integer(N - L - dofminus - sdofminus)
  }

  # Chi-sq
  chi2_stat <- if (is.na(wald)) NA_real_ else wald
  chi2_p <- if (is.na(wald)) NA_real_ else {
    stats::pchisq(wald, df = L1, lower.tail = FALSE)
  }

  # Classical Wald for partial R² (VCE-independent)
  wald_classical <- .rf_classical_wald(
    Z = Z, resid = resid, coefs = coefs, ZtWZ_inv = ZtWZ_inv,
    weights = weights, N = N, dofminus = dofminus, excl_idx = excl_idx
  )

  list(
    vcov           = vcov_raw,
    f_stat         = f_stat,
    f_p            = f_p,
    f_df1          = f_df1,
    f_df2          = f_df2,
    chi2_stat      = chi2_stat,
    chi2_p         = chi2_p,
    chi2_df        = as.integer(L1),
    wald_classical = wald_classical
  )
}


#' Classical Wald statistic (for partial R², VCE-independent)
#' @keywords internal
#' @noRd
.rf_classical_wald <- function(Z, resid, coefs, ZtWZ_inv, weights,
                                 N, dofminus, excl_idx) {
  rss <- if (is.null(weights)) {
    drop(crossprod(resid))
  } else {
    sum(weights * resid^2)
  }
  sigma2 <- rss / (N - dofminus)
  Rb <- coefs[excl_idx]
  RVR <- sigma2 * ZtWZ_inv[excl_idx, excl_idx, drop = FALSE]
  tryCatch({
    R_chol <- chol(RVR)
    z <- forwardsolve(t(R_chol), Rb)
    drop(crossprod(z))
  }, error = function(e) NA_real_)
}


#' Compute partial R² from classical Wald
#' @keywords internal
#' @noRd
.rf_partial_r2 <- function(wald_classical, N, dofminus) {
  if (is.na(wald_classical)) return(NA_real_)
  wald_scaled <- wald_classical / (N - dofminus)
  wald_scaled / (1 + wald_scaled)
}


#' F-test of excluded IVs from a VCV block (system mode helper)
#' @keywords internal
#' @noRd
.rf_excl_test_from_vcov <- function(coefs, vcov, excl_idx,
                                      L1, N, L, M, ss_df,
                                      cluster_vec) {
  Rb <- coefs[excl_idx]
  RVR <- vcov[excl_idx, excl_idx, drop = FALSE]

  wald <- tryCatch({
    R_chol <- chol(RVR)
    z <- forwardsolve(t(R_chol), Rb)
    drop(crossprod(z))
  }, error = function(e) NA_real_)

  # F = Wald / L1
  if (is.na(wald)) {
    f_stat <- NA_real_
    f_p    <- NA_real_
  } else {
    f_stat <- wald / L1
    f_df2 <- if (!is.null(cluster_vec)) {
      as.integer(M - 1L)
    } else {
      as.integer(ss_df)
    }
    f_p <- stats::pf(f_stat, df1 = L1, df2 = f_df2, lower.tail = FALSE)
  }

  f_df1 <- as.integer(L1)
  f_df2 <- if (!is.null(cluster_vec)) {
    as.integer(M - 1L)
  } else {
    as.integer(ss_df)
  }

  chi2_stat <- if (is.na(wald)) NA_real_ else wald
  chi2_p <- if (is.na(wald)) NA_real_ else {
    stats::pchisq(wald, df = L1, lower.tail = FALSE)
  }

  list(
    f_stat    = f_stat,
    f_p       = f_p,
    f_df1     = f_df1,
    f_df2     = f_df2,
    chi2_stat = chi2_stat,
    chi2_p    = chi2_p,
    chi2_df   = as.integer(L1)
  )
}
