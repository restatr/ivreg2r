# --------------------------------------------------------------------------
# .compute_endogeneity_test
# --------------------------------------------------------------------------
#' Endogeneity test (C-statistic / difference-in-J)
#'
#' Tests H0: specified endogenous regressors are actually exogenous.
#' Computes a difference-of-J (robust/cluster) or difference-of-Sargan (IID)
#' statistic by re-estimating a restricted model where the tested regressors
#' are treated as exogenous (and hence become their own instruments).
#'
#' The same-S-matrix constraint guarantees C >= 0: the S matrix comes from
#' the restricted model (larger instrument set), and the full model's J
#' is re-estimated using the appropriate submatrix of that S.
#'
#' @param Z N x L instrument matrix (full model).
#' @param X N x K regressor matrix.
#' @param y N x 1 response vector.
#' @param residuals N x 1 residual vector from the full model (unused directly,
#'   but included for API consistency).
#' @param rss Residual sum of squares from the full model (unused directly).
#' @param weights Normalized weights (sum to N), or NULL.
#' @param cluster_vec Cluster membership vector, or NULL.
#' @param vcov_type Character: "iid", "HC0", "HC1", or "CL".
#' @param N Number of observations.
#' @param K Number of regressors.
#' @param L Number of instruments.
#' @param K1 Number of endogenous regressors.
#' @param endo_names Character vector of endogenous regressor names.
#' @param endog_vars Character vector of variables to test (subset of
#'   `endo_names`), or NULL to test all.
#' @return Named list with `stat`, `p`, `df`, `test_name`, `tested_vars`,
#'   or NULL if this is not an IV model.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @note The C-statistic for IID models is computed using large-sample
#'   sigma-squared `e'e/(N-dofminus)` for both models. No small-sample
#'   correction is applied even when \code{small = TRUE}. This matches
#'   Stata's \code{ivreg2}.
#' @keywords internal
.compute_endogeneity_test <- function(Z, X, y, residuals, rss, weights,
                                      cluster_vec, vcov_type, N, K, L,
                                      K1, endo_names, endog_vars,
                                      dofminus = 0L) {
  # Default: test all endogenous regressors
  if (is.null(endog_vars)) endog_vars <- endo_names
  q <- length(endog_vars)

  # --- Build restricted model instrument matrix ---
  # Tested endogenous vars become exogenous → add them to instruments
  tested_idx <- match(endog_vars, colnames(X))
  Z_r <- cbind(Z, X[, tested_idx, drop = FALSE])
  L_r <- ncol(Z_r)

  # --- Estimate restricted model via 2SLS with Z_r ---
  # When all endogenous are tested (K1_r = 0), X ⊂ col(Z_r), so X_hat = X
  # and 2SLS reduces to OLS. The same code path handles both cases.
  e_r <- tryCatch({
    if (is.null(weights)) {
      stage1 <- lm.fit(Z_r, X)
    } else {
      stage1 <- lm.wfit(Z_r, X, weights)
    }
    X_hat_r <- as.matrix(stage1$fitted.values)
    if (is.null(weights)) {
      stage2 <- lm.fit(X_hat_r, y)
    } else {
      stage2 <- lm.wfit(X_hat_r, y, weights)
    }
    if (stage2$rank < K) NULL else drop(y - X %*% stage2$coefficients)
  }, error = function(e) {
    NULL
  })

  # Identification failure or rank deficiency in restricted model
  if (is.null(e_r)) {
    return(list(stat = 0, p = NA_real_, df = 0L,
                test_name = "Endogeneity",
                tested_vars = endog_vars))
  }

  rss_r <- if (is.null(weights)) sum(e_r^2) else sum(weights * e_r^2)

  # --- Compute Omega_r (from restricted model) ---
  if (vcov_type == "iid") {
    sigma_r_sq <- rss_r / (N - dofminus)
    if (is.null(weights)) {
      ZwZ_r <- crossprod(Z_r)
    } else {
      ZwZ_r <- crossprod(Z_r, weights * Z_r)
    }
    Omega_r <- sigma_r_sq * ZwZ_r / N
  } else {
    Omega_r <- .compute_omega(Z_r, e_r, weights, cluster_vec, N,
                               dofminus = dofminus)
  }

  # --- J_r: J statistic of restricted model ---
  J_r <- .compute_j_with_omega(Z_r, X, y, Omega_r, weights, N)
  if (is.na(J_r)) {
    return(list(stat = NA_real_, p = NA_real_, df = as.integer(q),
                test_name = "Endogeneity",
                tested_vars = endog_vars))
  }

  # --- Extract submatrix matching full model's Z columns ---
  # Z_r = [Z, tested_endo], so first L columns correspond to Z
  z_cols <- seq_len(L)
  Omega_sub <- Omega_r[z_cols, z_cols, drop = FALSE]

  # --- J_f: re-estimate full model with restricted S ---
  J_f <- .compute_j_with_omega(Z, X, y, Omega_sub, weights, N)
  if (is.na(J_f)) {
    return(list(stat = NA_real_, p = NA_real_, df = as.integer(q),
                test_name = "Endogeneity",
                tested_vars = endog_vars))
  }

  # --- C-statistic ---
  C <- J_r - J_f

  # Collinearity df check: overid_df_r should be overid_df + q
  overid_df <- L - K
  overid_df_r <- L_r - K  # same K (regressors unchanged)
  if (overid_df_r - overid_df != q) {
    return(list(stat = 0, p = NA_real_, df = 0L,
                test_name = "Endogeneity",
                tested_vars = endog_vars))
  }

  p <- stats::pchisq(C, df = q, lower.tail = FALSE)

  list(
    stat = C,
    p = p,
    df = as.integer(q),
    test_name = "Endogeneity",
    tested_vars = endog_vars
  )
}
