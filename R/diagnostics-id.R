# --------------------------------------------------------------------------
# Identification diagnostics (Ticket D2)
# --------------------------------------------------------------------------
# Anderson LM / Cragg-Donald F (IID) and Kleibergen-Paap rk LM / rk Wald F
# (robust/cluster) underidentification and weak identification tests.
#
# Ported from Stata's ranktest.ado (Kleibergen-Paap 2006) and ivreg2.ado.


# --------------------------------------------------------------------------
# .sym_sqrt
# --------------------------------------------------------------------------
#' Symmetric matrix square root via eigendecomposition
#'
#' Computes the symmetric square root of a PSD matrix via eigendecomposition.
#' Negative eigenvalues (numerical noise) are clamped to zero with a warning.
#'
#' @param A Symmetric matrix.
#' @return Symmetric square root matrix.
#' @keywords internal
.sym_sqrt <- function(A) {
  eig <- eigen(A, symmetric = TRUE)
  d <- eig$values
  if (any(d < 0)) {
    warning("Negative eigenvalues clamped to 0 in symmetric square root.",
            call. = FALSE)
    d[d < 0] <- 0
  }
  eig$vectors %*% (sqrt(d) * t(eig$vectors))
}


# --------------------------------------------------------------------------
# .partial_fwl
# --------------------------------------------------------------------------
#' FWL: partial out exogenous regressors from endogenous and excluded IVs
#'
#' Returns the residuals from projecting X1 and Z1 onto X2, where X2 is the exogenous regressor
#' matrix (including intercept if present).
#'
#' @param X1 N x K1 endogenous regressor matrix.
#' @param Z1 N x L1 excluded instrument matrix.
#' @param X2 N x K2 exogenous regressor matrix (including intercept).
#' @param weights Normalized weights (sum to N), or NULL.
#' @return List with `X1_perp` and `Z1_perp`.
#' @keywords internal
.partial_fwl <- function(X1, Z1, X2, weights) {
  if (is.null(weights)) {
    X1_perp <- X1 - X2 %*% lm.fit(X2, X1)$coefficients
    Z1_perp <- Z1 - X2 %*% lm.fit(X2, Z1)$coefficients
  } else {
    X1_perp <- X1 - X2 %*% lm.wfit(X2, X1, weights)$coefficients
    Z1_perp <- Z1 - X2 %*% lm.wfit(X2, Z1, weights)$coefficients
  }
  list(X1_perp = X1_perp, Z1_perp = Z1_perp)
}


# --------------------------------------------------------------------------
# .canonical_correlations
# --------------------------------------------------------------------------
#' SVD-based canonical correlations and intermediates
#'
#' Computes the Cholesky-standardized reduced-form coefficient matrix
#' theta = R_zz pihat irQyy and its full SVD, returning all intermediates
#' needed for both Anderson/CD (IID) and Kleibergen-Paap (robust) tests.
#'
#' @param X1_perp N x K1 partialled endogenous regressors.
#' @param Z1_perp N x L1 partialled excluded instruments.
#' @param N Number of observations.
#' @param K1 Number of endogenous regressors.
#' @param L1 Number of excluded instruments.
#' @param weights Normalized weights or NULL.
#' @return List with theta, U, cc, V, eval, pihat, irQyy, irQzz,
#'   or NULL if Cholesky fails.
#' @keywords internal
.canonical_correlations <- function(X1_perp, Z1_perp, N, K1, L1, weights) {
  # Cross-products (weighted if needed)
  if (is.null(weights)) {
    Qzz <- crossprod(Z1_perp) / N
    Qzy <- crossprod(Z1_perp, X1_perp) / N
    Qyy <- crossprod(X1_perp) / N
  } else {
    Qzz <- crossprod(Z1_perp, weights * Z1_perp) / N
    Qzy <- crossprod(Z1_perp, weights * X1_perp) / N
    Qyy <- crossprod(X1_perp, weights * X1_perp) / N
  }

  # Cholesky factorizations — R's chol() returns upper triangular R
  # such that R'R = A, i.e. A = R_zz' %*% R_zz
  R_zz <- tryCatch(chol(Qzz), error = function(e) NULL)
  R_yy <- tryCatch(chol(Qyy), error = function(e) NULL)
  if (is.null(R_zz) || is.null(R_yy)) {
    warning("Cholesky factorization failed in canonical correlations; ",
            "identification tests not computed.", call. = FALSE)
    return(NULL)
  }

  # Inverse of upper Cholesky factors
  irQzz <- backsolve(R_zz, diag(L1))
  irQyy <- backsolve(R_yy, diag(K1))

  # Reduced-form coefficient: pihat = Qzz^{-1} Qzy
  pihat <- chol2inv(R_zz) %*% Qzy

  # Standardized matrix: theta = R_zz %*% pihat %*% irQyy
  # Equivalently: t(R_zz) %*% pihat %*% irQyy (using Stata's lower convention)
  # Since R_zz is upper: R_zz' is lower. Stata uses cholesky() = lower.
  # theta = rQzz' * pihat * irQyy in Stata, where rQzz = cholesky(Qzz) = R_zz'
  # So: theta = R_zz %*% pihat %*% irQyy (R_zz is upper = Stata's rQzz')
  theta <- R_zz %*% pihat %*% irQyy

  # Full SVD
  sv <- svd(theta, nu = L1, nv = K1)
  cc <- sv$d[seq_len(min(K1, L1))]

  # Guard: max canonical correlation = 1 indicates collinearity
  if (any(cc >= 1)) {
    warning("Canonical correlation = 1 detected (collinearity); ",
            "identification tests not computed.", call. = FALSE)
    return(NULL)
  }

  eval <- cc^2

  list(
    theta  = theta,
    U      = sv$u,          # L1 x L1
    cc     = cc,            # min(K1,L1) vector
    V      = sv$v,          # K1 x K1
    eval   = eval,          # squared canonical correlations
    pihat  = pihat,         # L1 x K1
    irQyy  = irQyy,        # K1 x K1
    irQzz  = irQzz         # L1 x L1
  )
}


# --------------------------------------------------------------------------
# .anderson_lm_test
# --------------------------------------------------------------------------
#' Anderson canonical correlation LM test (IID underidentification)
#'
#' @param cc_result Output from `.canonical_correlations()`.
#' @param N Number of observations.
#' @param L1 Number of excluded instruments.
#' @param K1 Number of endogenous regressors.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @return Named list with stat, p, df, test_name.
#' @keywords internal
.anderson_lm_test <- function(cc_result, N, L1, K1, dofminus = 0L) {
  df <- as.integer(L1 - K1 + 1L)
  stat <- (N - dofminus) * min(cc_result$eval)
  p <- stats::pchisq(stat, df = df, lower.tail = FALSE)
  list(stat = stat, p = p, df = df, test_name = "Anderson canon. corr. LM statistic")
}


# --------------------------------------------------------------------------
# .cragg_donald_f
# --------------------------------------------------------------------------
#' Cragg-Donald Wald F statistic (weak identification)
#'
#' CD F = min(eval/(1-eval)) * (N - L) / L1
#' where eval are squared canonical correlations.
#' This equals the minimum eigenvalue of the CD matrix scaled by
#' the appropriate degrees of freedom.
#'
#' @param cc_result Output from `.canonical_correlations()`.
#' @param N Number of observations.
#' @param L Total number of instruments (including exogenous + intercept).
#' @param L1 Number of excluded instruments.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @return Named list with stat, test_name.
#' @keywords internal
.cragg_donald_f <- function(cc_result, N, L, L1, dofminus = 0L,
                             sdofminus = 0L) {
  min_eval <- min(cc_result$eval)
  cd <- min_eval / (1 - min_eval)
  stat <- cd * (N - L - dofminus - sdofminus) / L1
  list(stat = stat, test_name = "Cragg-Donald Wald F statistic")
}


# --------------------------------------------------------------------------
# .kp_omega
# --------------------------------------------------------------------------
#' Kronecker score covariance for Kleibergen-Paap test
#'
#' Computes the `(K1*L1) x (K1*L1)` score covariance matrix shat0.
#' For HC: shat0 = crossprod(scores) / N where row i of scores is
#'   kron(V_hat_i, Z1_perp_i) for each observation i.
#' For cluster: same but with rowsum over clusters first.
#'
#' @param Z1_perp N x L1 partialled excluded instruments.
#' @param V_hat N x K1 matrix (X1_perp for LM, residuals for Wald).
#' @param weights Normalized weights or NULL.
#' @param cluster_vec Cluster vector or NULL.
#' @param N Number of observations.
#' @param K1 Number of endogenous regressors.
#' @param L1 Number of excluded instruments.
#' @return `(K1*L1) x (K1*L1)` symmetric matrix.
#' @keywords internal
.kp_omega <- function(Z1_perp, V_hat, weights, cluster_vec, N, K1, L1,
                       weight_type = "aweight",
                       kernel = NULL, bw = NULL, time_index = NULL) {
  # Build N x (K1*L1) score matrix: row i = kron(V_hat[i,], Z1_perp[i,])
  # For K1=1: reduces to V_hat * Z1_perp (scalar broadcast)
  if (K1 == 1L) {
    scores <- drop(V_hat) * Z1_perp
  } else {
    scores <- matrix(0, nrow = N, ncol = K1 * L1)
    for (j in seq_len(K1)) {
      cols <- ((j - 1L) * L1 + 1L):(j * L1)
      scores[, cols] <- V_hat[, j] * Z1_perp
    }
  }

  # Cluster+kernel, cluster, HAC, or HC aggregation with weight-type dispatch
  if (!is.null(cluster_vec) && !is.null(kernel)) {
    # Cluster + kernel (DK or Thompson)
    if (!is.null(weights)) scores <- weights * scores
    if (is.list(cluster_vec)) {
      # Thompson: CGM decomposition — third term is HAC (not HC),
      # matching Stata livreg2.do line 336.
      shat1 <- crossprod(rowsum(scores, cluster_vec[[1L]], reorder = FALSE))
      shat1 <- (shat1 + t(shat1)) / 2
      shat2 <- .cluster_kernel_scores_meat(scores, time_index, kernel, bw)
      shat3 <- .hac_scores_meat(scores, time_index, kernel, bw)
      shat0 <- (shat1 + shat2 - shat3) / N
    } else {
      # DK: one-way cluster+kernel on tvar
      shat0 <- .cluster_kernel_scores_meat(scores, time_index, kernel, bw) / N
    }
  } else if (!is.null(cluster_vec)) {
    if (!is.null(weights)) scores <- weights * scores
    shat0 <- .cluster_meat(scores, cluster_vec) / N
  } else if (!is.null(kernel)) {
    # HAC path: use pre-formed scores with lag-loop autocovariances
    if (!is.null(weights)) scores <- weights * scores
    shat0 <- .hac_scores_meat(scores, time_index, kernel, bw) / N
  } else {
    if (is.null(weights)) {
      shat0 <- crossprod(scores) / N
    } else if (weight_type == "fweight") {
      shat0 <- crossprod(sqrt(weights) * scores) / N
    } else {
      shat0 <- crossprod(weights * scores) / N
    }
  }
  (shat0 + t(shat0)) / 2  # force symmetry
}


# --------------------------------------------------------------------------
# .kp_rk_stat
# --------------------------------------------------------------------------
#' Kleibergen-Paap rk statistic (LM or Wald variant)
#'
#' Computes the rk chi-squared statistic from the canonical correlation
#' intermediates and the Kronecker score covariance.
#'
#' @param cc_result Output from `.canonical_correlations()`.
#' @param shat0 Score covariance from `.kp_omega()`.
#' @param N Number of observations.
#' @param K1 Number of endogenous regressors.
#' @param L1 Number of excluded instruments.
#' @return Named list with chi2, df.
#' @keywords internal
.kp_rk_stat <- function(cc_result, shat0, N, K1, L1) {
  U <- cc_result$U        # L1 x L1
  V <- cc_result$V        # K1 x K1
  theta <- cc_result$theta # L1 x K1
  irQyy <- cc_result$irQyy
  irQzz <- cc_result$irQzz

  # Kronecker transform of score covariance
  # kpvar = (irQyy' kron irQzz') %*% shat0 %*% t(irQyy' kron irQzz')
  # In R: t(irQyy) is K1xK1, t(irQzz) is L1xL1
  kron_transform <- kronecker(t(irQyy), t(irQzz))
  kpvar <- kron_transform %*% shat0 %*% t(kron_transform)
  kpvar <- (kpvar + t(kpvar)) / 2  # force symmetry

  # SVD partition at rank kk = K1 - 1
  kk <- K1 - 1L

  if (kk == 0L) {
    # K1 = 1: no partitioning needed — u12/v12 are empty
    # u22 = full U (L1 x L1), v22 = full V (K1 x K1) = scalar 1x1
    u22 <- U
    v22 <- matrix(V, nrow = K1, ncol = K1)
    u12 <- matrix(0, nrow = 0, ncol = ncol(u22))
    v12 <- matrix(0, nrow = 0, ncol = ncol(v22))
  } else {
    # Extract blocks from U (L1 x L1) and V (K1 x K1)
    u12 <- U[seq_len(kk), (kk + 1L):L1, drop = FALSE]
    u22 <- U[(kk + 1L):L1, (kk + 1L):L1, drop = FALSE]
    v12 <- V[seq_len(kk), (kk + 1L):K1, drop = FALSE]
    v22 <- V[(kk + 1L):K1, (kk + 1L):K1, drop = FALSE]
  }

  # Symmetric square roots of u22 %*% t(u22) and v22 %*% t(v22)
  u22_half <- .sym_sqrt(u22 %*% t(u22))
  v22_half <- .sym_sqrt(v22 %*% t(v22))

  # Transformation matrices
  # aq = rbind(u12, u22) %*% solve(u22) %*% u22_half
  # bq = v22_half %*% solve(t(v22)) %*% t(rbind(v12, v22))
  u_full <- rbind(u12, u22)
  v_full <- rbind(v12, v22)
  aq <- u_full %*% solve(u22) %*% u22_half
  bq <- v22_half %*% solve(t(v22)) %*% t(v_full)

  # Test statistic
  vecthat <- as.numeric(theta)  # vec(theta), column-major
  kron_bq_aq <- kronecker(bq, t(aq))
  lab <- kron_bq_aq %*% vecthat
  vlab <- kron_bq_aq %*% kpvar %*% t(kron_bq_aq)
  vlab <- (vlab + t(vlab)) / 2  # force symmetry

  # Invert vlab — use pseudo-inverse if rank deficient
  vlab_inv <- tryCatch(solve(vlab), error = function(e) NULL)
  rrank <- 0L
  if (is.null(vlab_inv)) {
    eig <- eigen(vlab, symmetric = TRUE)
    tol <- max(eig$values) * .Machine$double.eps * max(dim(vlab))
    pos <- eig$values > tol
    rrank <- sum(!pos)
    if (rrank > 0L) {
      warning("vlab is rank-deficient (rank deficit = ", rrank,
              "); using pseudoinverse.", call. = FALSE)
    }
    vlab_inv <- eig$vectors[, pos, drop = FALSE] %*%
      (t(eig$vectors[, pos, drop = FALSE]) / eig$values[pos])
  }

  chi2 <- drop(N * crossprod(lab, vlab_inv %*% lab))
  df_base <- (L1 - kk) * 1L  # single equation: ii = 1
  df <- as.integer(df_base - rrank)

  list(chi2 = chi2, df = df)
}


# --------------------------------------------------------------------------
# .compute_id_tests
# --------------------------------------------------------------------------
#' Dispatcher for underidentification and weak identification tests
#'
#' Computes Anderson LM + Cragg-Donald F (always), plus Kleibergen-Paap
#' rk LM and rk Wald F when VCE is robust or clustered.
#'
#' @param X N x K regressor matrix.
#' @param Z N x L instrument matrix.
#' @param y N x 1 response vector.
#' @param residuals N x 1 residual vector from 2SLS.
#' @param weights Normalized weights or NULL.
#' @param cluster_vec Cluster vector or NULL.
#' @param vcov_type Character: "iid", "HC0", "HC1", or "CL".
#' @param N,K,L,K1,L1 Integer dimensions.
#' @param endo_names Character vector of endogenous variable names.
#' @param excluded_names Character vector of excluded instrument names.
#' @param has_intercept Logical.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @return List with underid, weak_id, weak_id_robust (or NULL).
#' @keywords internal
.compute_id_tests <- function(X, Z, y, residuals, weights, cluster_vec,
                              vcov_type, N, K, L, K1, L1, M = NULL,
                              endo_names, excluded_names, has_intercept,
                              dofminus = 0L, sdofminus = 0L,
                              weight_type = "aweight",
                              kernel = NULL, bw = NULL,
                              time_index = NULL) {

  # Top-level guard: catch unexpected errors
  result <- tryCatch({

    # --- Extract sub-matrices ---
    # X1 = endogenous columns of X (N x K1)
    # X2 = exogenous columns of X (N x K2)
    # Z1 = excluded instruments from Z (N x L1)
    X_cols <- colnames(X)
    Z_cols <- colnames(Z)

    endo_idx <- match(endo_names, X_cols)
    exog_idx <- setdiff(seq_len(ncol(X)), endo_idx)
    excl_idx <- match(excluded_names, Z_cols)

    X1 <- X[, endo_idx, drop = FALSE]
    X2 <- X[, exog_idx, drop = FALSE]
    Z1 <- Z[, excl_idx, drop = FALSE]

    # --- FWL partialling ---
    fwl <- .partial_fwl(X1, Z1, X2, weights)
    X1_perp <- fwl$X1_perp
    Z1_perp <- fwl$Z1_perp

    # --- Canonical correlations ---
    cc_result <- .canonical_correlations(X1_perp, Z1_perp, N, K1, L1, weights)
    if (is.null(cc_result)) {
      return(list(
        underid = list(stat = NA_real_, p = NA_real_,
                       df = as.integer(L1 - K1 + 1L),
                       test_name = if (vcov_type == "iid")
                         "Anderson canon. corr. LM statistic"
                       else "Kleibergen-Paap rk LM statistic"),
        weak_id = list(stat = NA_real_,
                       test_name = "Cragg-Donald Wald F statistic"),
        weak_id_robust = if (vcov_type != "iid")
          list(stat = NA_real_,
               test_name = "Kleibergen-Paap rk Wald F statistic")
        else NULL
      ))
    }

    # --- Always compute Cragg-Donald F ---
    weak_id <- .cragg_donald_f(cc_result, N, L, L1,
                                dofminus = dofminus, sdofminus = sdofminus)

    # --- IID path (no kernel) ---
    # AC with kernel falls through to the robust/HAC path below because
    # Stata's ivreg2 passes bwopt to ranktest (which triggers the non-IID
    # path in ranktest), so KP stats are computed even for AC.
    if (vcov_type == "iid" && is.null(kernel)) {
      underid <- .anderson_lm_test(cc_result, N, L1, K1,
                                    dofminus = dofminus)
      return(list(
        underid        = underid,
        weak_id        = weak_id,
        weak_id_robust = NULL
      ))
    }

    # --- Robust / cluster / HAC path ---
    # Reduced-form residuals for Wald variant:
    # V_wald = X1_perp - Z1_perp %*% pihat
    V_lm <- X1_perp
    V_wald <- X1_perp - Z1_perp %*% cc_result$pihat

    # Stata's ivreg2 has a bug where `kernopt` (the kernel option passed to
    # ranktest) is never captured from the vkernel subroutine — only `kernel`
    # (for VCV) and `bwopt` (for bandwidth) are captured. So ranktest always
    # defaults to Bartlett kernel regardless of the user-specified kernel.
    # We match this behavior for Stata parity: KP tests always use Bartlett.
    kp_kernel <- if (!is.null(kernel)) "Bartlett" else NULL

    # AC vs HAC/HC: different omega computation for KP test.
    # For AC (vcov_type=="AC"), Stata's ranktest enters the homoskedastic
    # block of m_omega (robust="" but kernel is set), using the Kronecker
    # structure sigma_tau # ZZ_tau. For HAC/HC, it uses score cross-products.
    if (vcov_type == "AC") {
      # AC omega via Kronecker structure (K1=1 case)
      ZwZ_perp <- if (is.null(weights)) {
        crossprod(Z1_perp)
      } else {
        crossprod(Z1_perp, weights * Z1_perp)
      }
      shat0_lm <- .ac_meat(Z1_perp, drop(V_lm), time_index, kp_kernel, bw,
                            N, dofminus, weights, weight_type, ZwZ_perp)
      shat0_wald <- .ac_meat(Z1_perp, drop(V_wald), time_index, kp_kernel, bw,
                              N, dofminus, weights, weight_type, ZwZ_perp)
    } else {
      # HAC / HC / cluster path: use score cross-products
      shat0_lm <- .kp_omega(Z1_perp, V_lm, weights, cluster_vec, N, K1, L1,
                              weight_type = weight_type,
                              kernel = kp_kernel, bw = bw,
                              time_index = time_index)
      shat0_wald <- .kp_omega(Z1_perp, V_wald, weights, cluster_vec, N, K1, L1,
                                weight_type = weight_type,
                                kernel = kp_kernel, bw = bw,
                                time_index = time_index)
    }

    # KP rk LM
    kp_lm <- .kp_rk_stat(cc_result, shat0_lm, N, K1, L1)

    # Underid: KP rk LM chi-squared
    # Non-cluster: stat = chi2/N * (N - dofminus) (Stata line 1671)
    # Cluster: stat = chi2 (no adjustment, Stata line 1675)
    if (is.null(cluster_vec)) {
      underid_stat <- kp_lm$chi2 / N * (N - dofminus)
    } else {
      underid_stat <- kp_lm$chi2
    }
    underid_df <- kp_lm$df
    underid_p <- stats::pchisq(underid_stat, df = underid_df, lower.tail = FALSE)
    underid <- list(stat = underid_stat, p = underid_p, df = underid_df,
                    test_name = "Kleibergen-Paap rk LM statistic")

    # KP rk Wald
    kp_wald <- .kp_rk_stat(cc_result, shat0_wald, N, K1, L1)

    # Convert KP Wald chi-sq to F:
    # Non-cluster: F = chi2/N * (N - L - dofminus - sdofminus) / L1
    # Cluster:     F = chi2/(N-1) * (N - L - sdofminus) * (M-1)/M / L1
    if (is.null(cluster_vec)) {
      wald_f <- kp_wald$chi2 / N * (N - L - dofminus - sdofminus) / L1
    } else {
      wald_f <- kp_wald$chi2 / (N - 1) * (N - L - sdofminus) *
        (M - 1) / M / L1
    }
    weak_id_robust <- list(stat = wald_f,
                           test_name = "Kleibergen-Paap rk Wald F statistic")

    list(
      underid        = underid,
      weak_id        = weak_id,
      weak_id_robust = weak_id_robust
    )

  }, error = function(e) {
    warning("Identification test computation failed: ", conditionMessage(e),
            call. = FALSE)
    list(
      underid = list(stat = NA_real_, p = NA_real_,
                     df = as.integer(L1 - K1 + 1L),
                     test_name = if (vcov_type == "iid")
                       "Anderson canon. corr. LM statistic"
                     else "Kleibergen-Paap rk LM statistic"),
      weak_id = list(stat = NA_real_,
                     test_name = "Cragg-Donald Wald F statistic"),
      weak_id_robust = if (vcov_type != "iid")
        list(stat = NA_real_,
             test_name = "Kleibergen-Paap rk Wald F statistic")
      else NULL
    )
  })

  result
}
