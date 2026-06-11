# --------------------------------------------------------------------------
# Instrument redundancy test (Ticket P1)
# --------------------------------------------------------------------------
# Tests H0: specified excluded instruments are redundant (zero explanatory
# power for endogenous regressors) conditional on maintained instruments.
# Equivalent to Stata's `redundant(varlist)` option in ivreg2.
#
# Uses the KP rk LM test of H0: rank=0 on the matrix
# Z_tested_perp' * X1_perp, where both are partialled of exogenous
# regressors AND maintained (non-tested) excluded instruments.
#
# For H0: rank=0 the full KP rk machinery simplifies:
#   IID: stat = (N - dofminus) * sum(ALL squared canonical correlations)
#   Robust: direct GMM quadratic form on vec(Z_tested' * X1 / N)


# --------------------------------------------------------------------------
# .compute_redundancy_test
# --------------------------------------------------------------------------
#' Instrument redundancy test (LM)
#'
#' Computes the KP rk LM statistic testing H0: rank=0 for the
#' first-stage relationship between endogenous regressors and tested
#' excluded instruments, conditional on maintained instruments.
#'
#' @param X N x K regressor matrix.
#' @param Z N x L instrument matrix.
#' @param weights Normalized weights (sum to N), or NULL.
#' @param cluster_vec Cluster membership vector, or NULL.
#' @param vcov_type Character: "iid", "robust", "HAC", "AC", or "CL".
#' @param N Number of observations.
#' @param K1 Number of endogenous regressors.
#' @param endo_colnames Character vector of endogenous regressor column names.
#' @param excluded_colnames Character vector of all excluded instrument column names.
#' @param redundant_vars Character vector of excluded instrument names to test.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param weight_type Character: weight type.
#' @param kernel Canonical kernel name or NULL.
#' @param bw Numeric bandwidth or NULL.
#' @param time_index Time index list or NULL.
#' @param center Logical: center scores (default FALSE).
#' @return Named list with `stat`, `p`, `df`, `test_name`, `tested_vars`.
#' @keywords internal
.compute_redundancy_test <- function(X, Z, weights, cluster_vec, vcov_type,
                                      N, K1, endo_colnames, excluded_colnames,
                                      redundant_vars, dofminus = 0L,
                                      weight_type = "aweight",
                                      kernel = NULL, bw = NULL,
                                      time_index = NULL,
                                      center = FALSE) {

  L_tested <- length(redundant_vars)
  df <- as.integer(K1 * L_tested)

  na_result <- list(stat = NA_real_, p = NA_real_, df = df,
                    test_name = "Redundancy test (LM)",
                    tested_vars = redundant_vars)

  result <- tryCatch({

    # --- Extract sub-matrices ---
    X_cols <- colnames(X)
    Z_cols <- colnames(Z)

    endo_idx <- match(endo_colnames, X_cols)
    exog_idx <- setdiff(seq_len(ncol(X)), endo_idx)

    X1 <- X[, endo_idx, drop = FALSE]   # N x K1 endogenous
    X2 <- X[, exog_idx, drop = FALSE]   # N x K2 exogenous

    tested_idx <- match(redundant_vars, Z_cols)
    excl_idx <- match(excluded_colnames, Z_cols)
    maintained_idx <- setdiff(excl_idx, tested_idx)

    Z_tested <- Z[, tested_idx, drop = FALSE]     # N x L_tested

    # Partial vars: exogenous regressors + maintained excluded instruments
    if (length(maintained_idx) > 0L) {
      Z_maintained <- Z[, maintained_idx, drop = FALSE]
      partial_vars <- cbind(X2, Z_maintained)
    } else {
      partial_vars <- X2
    }

    # --- FWL: partial out [X2, Z_maintained] from X1 and Z_tested ---
    fwl <- .partial_fwl(X1, Z_tested, partial_vars, weights)
    X1_perp <- fwl$X1_perp
    Z_tested_perp <- fwl$Z1_perp

    # --- IID path (vcov_type == "iid" and no kernel) ---
    # AC with kernel falls through to robust path (matching KP id test)
    if (vcov_type == "iid" && is.null(kernel)) {
      # Canonical correlations between X1_perp and Z_tested_perp
      cc_result <- .canonical_correlations(X1_perp, Z_tested_perp,
                                            N, K1, L_tested, weights)
      if (is.null(cc_result)) return(na_result)

      # H0: rank=0 → sum ALL squared canonical correlations
      stat <- (N - dofminus) * sum(cc_result$eval)
      p <- stats::pchisq(stat, df = df, lower.tail = FALSE)

      return(list(stat = stat, p = p, df = df,
                  test_name = "Redundancy test (LM)",
                  tested_vars = redundant_vars))
    }

    # --- Robust / cluster / HAC path ---
    # Direct GMM-style quadratic form (KP rk with rank=0 simplifies to this)

    # Moment condition: gbar = vec(Z_tested_perp' * X1_perp / N)
    if (is.null(weights)) {
      Qzy <- crossprod(Z_tested_perp, X1_perp) / N
    } else {
      Qzy <- crossprod(Z_tested_perp, weights * X1_perp) / N
    }
    gbar <- as.numeric(Qzy)  # vec(), column-major: (L_tested*K1) x 1

    # Score covariance: shat0 via .kp_omega()
    # V_hat = X1_perp for LM variant (testing H0: rank=0)
    # KP tests always use Bartlett kernel internally (matching Stata bug/behavior)
    kp_kernel <- if (!is.null(kernel)) "Bartlett" else NULL

    if (vcov_type == "AC") {
      # AC omega via Kronecker structure (same as existing KP test AC path)
      # For K1 == 1: use .ac_meat() directly
      # For K1 > 1: route through .kp_omega() with kernel (score-based HAC;
      #   algebraically equivalent to Kronecker AC for homoskedastic errors)
      if (K1 == 1L) {
        ZwZ_perp <- if (is.null(weights)) {
          crossprod(Z_tested_perp)
        } else {
          crossprod(Z_tested_perp, weights * Z_tested_perp)
        }
        shat0 <- .ac_meat(Z_tested_perp, drop(X1_perp), time_index,
                           kp_kernel, bw, N, dofminus, weights, weight_type,
                           ZwZ_perp)
      } else {
        shat0 <- .kp_omega(Z_tested_perp, X1_perp, weights, cluster_vec,
                            N, K1, L_tested,
                            weight_type = weight_type,
                            kernel = kp_kernel, bw = bw,
                            time_index = time_index,
                            center = center)
      }
    } else {
      shat0 <- .kp_omega(Z_tested_perp, X1_perp, weights, cluster_vec,
                          N, K1, L_tested,
                          weight_type = weight_type,
                          kernel = kp_kernel, bw = bw,
                          time_index = time_index,
                          center = center)
    }

    # Invert shat0 — use pseudo-inverse if rank deficient
    shat0_inv <- tryCatch(solve(shat0), error = function(e) NULL)
    rrank <- 0L
    if (is.null(shat0_inv)) {
      eig <- eigen(shat0, symmetric = TRUE)
      tol <- max(eig$values) * .Machine$double.eps * max(dim(shat0))
      pos <- eig$values > tol
      rrank <- sum(!pos)
      if (rrank > 0L) {
        warning("shat0 is rank-deficient (rank deficit = ", rrank,
                "); using pseudoinverse for redundancy test.", call. = FALSE)
      }
      shat0_inv <- eig$vectors[, pos, drop = FALSE] %*%
        (t(eig$vectors[, pos, drop = FALSE]) / eig$values[pos])
    }

    # GMM quadratic form: j = N * gbar' * inv(shat0) * gbar
    j <- drop(N * crossprod(gbar, shat0_inv %*% gbar))

    # DoF and stat adjustment
    if (is.null(cluster_vec)) {
      stat <- j / N * (N - dofminus)
    } else {
      stat <- j
    }

    df_adj <- as.integer(df - rrank)
    p <- stats::pchisq(stat, df = df_adj, lower.tail = FALSE)

    list(stat = stat, p = p, df = df_adj,
         test_name = "Redundancy test (LM)",
         tested_vars = redundant_vars)

  }, error = function(e) {
    warning("Redundancy test computation failed: ", conditionMessage(e),
            call. = FALSE)
    na_result
  })

  result
}
