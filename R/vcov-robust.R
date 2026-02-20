# --------------------------------------------------------------------------
# .compute_hc_vcov
# --------------------------------------------------------------------------
#' Compute HC0 or HC1 heteroskedasticity-consistent VCV
#'
#' Sandwich estimator using the score-based approach: meat = X' diag(e^2) X,
#' where X is `X_hat` for IV or `X` for OLS. HC1 always applies the N/(N-K)
#' finite-sample correction; HC0 never does.
#'
#' @param bread K x K bread matrix: \eqn{(X'P_Z X)^{-1}} for IV,
#'   \eqn{(X'X)^{-1}} for OLS.
#' @param X_hat N x K matrix: projected regressors for IV, original X for OLS.
#' @param resid Length-N residual vector.
#' @param N Integer: number of observations.
#' @param K Integer: number of regressors.
#' @param vcov_type Character: `"HC0"` or `"HC1"`. HC1 applies
#'   `N/(N-K-dofminus-sdofminus)`. HC0 applies `N/(N-dofminus)`.
#' @param small Logical: whether small-sample corrections are applied.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @return K x K variance-covariance matrix.
#' @keywords internal
.compute_hc_vcov <- function(bread, X_hat, resid, N, K, vcov_type,
                              small = FALSE, dofminus = 0L, sdofminus = 0L) {
  scores <- X_hat * resid              # N x K
  omega <- crossprod(scores)           # K x K meat
  V <- bread %*% omega %*% bread      # K x K sandwich
  V <- (V + t(V)) / 2                 # enforce symmetry

  if (vcov_type == "HC1") {
    V <- V * (N / (N - K - dofminus - sdofminus))
  } else {
    # HC0: omega normalization adjustment (Stata divides meat by N-dofminus)
    V <- V * (N / (N - dofminus))
  }

  colnames(V) <- rownames(V) <- colnames(bread)
  V
}


# --------------------------------------------------------------------------
# .cluster_meat
# --------------------------------------------------------------------------
#' Compute clustered meat matrix (one-way or two-way)
#'
#' Shared helper for the `rowsum() -> crossprod()` pattern used in VCV
#' computation and all diagnostic score covariance computations.
#'
#' For one-way clustering, computes `crossprod(rowsum(scores, cluster_vec))`.
#' For two-way clustering (Cameron-Gelbach-Miller 2006), computes
#' `meat1 + meat2 - meat_intersection` where the intersection uses
#' `interaction(cv1, cv2, drop = TRUE)`.
#'
#' @param scores N x P score matrix (P = K for VCV, P = L for diagnostics).
#' @param cluster_vec Length-N vector (one-way) or list of 2 vectors (two-way)
#'   of cluster identifiers.
#' @return P x P symmetric meat matrix (unscaled).
#' @keywords internal
.cluster_meat <- function(scores, cluster_vec) {
  if (is.list(cluster_vec)) {
    int_vec <- interaction(cluster_vec[[1]], cluster_vec[[2]], drop = TRUE)
    meat1 <- crossprod(rowsum(scores, cluster_vec[[1]], reorder = FALSE))
    meat2 <- crossprod(rowsum(scores, cluster_vec[[2]], reorder = FALSE))
    meat3 <- crossprod(rowsum(scores, int_vec, reorder = FALSE))
    meat <- meat1 + meat2 - meat3
  } else {
    meat <- crossprod(rowsum(scores, cluster_vec, reorder = FALSE))
  }
  (meat + t(meat)) / 2
}


# --------------------------------------------------------------------------
# .compute_cl_vcov
# --------------------------------------------------------------------------
#' Compute cluster-robust VCV (one-way or two-way)
#'
#' Sandwich estimator with scores aggregated within clusters via
#' [.cluster_meat()].  When `small = TRUE`, applies the Stata `ivreg2`
#' correction `(N-1)/(N-K-sdofminus) * M/(M-1)`.
#'
#' @param bread K x K bread matrix: \eqn{(X'P_Z X)^{-1}} for IV,
#'   \eqn{(X'X)^{-1}} for OLS.
#' @param X_hat N x K matrix: projected regressors for IV, original X for OLS.
#' @param resid Length-N residual vector.
#' @param cluster_vec Length-N vector (one-way) or list of 2 vectors (two-way)
#'   of cluster identifiers.
#' @param N Integer: number of observations.
#' @param K Integer: number of regressors.
#' @param M Integer: number of clusters (min(M1, M2) for two-way).
#' @param small Logical: if `TRUE`, apply the finite-sample correction
#'   `(N-1)/(N-K-sdofminus) * M/(M-1)`.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#'   Note: dofminus does NOT appear in cluster VCV scaling (Stata convention).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @return K x K variance-covariance matrix.
#' @keywords internal
.compute_cl_vcov <- function(bread, X_hat, resid, cluster_vec, N, K, M, small,
                              dofminus = 0L, sdofminus = 0L) {
  scores <- X_hat * resid                       # N x K
  omega <- .cluster_meat(scores, cluster_vec)    # K x K (unscaled)
  V <- bread %*% omega %*% bread                # sandwich
  V <- (V + t(V)) / 2                           # enforce symmetry

  if (small) {
    V <- V * ((N - 1) / (N - K - sdofminus)) * (M / (M - 1))
  }

  colnames(V) <- rownames(V) <- colnames(bread)
  V
}
