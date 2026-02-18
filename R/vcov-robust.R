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
#' @param vcov_type Character: `"HC0"` or `"HC1"`. HC1 always applies N/(N-K).
#' @return K x K variance-covariance matrix.
#' @keywords internal
.compute_hc_vcov <- function(bread, X_hat, resid, N, K, vcov_type) {
  scores <- X_hat * resid              # N x K
  omega <- crossprod(scores)           # K x K meat
  V <- bread %*% omega %*% bread      # K x K sandwich
  V <- (V + t(V)) / 2                 # enforce symmetry

  if (vcov_type == "HC1") {
    V <- V * (N / (N - K))
  }

  colnames(V) <- rownames(V) <- colnames(bread)
  V
}


# --------------------------------------------------------------------------
# .compute_cl_vcov
# --------------------------------------------------------------------------
#' Compute one-way cluster-robust VCV
#'
#' Sandwich estimator with scores aggregated within clusters via
#' [rowsum()].  When `small = TRUE`, applies the Stata `ivreg2` correction
#' `(N-1)/(N-K) * M/(M-1)`.
#'
#' @param bread K x K bread matrix: \eqn{(X'P_Z X)^{-1}} for IV,
#'   \eqn{(X'X)^{-1}} for OLS.
#' @param X_hat N x K matrix: projected regressors for IV, original X for OLS.
#' @param resid Length-N residual vector.
#' @param cluster_vec Length-N vector of cluster identifiers (need not be
#'   sequential integers).
#' @param N Integer: number of observations.
#' @param K Integer: number of regressors.
#' @param M Integer: number of clusters.
#' @param small Logical: if `TRUE`, apply the finite-sample correction
#'   `(N-1)/(N-K) * M/(M-1)`.
#' @return K x K variance-covariance matrix.
#' @keywords internal
.compute_cl_vcov <- function(bread, X_hat, resid, cluster_vec, N, K, M, small) {
  scores <- X_hat * resid                       # N x K
  cluster_scores <- rowsum(scores, cluster_vec)  # M x K
  omega <- crossprod(cluster_scores)             # K x K (unscaled)
  V <- bread %*% omega %*% bread                # sandwich
  V <- (V + t(V)) / 2                           # enforce symmetry

  if (small) {
    V <- V * ((N - 1) / (N - K)) * (M / (M - 1))
  }

  colnames(V) <- rownames(V) <- colnames(bread)
  V
}
