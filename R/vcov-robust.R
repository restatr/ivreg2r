# --------------------------------------------------------------------------
# .psd_correct
# --------------------------------------------------------------------------
#' Apply PSD correction to a symmetric matrix
#'
#' Performs eigenvalue-based correction to ensure a matrix is positive
#' semi-definite.  Two modes are available:
#' \itemize{
#'   \item `"psd0"`: set negative eigenvalues to zero (Politis 2007).
#'   \item `"psda"`: replace negative eigenvalues with their absolute values
#'     (Stock & Watson 2008, Remark 8).
#' }
#'
#' When correction is applied, a warning is emitted indicating how many
#' eigenvalues were corrected and which mode was used.
#'
#' @param mat Symmetric matrix to correct.
#' @param psd Character: `"psd0"` or `"psda"`, or `NULL` (no correction).
#' @return Corrected symmetric matrix (or `mat` unchanged if `psd` is `NULL`
#'   or no negative eigenvalues are found).
#' @keywords internal
.psd_correct <- function(mat, psd) {
  if (is.null(psd)) return(mat)
  eig <- eigen(mat, symmetric = TRUE)
  # Use a relative tolerance to avoid triggering on floating-point noise
  tol <- .Machine$double.eps^0.5 * max(abs(eig$values))
  neg_idx <- eig$values < -tol
  if (!any(neg_idx)) return(mat)
  n_neg <- sum(neg_idx)
  warning("The covariance matrix was not positive semidefinite; ", n_neg,
          " negative eigenvalue", if (n_neg > 1L) "s",
          " corrected via the '", psd, "' method.", call. = FALSE)
  if (psd == "psd0") {
    eig$values <- pmax(eig$values, 0)
  } else {
    eig$values <- abs(eig$values)
  }
  mat_corrected <- eig$vectors %*% (eig$values * t(eig$vectors))
  (mat_corrected + t(mat_corrected)) / 2
}


# --------------------------------------------------------------------------
# .vcov_from_omega
# --------------------------------------------------------------------------
#' Assemble a sandwich VCV from the L x L moment covariance Omega
#'
#' Stata assembles every plain (non-GMM) robust-family VCV from the L x L
#' instrument-space moment covariance returned by `m_omega`, by congruence
#' with the first-stage map A (`X_hat = Z A`): `V = 1/N * aux5 * omega * aux5'`
#' (s_iegmm, ivreg2.ado:5556-5558). In unscaled R quantities that is
#' `V = N * bread * (A' Omega A) * bread`, followed by the small-sample
#' factor. This route is used when `psd` is set, so the psd correction is
#' applied to Omega *before* the congruence — `psd(A' S A) != A' psd(S) A`,
#' which is why correcting the K x K meat (or the final VCV) diverges from
#' Stata whenever the correction binds.
#'
#' The single formula reproduces every plain path exactly because each
#' VCE type's Omega normalization cancels its sandwich scaling:
#' HC and HAC use `Omega = meat/(N - dofminus)` with
#' `V = bread meat bread * N/(N - dofminus)`; cluster, two-way, DK, and
#' cluster+kernel use `Omega = meat/N` with `V = bread meat bread`; AC and
#' Stock-Watson Omegas are already per-observation with
#' `V = bread Omega bread * N`.
#'
#' @param bread K x K bread matrix (k-class bread for LIML/k-class).
#' @param A L x K first-stage coefficient matrix, or `NULL` for OLS
#'   (where Z = X and the congruence is the identity).
#' @param omega L x L moment covariance from [.compute_moment_cov()]
#'   (psd-corrected there when `psd` is set).
#' @param N Integer: number of observations.
#' @param K Integer: number of regressors.
#' @param M Integer: number of clusters (required when `cluster = TRUE`).
#' @param small Logical: apply the finite-sample correction.
#' @param dofminus,sdofminus Integer DoF adjustments.
#' @param cluster Logical: `TRUE` for the cluster-family small-sample factor
#'   `(N-1)/(N-K-sdofminus) * M/(M-1)`; `FALSE` for the robust-family factor
#'   `(N-dofminus)/(N-K-dofminus-sdofminus)`.
#' @return K x K variance-covariance matrix.
#' @keywords internal
.vcov_from_omega <- function(bread, A, omega, N, K, M = NULL, small = FALSE,
                             dofminus = 0L, sdofminus = 0L, cluster = FALSE) {
  meat <- if (is.null(A)) omega else crossprod(A, omega %*% A)
  V <- bread %*% meat %*% bread * N
  V <- (V + t(V)) / 2

  if (small) {
    V <- if (cluster) {
      V * ((N - 1) / (N - K - sdofminus)) * (M / (M - 1))
    } else {
      V * ((N - dofminus) / (N - K - dofminus - sdofminus))
    }
  }

  colnames(V) <- rownames(V) <- colnames(bread)
  V
}


# --------------------------------------------------------------------------
# .hc_meat
# --------------------------------------------------------------------------
#' Compute HC meat matrix with weight-type dispatch
#'
#' For aweight/pweight: meat = (w * X * e)' (w * X * e) = X' diag(w^2 e^2) X.
#' For fweight: meat = X' diag(w * e^2) X (linear in weights, not quadratic).
#' For unweighted: meat = X' diag(e^2) X.
#'
#' @param basis N x K matrix (X_hat for IV, X for OLS).
#' @param resid N-vector of residuals.
#' @param weights Normalized weights or NULL.
#' @param weight_type Character: `"aweight"`, `"fweight"`, or `"pweight"`.
#' @param center Logical: if `TRUE`, subtract the mean of moment conditions
#'   before computing the outer product. The mean computation depends on
#'   weight type, matching Stata's `m_omega()` (livreg2.do lines 179-188).
#' @return K x K symmetric meat matrix.
#' @keywords internal
.hc_meat <- function(basis, resid, weights = NULL, weight_type = "aweight",
                     center = FALSE) {
  if (!center) {
    # Uncentered paths (unchanged)
    if (is.null(weights)) return(crossprod(basis * resid))
    if (weight_type == "fweight") {
      return(crossprod(basis, weights * resid^2 * basis))
    }
    return(crossprod(weights * basis * resid))
  }
  # Centered path: form unweighted scores, compute weighted mean, center
  g <- basis * resid          # N x P unweighted scores
  N <- nrow(g)
  if (is.null(weights)) {
    gbar <- colMeans(g)
    g_c <- sweep(g, 2L, gbar)
    return(crossprod(g_c))
  }
  if (weight_type == "fweight") {
    N_eff <- sum(weights)                 # effective N = sum(w) for fweight
    gbar <- colSums(weights * g) / N_eff  # Stata: wv = wvar (linear)
    g_c <- sweep(g, 2L, gbar)
    return(crossprod(g_c, weights * g_c)) # sum w_i (g_c)(g_c)'
  }
  # aweight/pweight: Stata uses wv = (w*wf)^2; N = nrow (weights sum to N)
  gbar <- colSums(weights^2 * g) / N
  g_c <- sweep(g, 2L, gbar)
  crossprod(weights * g_c)               # sum w_i^2 (g_c)(g_c)'
}


# --------------------------------------------------------------------------
# .cl_scores
# --------------------------------------------------------------------------
#' Compute cluster scores (weight-type-agnostic), optionally centered
#'
#' Cluster scores are `weights * basis * resid` for all weight types.
#' The definition of `weights` (normalized for aweight/pweight, raw for fweight)
#' makes this expression correct for all types.
#'
#' When `center = TRUE`, subtracts the (weighted) mean of scores before
#' returning, so all downstream cluster/kernel aggregation operates on
#' already-centered scores.
#'
#' @param basis N x K matrix.
#' @param resid N-vector of residuals.
#' @param weights Normalized weights or NULL.
#' @param center Logical: if `TRUE`, center scores by subtracting their mean.
#' @param weight_type Character: `"aweight"`, `"fweight"`, or `"pweight"`.
#' @return N x K score matrix.
#' @keywords internal
.cl_scores <- function(basis, resid, weights = NULL,
                       center = FALSE, weight_type = "aweight") {
  scores <- if (is.null(weights)) basis * resid else weights * basis * resid
  if (center) {
    if (is.null(weights)) {
      gbar <- colMeans(scores)
    } else if (weight_type == "fweight") {
      gbar <- colSums(scores) / sum(weights)  # effective N = sum(w) for fweight
    } else {
      N <- nrow(scores)
      gbar <- colSums(weights * scores) / N   # sum(w^2 * g) / N
    }
    scores <- sweep(scores, 2L, gbar)
  }
  scores
}


# --------------------------------------------------------------------------
# .compute_hc_vcov
# --------------------------------------------------------------------------
#' Compute heteroskedasticity-consistent (robust) VCV
#'
#' Sandwich estimator using the score-based approach: meat = X' diag(e^2) X,
#' where X is `X_hat` for IV or `X` for OLS. The base normalization is
#' `N/(N-dofminus)`. When `small = TRUE`, applies the additional
#' finite-sample correction `(N-dofminus)/(N-K-dofminus-sdofminus)`,
#' giving a net scaling of `N/(N-K-dofminus-sdofminus)`.
#'
#' @param bread K x K bread matrix: \eqn{(X'P_Z X)^{-1}} for IV,
#'   \eqn{(X'X)^{-1}} for OLS.
#' @param X_hat N x K matrix: projected regressors for IV, original X for OLS.
#' @param resid Length-N residual vector.
#' @param N Integer: number of observations.
#' @param K Integer: number of regressors.
#' @param small Logical: if `TRUE`, apply the finite-sample correction
#'   `(N-dofminus)/(N-K-dofminus-sdofminus)`.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @return K x K variance-covariance matrix.
#' @keywords internal
.compute_hc_vcov <- function(bread, X_hat, resid, N, K,
                              small = FALSE, dofminus = 0L, sdofminus = 0L,
                              weights = NULL, weight_type = "aweight",
                              center = FALSE) {
  omega <- .hc_meat(X_hat, resid, weights, weight_type, center = center)
  V <- bread %*% omega %*% bread      # K x K sandwich
  V <- (V + t(V)) / 2                 # enforce symmetry

  # Base normalization (always applied)
  V <- V * (N / (N - dofminus))

  # Small-sample correction (Stata ivreg2.ado lines 1182-1189)
  if (small) {
    V <- V * ((N - dofminus) / (N - K - dofminus - sdofminus))
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
                              dofminus = 0L, sdofminus = 0L,
                              weights = NULL,
                              center = FALSE, weight_type = "aweight") {
  scores <- .cl_scores(X_hat, resid, weights,
                       center = center, weight_type = weight_type)
  omega <- .cluster_meat(scores, cluster_vec)    # K x K (unscaled)
  V <- bread %*% omega %*% bread                # sandwich
  V <- (V + t(V)) / 2                           # enforce symmetry

  if (small) {
    V <- V * ((N - 1) / (N - K - sdofminus)) * (M / (M - 1))
  }

  colnames(V) <- rownames(V) <- colnames(bread)
  V
}
