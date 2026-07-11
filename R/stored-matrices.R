# --------------------------------------------------------------------------
# ivreg2 pipeline: stored moment matrices (e(S)/e(W) analogues)
# Called only by ivreg2(); see R/ivreg2.R for the orchestrator.
# --------------------------------------------------------------------------

#' Invert a symmetric positive-definite matrix with graceful failure.
#'
#' Tries Cholesky (chol2inv), falls back to qr.solve, and returns NULL with
#' a warning if the matrix is computationally singular. Used for the stored
#' weighting matrix W, where Stata's hard error (exit 506) is relaxed to the
#' package's warn-and-skip convention.
#'
#' @return Inverse matrix, or NULL if singular.
#' @noRd
.safe_inverse <- function(A, what = "matrix") {
  inv <- tryCatch(chol2inv(chol(A)), error = function(e) NULL)
  if (is.null(inv)) {
    inv <- tryCatch(qr.solve(A), error = function(e) NULL)
  }
  if (is.null(inv)) {
    warning("The GMM weighting matrix (`$W`) could not be computed because ",
            what, " is singular; `$W` is set to NULL.", call. = FALSE)
  }
  inv
}

#' Compute the L x L moment covariance S from the final residuals.
#'
#' Thin wrapper that assembles the [.compute_moment_cov()] arguments from the
#' fit/parsed/prep/opts objects (including the AC path's Z'WZ precompute).
#' Shared by the stored `fit$S` and the psd-corrected VCV assembly in
#' `.compute_vcov()`, so the two are the same matrix by construction.
#' psd correction happens inside [.compute_moment_cov()] when `opts$psd`
#' is set.
#'
#' Must be called BEFORE the unsort step: HAC/AC paths index residuals by the
#' sorted `prep$time_index`.
#'
#' @return L x L symmetric matrix S (no dimnames).
#' @noRd
.fresh_moment_cov <- function(fit, parsed, prep, opts, bw, center) {
  Zmat <- if (parsed$is_iv) parsed$Z else parsed$X
  is_iid <- is.null(prep$cluster_vec) && opts$vcov == "iid" &&
    is.null(opts$kernel)
  ZwZ <- NULL
  if (opts$vcov == "AC") {
    ZwZ <- if (is.null(parsed$weights)) {
      crossprod(Zmat)
    } else {
      crossprod(Zmat, parsed$weights * Zmat)
    }
  }
  .compute_moment_cov(
    Z = Zmat, residuals = fit$residuals, weights = parsed$weights,
    cluster_vec = prep$cluster_vec, N = parsed$N, iid = is_iid,
    dofminus = opts$dofminus, weight_type = opts$weight_type,
    kernel = opts$kernel, bw = bw, time_index = prep$time_index,
    center = center, psd = opts$psd, vcov_type = opts$vcov,
    ZwZ = ZwZ, sw = isTRUE(opts$sw), ivar_vec = prep$ivar_vec
  )
}

#' Compute the stored moment-covariance S and weighting matrix W.
#'
#' Builds the matrices posted on the fitted object as `fit$S` and `fit$W`,
#' matching Stata's `e(S)` and `e(W)` (ivreg2.ado:1905-1918):
#'
#' - S: user `smatrix` is echoed (never recomputed); GMM methods reuse the
#'   defining Omega already computed by the fitter (`fit$omega`); all other
#'   paths (OLS, 2SLS, LIML/k-class, any VCE) compute S fresh from the final
#'   residuals via `.fresh_moment_cov()`.
#' - W: NULL for LIML/k-class ("No weighting matrix defined", ivreg2.ado:1912);
#'   user `wmatrix` is echoed (even under gmm2s, where it is the first-step
#'   W); gmm2s/CUE/b0/smatrix-implied-GMM use `solve(S)`; 1-step OLS/2SLS use
#'   `(1/sigma^2) (Z'WZ/N)^{-1}` from the final residuals with the
#'   large-sample `sigma^2 = RSS/(N - dofminus)` (s_gmm1s, ivreg2.ado:5287-5371;
#'   independent of `small`).
#'
#' Must be called BEFORE the unsort step: HAC/AC paths index residuals by the
#' sorted `prep$time_index`.
#'
#' @return list(S = matrix, W = matrix or NULL), dimnames = instrument names.
#' @noRd
.compute_stored_sw <- function(fit, parsed, prep, opts, est, bw, center) {
  method <- est$method
  Zmat <- if (parsed$is_iv) parsed$Z else parsed$X
  inames <- colnames(Zmat)

  # --- S ---
  if (!is.null(prep$smatrix)) {
    S <- prep$smatrix
  } else if (method %in% c("gmm2s", "gmmw", "cue")) {
    S <- fit$omega
  } else if (!is.null(fit$psd_moment_cov)) {
    # Reuse the corrected Omega the psd VCV assembly already computed
    # (.compute_vcov): same residuals, same correction — recomputing would
    # only duplicate the psd warning.
    S <- fit$psd_moment_cov
  } else {
    S <- .fresh_moment_cov(fit, parsed, prep, opts, bw, center)
  }
  dimnames(S) <- list(inames, inames)

  # --- W ---
  if (method %in% c("liml", "kclass")) {
    W <- NULL
  } else if (!is.null(est$wmatrix)) {
    W <- est$wmatrix
  } else if (method %in% c("gmm2s", "gmmw", "cue")) {
    W <- .safe_inverse(S, what = "the moment covariance")
  } else {
    ZWZ <- if (is.null(parsed$weights)) {
      crossprod(Zmat)
    } else {
      crossprod(Zmat, parsed$weights * Zmat)
    }
    sigma2 <- fit$rss / (parsed$N - opts$dofminus)
    QZZinv <- .safe_inverse(ZWZ / parsed$N,
                            what = "the instrument cross-product matrix")
    W <- if (is.null(QZZinv)) NULL else QZZinv / sigma2
  }
  if (!is.null(W)) dimnames(W) <- list(inames, inames)

  list(S = S, W = W)
}

#' Condition numbers of X'X and Z'Z, and the Gaussian log-likelihood.
#'
#' Extracted verbatim from the ivreg2() orchestrator; expressions and their
#' order are unchanged so the stored values are bit-identical.
#' @noRd
.compute_cond_ll <- function(parsed, fit) {
  # Condition numbers: Stata's cond(XX) = max(eigenvalue)/min(eigenvalue)
  # which equals kappa(XX, exact = TRUE) in R.
  if (is.null(parsed$weights)) {
    XX <- crossprod(parsed$X)
  } else {
    XX <- crossprod(sqrt(parsed$weights) * parsed$X)
  }
  condxx <- kappa(XX, exact = TRUE)

  if (parsed$is_iv) {
    if (is.null(parsed$weights)) {
      ZZ <- crossprod(parsed$Z)
    } else {
      ZZ <- crossprod(sqrt(parsed$weights) * parsed$Z)
    }
    condzz <- kappa(ZZ, exact = TRUE)
  } else {
    # OLS: Z = X, so condzz = condxx (matches Stata behavior)
    condzz <- condxx
  }

  # Gaussian log-likelihood (Stata line 2009)
  ll <- -0.5 * (parsed$N * log(2 * pi) + parsed$N * log(fit$rss / parsed$N) +
                 parsed$N)

  list(condxx = condxx, condzz = condzz, ll = ll)
}
