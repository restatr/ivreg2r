# --------------------------------------------------------------------------
# ivreg2 pipeline: covariance construction and small-sample corrections
# Called only by ivreg2(); see R/ivreg2.R for the orchestrator.
# --------------------------------------------------------------------------

#' Apply VCV corrections and compute final variance-covariance matrix.
#'
#' GMM path: apply small-sample corrections to fit$vcov, recompute sigma.
#' Non-GMM path: select bread (k-class vs 2SLS for COVIV), dispatch to VCV
#' helpers. When `psd` is set, plain robust-family paths assemble the VCV
#' from the psd-corrected L x L moment covariance via [.vcov_from_omega()]
#' (Stata m_omega parity); the iid VCV is never psd-corrected.
#'
#' @return Updated fit with final vcov, sigma, df.residual.
#' @noRd
.compute_vcov <- function(fit, parsed, prep, opts, method, bw) {
  vcov <- opts$vcov
  kernel <- opts$kernel
  small <- opts$small
  dofminus <- opts$dofminus
  coviv <- opts$coviv
  weight_type <- opts$weight_type
  center <- opts$center
  psd <- opts$psd

  sdofminus   <- prep$sdofminus
  cluster_vec <- prep$cluster_vec
  time_index  <- prep$time_index
  M           <- prep$M

  if (method %in% c("gmm2s", "gmmw", "cue")) {
    needs_vcov_correction <- small
    if (needs_vcov_correction) {
      if (!is.null(cluster_vec)) {
        fit$vcov <- fit$vcov * ((parsed$N - 1) / (parsed$N - parsed$K - sdofminus)) *
          (M / (M - 1))
        fit$df.residual <- as.integer(M - 1L)
      } else {
        fit$vcov <- fit$vcov * ((parsed$N - dofminus) /
                                  (parsed$N - parsed$K - dofminus - sdofminus))
      }
    } else if (!is.null(cluster_vec)) {
      fit$df.residual <- as.integer(M - 1L)
    }
    if (small) {
      fit$sigma <- sqrt(fit$rss / (parsed$N - parsed$K - dofminus - sdofminus))
    }
    colnames(fit$vcov) <- rownames(fit$vcov) <- names(fit$coefficients)
  } else {

  # For LIML/kclass, select bread: k-class bread by default, 2SLS bread if coviv
  bread_vcov <- if (method %in% c("liml", "kclass") && !coviv) {
    fit$bread_kclass
  } else {
    fit$bread
  }

  # COVIV + IID: override the k-class IID VCV with 2SLS-bread IID VCV
  if (method %in% c("liml", "kclass") && coviv &&
      vcov == "iid" && is.null(cluster_vec)) {
    fit$vcov <- fit$sigma^2 * fit$bread
    colnames(fit$vcov) <- rownames(fit$vcov) <- names(fit$coefficients)
  }

  if (isTRUE(opts$sw) || !is.null(cluster_vec) ||
      vcov %in% c("robust", "HAC", "AC")) {
    # K1 = 0 (empty-endogenous) models are estimated by OLS: the meat basis
    # is X itself, like the 1-part OLS path (Stata's robust VCV for `(=z)`
    # models never involves Z).
    X_hat_vcov <- if (parsed$is_iv && parsed$K1 > 0L) fit$X_hat else parsed$X
    resid_vcov <- fit$residuals
  }

  if (!is.null(psd) &&
      (isTRUE(opts$sw) || !is.null(cluster_vec) ||
       vcov %in% c("robust", "HAC", "AC"))) {
    # F2 (D9b): Stata applies the psd0/psda correction to the L x L moment
    # covariance S inside m_omega (livreg2.do:607-617) and assembles the
    # plain non-GMM VCV from the corrected S by congruence with the
    # first-stage map A (s_iegmm, ivreg2.ado:5556). Correcting the K x K
    # meat or the final VCV diverges whenever the correction binds, so with
    # psd set we route every plain robust-family path (HC, cluster, two-way,
    # HAC, AC, DK, cluster+kernel, SW) through the corrected S. This S is
    # the same computation as the stored fit$S.
    omega_vcov <- .fresh_moment_cov(fit, parsed, prep, opts, bw, center)
    is_cluster_family <- !is.null(cluster_vec)
    # First-stage map A with X_hat = Z A. For K1 = 0 (OLS estimation with
    # surplus instruments), X is a subset of Z's columns, so A is the exact
    # 0/1 selection matrix — no solve needed.
    A_vcov <- if (!parsed$is_iv) {
      NULL
    } else if (parsed$K1 > 0L) {
      fit$proj_coef
    } else {
      A_sel <- matrix(0, nrow = parsed$L, ncol = parsed$K,
                      dimnames = list(colnames(parsed$Z), colnames(parsed$X)))
      A_sel[cbind(match(colnames(parsed$X), colnames(parsed$Z)),
                  seq_len(parsed$K))] <- 1
      A_sel
    }
    fit$vcov <- .vcov_from_omega(
      bread = bread_vcov,
      A = A_vcov,
      omega = omega_vcov,
      N = parsed$N, K = parsed$K, M = M, small = small,
      dofminus = dofminus, sdofminus = sdofminus,
      cluster = is_cluster_family
    )
    if (is_cluster_family) {
      fit$df.residual <- as.integer(M - 1L)
    }
    # Stash for .compute_stored_sw(): fit$S must be this exact matrix, and
    # recomputing it would emit a duplicate psd-correction warning.
    fit$psd_moment_cov <- omega_vcov
  } else if (isTRUE(opts$sw)) {
    fit$vcov <- .compute_sw_vcov(
      bread = bread_vcov, X_hat = X_hat_vcov, resid = resid_vcov,
      ivar_vec = prep$ivar_vec, N = parsed$N, K = parsed$K,
      small = small, dofminus = dofminus, sdofminus = sdofminus,
      weights = parsed$weights, weight_type = weight_type,
      center = center
    )
  } else if (!is.null(cluster_vec) && !is.null(kernel)) {
    is_twoway <- is.list(cluster_vec)
    fit$vcov <- .compute_cluster_kernel_vcov(
      bread = bread_vcov, X_hat = X_hat_vcov, resid = resid_vcov,
      cluster_vec = cluster_vec, time_index = time_index,
      kernel = kernel, bw = bw,
      N = parsed$N, K = parsed$K, M = M, small = small,
      dofminus = dofminus, sdofminus = sdofminus,
      weights = parsed$weights, weight_type = weight_type,
      is_twoway = is_twoway,
      center = center
    )
    fit$df.residual <- as.integer(M - 1L)
  } else if (vcov == "HAC") {
    fit$vcov <- .compute_hac_vcov(bread_vcov, X_hat_vcov, resid_vcov,
                                   time_index, kernel, bw,
                                   parsed$N, parsed$K,
                                   dofminus = dofminus, sdofminus = sdofminus,
                                   small = small,
                                   weights = parsed$weights,
                                   weight_type = weight_type,
                                   center = center)
  } else if (vcov == "AC") {
    fit$vcov <- .compute_ac_vcov(bread_vcov, X_hat_vcov, resid_vcov,
                                  time_index, kernel, bw,
                                  parsed$N, parsed$K,
                                  dofminus = dofminus, sdofminus = sdofminus,
                                  small = small,
                                  weights = parsed$weights,
                                  weight_type = weight_type)
  } else if (!is.null(cluster_vec)) {
    fit$vcov <- .compute_cl_vcov(bread_vcov, X_hat_vcov, resid_vcov,
                                  cluster_vec, parsed$N, parsed$K, M, small,
                                  dofminus = dofminus, sdofminus = sdofminus,
                                  weights = parsed$weights,
                                  weight_type = weight_type,
                                  center = center)
    fit$df.residual <- as.integer(M - 1L)
  } else if (vcov == "robust") {
    fit$vcov <- .compute_hc_vcov(bread_vcov, X_hat_vcov, resid_vcov,
                                  parsed$N, parsed$K,
                                  small = small, dofminus = dofminus,
                                  sdofminus = sdofminus,
                                  weights = parsed$weights,
                                  weight_type = weight_type,
                                  center = center)
  }

  }  # end of non-GMM VCV block

  # No final-VCV psd pass: the correction lives at the S level (the
  # psd branch above for plain paths; omega_fn for GMM paths), and the
  # iid VCV is never corrected, matching Stata (m_omega is not involved
  # in the iid VCV; psd is silently inert there).

  fit
}
