# --------------------------------------------------------------------------
# ivreg2 pipeline: estimator dispatch
# Called only by ivreg2(); see R/ivreg2.R for the orchestrator.
# --------------------------------------------------------------------------

#' Determine estimation path, build omega_fn, dispatch to fitter.
#'
#' Handles wmatrix/smatrix routing, auto-bandwidth (before estimation for GMM,
#' after for non-GMM), omega_fn closure construction, and estimation dispatch
#' (7 branches: OLS, 2SLS, LIML/kclass, GMM2S, CUE, gmmw, smatrix-gmm2s).
#'
#' @return Named list with fit, method (possibly promoted), bw (resolved),
#'   wmatrix (possibly nulled), omega_fn.
#' @noRd
.resolve_and_fit <- function(parsed, prep, opts) {
  method <- opts$method
  vcov <- opts$vcov
  kernel <- opts$kernel
  small <- opts$small
  dofminus <- opts$dofminus
  fuller <- opts$fuller
  kclass <- opts$kclass
  weight_type <- opts$weight_type
  center <- opts$center
  psd <- opts$psd
  ivar <- opts$ivar
  # name-normalized in .prepare_model (.match_user_matrix)
  wmatrix <- prep$wmatrix
  smatrix <- prep$smatrix

  sdofminus   <- prep$sdofminus
  b0          <- prep$b0
  cluster_vec <- prep$cluster_vec
  time_index  <- prep$time_index
  bw          <- prep$bw

  # --- Determine effective estimation path for wmatrix/smatrix ---
  use_wmatrix <- !is.null(wmatrix)
  use_smatrix <- !is.null(smatrix)

  # wmatrix without robust/cluster/kernel/gmm2s → warn & ignore (Stata behavior)
  if (use_wmatrix && is.null(cluster_vec) && vcov == "iid" &&
      is.null(kernel) && method != "gmm2s") {
    warning("`wmatrix` is ignored without robust VCE, clustering, or HAC; ",
            "using standard IV estimation.", call. = FALSE)
    use_wmatrix <- FALSE
    wmatrix <- NULL
  }

  # smatrix alone implies efficient GMM if method is default.
  if (use_smatrix && !use_wmatrix && method == "2sls") {
    method <- "gmm2s"
  }

  # --- Build omega_fn closure ---
  omega_fn <- NULL
  needs_omega_fn <- method %in% c("gmm2s", "cue") || (use_wmatrix && !use_smatrix)

  if (needs_omega_fn) {
    # Resolve auto-bw using step-1 residuals (Stata ivreg2.ado:970-980)
    if (is.character(bw) && tolower(bw) == "auto") {
      if (!is.null(ivar)) {
        stop("Automatic bandwidth selection not available for panel data.",
             call. = FALSE)
      }
      if (use_wmatrix) {
        abw_resid <- .wmatrix_first_step_resid(parsed, wmatrix)
      } else {
        fit_step1 <- .fit_2sls(parsed, small = FALSE, dofminus = dofminus,
                               sdofminus = sdofminus)
        abw_resid <- fit_step1$residuals
      }
      bw <- .auto_bandwidth(
        resid = abw_resid,
        Z = parsed$Z,
        time_index = time_index,
        kernel = kernel,
        has_intercept = parsed$has_intercept,
        N = parsed$N
      )
      max_bw <- (time_index$T_span - 1) / time_index$tdelta
      if (bw > max_bw) {
        warning("Automatic bandwidth (", formatC(bw, format = "g"),
                ") exceeds time-span limit (", formatC(max_bw, format = "g"),
                "); capping at ", formatC(max_bw, format = "g"), ".",
                call. = FALSE)
        bw <- max_bw
      }
    }

    # Capture VCE parameters in closure
    gmm_Z <- parsed$Z
    gmm_w <- parsed$weights
    gmm_N <- parsed$N
    gmm_cv <- cluster_vec
    gmm_dm <- dofminus
    gmm_wt <- weight_type
    gmm_k  <- kernel
    gmm_bw <- bw
    gmm_ti <- time_index
    gmm_is_iid <- is.null(cluster_vec) && vcov == "iid" && is.null(kernel)
    gmm_center <- center
    gmm_psd <- psd
    gmm_vcov_type <- vcov
    gmm_sw <- isTRUE(opts$sw)
    gmm_ivar_vec <- prep$ivar_vec
    if (vcov == "AC") {
      if (!is.null(parsed$weights)) {
        gmm_ZwZ <- crossprod(parsed$Z, parsed$weights * parsed$Z)
      } else {
        gmm_ZwZ <- crossprod(parsed$Z)
      }
    } else {
      gmm_ZwZ <- NULL
    }
    omega_fn <- function(resid) {
      .compute_moment_cov(gmm_Z, resid, gmm_w, gmm_cv, gmm_N,
                          iid = gmm_is_iid,
                          dofminus = gmm_dm, weight_type = gmm_wt,
                          kernel = gmm_k, bw = gmm_bw, time_index = gmm_ti,
                          center = gmm_center, psd = gmm_psd,
                          vcov_type = gmm_vcov_type, ZwZ = gmm_ZwZ,
                          sw = gmm_sw, ivar_vec = gmm_ivar_vec)
    }
    # Structural rank bound of the Omega the GMM fitters must invert; only
    # meaningful for this closure (user-supplied smatrix paths pass Inf and
    # keep the numeric checks). See .cluster_rank_bound.
    omega_rank_bound <- .cluster_rank_bound(cluster_vec, kernel, center,
                                            parsed$weights)
  }

  # --- Estimation dispatch ---
  if (use_smatrix && !use_wmatrix) {
    omega_fn <- function(resid) smatrix
    fit <- .fit_gmm2s(parsed, small = small, dofminus = dofminus,
                      sdofminus = sdofminus, omega_fn = omega_fn)

  } else if (use_wmatrix && method == "gmm2s") {
    wstep_resid <- .wmatrix_first_step_resid(parsed, wmatrix)
    if (use_smatrix) {
      omega_fn <- function(resid) smatrix
      omega_rank_bound <- Inf
    } else {
      omega_fn_base <- omega_fn
      omega_fn <- function(resid) omega_fn_base(wstep_resid)
    }
    fit <- .fit_gmm2s(parsed, small = small, dofminus = dofminus,
                      sdofminus = sdofminus, omega_fn = omega_fn,
                      omega_rank_bound = omega_rank_bound)

  } else if (use_wmatrix) {
    if (use_smatrix) {
      omega_fn <- function(resid) smatrix
      omega_rank_bound <- Inf
    }
    fit <- .fit_gmm_wmatrix(parsed, small = small, dofminus = dofminus,
                             sdofminus = sdofminus,
                             W = wmatrix, omega_fn = omega_fn,
                             omega_rank_bound = omega_rank_bound)
    method <- "gmmw"

  } else if (method == "cue") {
    fit <- .fit_cue(parsed, small = small, dofminus = dofminus,
                    sdofminus = sdofminus, omega_fn = omega_fn, b0 = b0,
                    iid = gmm_is_iid,
                    omega_rank_bound = omega_rank_bound)

  } else if (method == "gmm2s") {
    fit <- .fit_gmm2s(parsed, small = small, dofminus = dofminus,
                      sdofminus = sdofminus, omega_fn = omega_fn,
                      omega_rank_bound = omega_rank_bound)

  } else if (method %in% c("liml", "kclass")) {
    fit <- .fit_kclass(parsed, method = method, kclass = kclass, fuller = fuller,
                small = small, dofminus = dofminus, sdofminus = sdofminus)
  } else if (parsed$is_iv && parsed$K1 > 0L) {
    fit <- .fit_2sls(parsed, small = small, dofminus = dofminus,
              sdofminus = sdofminus)
  } else {
    # 1-part OLS, or empty-endogenous (K1 = 0): estimation is exact OLS —
    # the excluded instruments enter only the diagnostics (Stata posts
    # e(model) = "ols" for `(=z)` models, ivreg2.ado:2101-2106).
    fit <- .fit_ols(parsed, small = small, dofminus = dofminus,
             sdofminus = sdofminus)
  }

  # --- Resolve bw = "auto" (Newey-West 1994) ---
  # For non-GMM: resolve AFTER estimation (need residuals).
  # For GMM/CUE: already resolved above.
  if (!method %in% c("gmm2s", "gmmw", "cue") && is.character(bw) && tolower(bw) == "auto") {
    if (!is.null(ivar)) {
      stop("Automatic bandwidth selection not available for panel data.",
           call. = FALSE)
    }
    abw_Z <- if (is.null(parsed$Z)) parsed$X else parsed$Z
    bw <- .auto_bandwidth(
      resid = fit$residuals,
      Z = abw_Z,
      time_index = time_index,
      kernel = kernel,
      has_intercept = parsed$has_intercept,
      N = parsed$N
    )
    max_bw <- (time_index$T_span - 1) / time_index$tdelta
    if (bw > max_bw) {
      warning("Automatic bandwidth (", formatC(bw, format = "g"),
              ") exceeds time-span limit (", formatC(max_bw, format = "g"),
              "); capping at ", formatC(max_bw, format = "g"), ".",
              call. = FALSE)
      bw <- max_bw
    }
  }

  list(fit = fit, method = method, bw = bw, wmatrix = wmatrix,
       omega_fn = omega_fn)
}
