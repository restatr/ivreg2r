# --------------------------------------------------------------------------
# HAC and AC VCV estimation (Ticket M2)
# --------------------------------------------------------------------------
# Implements HAC (heteroskedasticity and autocorrelation consistent) and AC
# (autocorrelation consistent) VCE types for time-series data.
#
# Reference: livreg2.do — m_omega() (lines 142–350).
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# .build_time_index
# --------------------------------------------------------------------------
#' Build time-index structure for HAC/AC estimation
#'
#' Computes sort order, time span, minimum time gap, and gap count from
#' time and panel variables.
#'
#' @param tvar_vec Numeric vector of time values (aligned to model frame rows).
#' @param ivar_vec Factor/character/numeric panel identifier vector, or NULL
#'   for pure time series.
#' @return List with components:
#'   \describe{
#'     \item{sort_order}{Integer permutation that sorts data by (ivar, tvar).}
#'     \item{tvar_sorted}{Time values in sorted order.}
#'     \item{ivar_sorted}{Panel IDs in sorted order, or NULL.}
#'     \item{T_span}{Integer: max(tvar) - min(tvar) + 1.}
#'     \item{tdelta}{Numeric: minimum nonzero gap between consecutive sorted
#'       time values (within panels if panel data).}
#'     \item{n_gaps}{Integer: number of missing time periods in the grid.}
#'     \item{panel_info}{List of per-panel start/end indices in sorted data,
#'       or NULL for pure time series.}
#'   }
#' @keywords internal
.build_time_index <- function(tvar_vec, ivar_vec = NULL) {
  N <- length(tvar_vec)

  # Validate unique time keys (Stata enforces this via tsset)
  if (is.null(ivar_vec)) {
    if (anyDuplicated(tvar_vec)) {
      stop("Time variable contains duplicate values. ",
           "HAC/AC estimation requires unique time values.",
           call. = FALSE)
    }
  } else {
    if (anyDuplicated(data.frame(ivar_vec, tvar_vec))) {
      stop("Duplicate (panel, time) combinations found. ",
           "HAC/AC estimation requires unique time values within each panel.",
           call. = FALSE)
    }
  }

  # Sort by (ivar, tvar)
  if (!is.null(ivar_vec)) {
    sort_order <- order(ivar_vec, tvar_vec)
  } else {
    sort_order <- order(tvar_vec)
  }

  tvar_sorted <- tvar_vec[sort_order]
  ivar_sorted <- if (!is.null(ivar_vec)) ivar_vec[sort_order] else NULL

  # T_span: max - min + 1 (global, matching Stata T=max(t)-min(t)+1)
  T_span <- as.integer(max(tvar_sorted) - min(tvar_sorted) + 1L)

  # tdelta: minimum nonzero gap between consecutive time values
  # Within panels if panel data
  if (!is.null(ivar_sorted)) {
    panels <- split(seq_len(N), ivar_sorted[seq_len(N)])
    all_gaps <- numeric(0)
    for (p in panels) {
      if (length(p) > 1L) {
        g <- diff(tvar_sorted[p])
        g <- g[g > 0]
        if (length(g) > 0L) all_gaps <- c(all_gaps, g)
      }
    }
  } else {
    all_gaps <- diff(tvar_sorted)
    all_gaps <- all_gaps[all_gaps > 0]
  }
  tdelta <- if (length(all_gaps) > 0L) min(all_gaps) else 1

  # n_gaps: count gaps in the time grid
  if (!is.null(ivar_sorted)) {
    n_gaps <- 0L
    for (p in panels) {
      if (length(p) > 1L) {
        tv <- tvar_sorted[p]
        expected <- (max(tv) - min(tv)) / tdelta
        actual <- length(p) - 1L
        n_gaps <- n_gaps + as.integer(expected - actual)
      }
    }
  } else {
    expected <- (max(tvar_sorted) - min(tvar_sorted)) / tdelta
    n_gaps <- as.integer(expected - (N - 1L))
  }

  # Panel info: list of start/end indices per panel
  panel_info <- NULL
  if (!is.null(ivar_sorted)) {
    panel_ids <- unique(ivar_sorted)
    panel_info <- vector("list", length(panel_ids))
    names(panel_info) <- as.character(panel_ids)
    for (i in seq_along(panel_ids)) {
      idx <- which(ivar_sorted == panel_ids[i])
      panel_info[[i]] <- list(start = min(idx), end = max(idx))
    }
  }

  list(
    sort_order  = sort_order,
    tvar_sorted = tvar_sorted,
    ivar_sorted = ivar_sorted,
    T_span      = T_span,
    tdelta      = tdelta,
    n_gaps      = n_gaps,
    panel_info  = panel_info
  )
}


# --------------------------------------------------------------------------
# .lag_pairs
# --------------------------------------------------------------------------
#' Find observation pairs separated by tau time periods
#'
#' For each observation at time t, finds the observation at time t - tau*tdelta.
#' Uses `match()` for O(N) lookup within panels.
#'
#' @param time_index List from `.build_time_index()`.
#' @param tau Integer lag (positive).
#' @return Two-column integer matrix `[i_now, i_lag]` with row indices into
#'   the **sorted** data. Zero rows if no matches exist.
#' @keywords internal
.lag_pairs <- function(time_index, tau) {
  tvar <- time_index$tvar_sorted
  tdelta <- time_index$tdelta
  N <- length(tvar)
  lag_gap <- tau * tdelta

  if (is.null(time_index$panel_info)) {
    # Pure time series: simple match on sorted time values
    target_times <- tvar - lag_gap
    lag_idx <- match(target_times, tvar)
    valid <- which(!is.na(lag_idx))
    if (length(valid) == 0L) {
      return(matrix(integer(0), ncol = 2L,
                    dimnames = list(NULL, c("i_now", "i_lag"))))
    }
    cbind(i_now = valid, i_lag = lag_idx[valid])
  } else {
    # Panel data: match within panels
    results <- vector("list", length(time_index$panel_info))
    for (j in seq_along(time_index$panel_info)) {
      pi <- time_index$panel_info[[j]]
      idx <- pi$start:pi$end
      tv <- tvar[idx]
      target_times <- tv - lag_gap
      lag_pos <- match(target_times, tv)
      valid <- which(!is.na(lag_pos))
      if (length(valid) > 0L) {
        results[[j]] <- cbind(
          i_now = idx[valid],
          i_lag = idx[lag_pos[valid]]
        )
      }
    }
    result <- do.call(rbind, results)
    if (is.null(result) || nrow(result) == 0L) {
      return(matrix(integer(0), ncol = 2L,
                    dimnames = list(NULL, c("i_now", "i_lag"))))
    }
    result
  }
}


# --------------------------------------------------------------------------
# .hac_meat
# --------------------------------------------------------------------------
#' Compute HAC meat matrix (heteroskedasticity and autocorrelation consistent)
#'
#' Returns the **unscaled** P x P meat matrix. The caller divides by
#' `(N - dofminus)` to get the HAC omega (Stata livreg2.do line 326).
#'
#' @param basis N x P matrix (X_hat for VCV, Z for diagnostics). Must be
#'   in sorted time order (matching `time_index`).
#' @param resid N-vector of residuals, in sorted time order.
#' @param time_index List from `.build_time_index()`.
#' @param kernel Canonical kernel name (from `.validate_kernel()`).
#' @param bw Numeric bandwidth.
#' @param weights Normalized weights (sorted) or NULL.
#' @param weight_type Character: `"aweight"` or `"pweight"` (fweight blocked).
#' @return P x P symmetric meat matrix (unscaled).
#' @keywords internal
.hac_meat <- function(basis, resid, time_index, kernel, bw,
                      weights = NULL, weight_type = "aweight") {
  P <- ncol(basis)

  # Diagonal (tau = 0): reuse existing HC meat
  shat <- .hc_meat(basis, resid, weights, weight_type)

  # Determine max lag TAU
  ktype <- .kernel_type(kernel)
  if (ktype == "spectral") {
    TAU <- as.integer(time_index$T_span / time_index$tdelta - 1)
  } else {
    TAU <- as.integer(floor(bw))
  }

  if (TAU < 1L) return(shat)

  # Off-diagonal loop (tau = 1..TAU)
  for (tau in seq_len(TAU)) {
    kw <- .kernel_weights(tau, bw, kernel)

    # Zero kernel weight: skip (Stata livreg2.do:209)
    if (kw == 0) next

    # Find observation pairs at this lag
    pairs <- .lag_pairs(time_index, tau)

    # Zero-row lag check (livreg2.do:218)
    if (nrow(pairs) == 0L) next

    i_now <- pairs[, 1L]
    i_lag <- pairs[, 2L]

    # Compute cross-product: sum of (w_t * e_t) * (w_{t-tau} * e_{t-tau}) * Z_t' Z_{t-tau}
    if (is.null(weights)) {
      wv <- resid[i_now] * resid[i_lag]
    } else {
      # aweight/pweight: quadratic weight product w_t * w_{t-tau}
      # Stata: wv = wvar[tmatrix[.,1]] :* wvar[tmatrix[.,2]] * (wf^2)
      # Our weights are already normalized (wf absorbed), so product is w_t * w_{t-tau}
      wv <- weights[i_now] * weights[i_lag] * resid[i_now] * resid[i_lag]
    }

    ghat <- crossprod(basis[i_now, , drop = FALSE],
                      wv * basis[i_lag, , drop = FALSE])

    shat <- shat + kw * (ghat + t(ghat))
  }

  shat
}


# --------------------------------------------------------------------------
# .hac_scores_meat
# --------------------------------------------------------------------------
#' Compute HAC meat from pre-formed score vectors
#'
#' Like `.hac_meat()` but takes pre-formed N x P scores directly, without
#' separate basis and residual vectors. Used by the Kleibergen-Paap path
#' where scores are Kronecker products that cannot be factored into
#' basis * resid form.
#'
#' @param scores N x P score matrix, in sorted time order.
#' @param time_index List from `.build_time_index()`.
#' @param kernel Canonical kernel name.
#' @param bw Numeric bandwidth.
#' @param weights Normalized weights (sorted) or NULL. For the KP path,
#'   weights are already incorporated into the scores.
#' @return P x P symmetric meat matrix (unscaled).
#' @keywords internal
.hac_scores_meat <- function(scores, time_index, kernel, bw,
                             weights = NULL) {
  P <- ncol(scores)

  # Diagonal (tau = 0)
  shat <- crossprod(scores)

  # Determine max lag TAU
  ktype <- .kernel_type(kernel)
  if (ktype == "spectral") {
    TAU <- as.integer(time_index$T_span / time_index$tdelta - 1)
  } else {
    TAU <- as.integer(floor(bw))
  }

  if (TAU < 1L) return((shat + t(shat)) / 2)

  # Off-diagonal loop
  for (tau in seq_len(TAU)) {
    kw <- .kernel_weights(tau, bw, kernel)
    if (kw == 0) next

    pairs <- .lag_pairs(time_index, tau)
    if (nrow(pairs) == 0L) next

    i_now <- pairs[, 1L]
    i_lag <- pairs[, 2L]

    ghat <- crossprod(scores[i_now, , drop = FALSE],
                      scores[i_lag, , drop = FALSE])
    shat <- shat + kw * (ghat + t(ghat))
  }

  (shat + t(shat)) / 2
}


# --------------------------------------------------------------------------
# .ac_meat
# --------------------------------------------------------------------------
#' Compute AC meat matrix (autocorrelation consistent, homoskedastic)
#'
#' Returns `shat / N` (the per-observation AC omega). The diagonal is
#' `sigma^2 * Z'WZ / N` where `sigma^2 = e'We / (N - dofminus)`.
#' The VCV wrapper multiplies by N: `V = N * bread * (shat/N) * bread`.
#' Diagnostics use `shat/N` directly as the omega (matching the IID
#' convention `sigma^2 * ZwZ / N`).
#'
#' The AC meat uses the Kronecker structure: `sigma_tau * Z'WZ_tau` at each
#' lag, where `sigma_tau` is a scalar autocovariance of the residuals and
#' `Z'WZ_tau` is the cross-product of instruments at lag tau.
#'
#' @param basis N x P matrix (X_hat for VCV, Z for diagnostics), sorted.
#' @param resid N-vector of residuals, sorted.
#' @param time_index List from `.build_time_index()`.
#' @param kernel Canonical kernel name.
#' @param bw Numeric bandwidth.
#' @param N Integer: number of observations.
#' @param dofminus Integer: large-sample DoF adjustment.
#' @param weights Normalized weights (sorted) or NULL.
#' @param weight_type Character: weight type.
#' @param ZwZ P x P precomputed instrument cross-product `Z'WZ`.
#' @return P x P matrix `shat / N` (normalized).
#' @keywords internal
.ac_meat <- function(basis, resid, time_index, kernel, bw,
                     N, dofminus, weights = NULL,
                     weight_type = "aweight", ZwZ) {
  # Diagonal (tau = 0): sigma^2 * ZwZ
  # sigma^2 = e'We / (N - dofminus)
  if (is.null(weights)) {
    sigma2 <- sum(resid^2) / (N - dofminus)
  } else {
    sigma2 <- sum(weights * resid^2) / (N - dofminus)
  }
  shat <- sigma2 * ZwZ

  # Determine max lag TAU
  ktype <- .kernel_type(kernel)
  if (ktype == "spectral") {
    TAU <- as.integer(time_index$T_span / time_index$tdelta - 1)
  } else {
    TAU <- as.integer(floor(bw))
  }

  if (TAU >= 1L) {
    for (tau in seq_len(TAU)) {
      kw <- .kernel_weights(tau, bw, kernel)
      if (kw == 0) next

      pairs <- .lag_pairs(time_index, tau)
      if (nrow(pairs) == 0L) next

      i_now <- pairs[, 1L]
      i_lag <- pairs[, 2L]

      # Scalar autocovariance: sigma_tau = sum(w_t * w_{t-tau} * e_t * e_{t-tau}) / (N - dofminus)
      # Stata: sigmahat = quadcross(e[tmatrix[.,1],.], wv, e[tmatrix[.,2],.]) / (N-dofminus)
      if (is.null(weights)) {
        sigma_tau <- sum(resid[i_now] * resid[i_lag]) / (N - dofminus)
      } else {
        wv <- weights[i_now] * weights[i_lag]
        sigma_tau <- sum(wv * resid[i_now] * resid[i_lag]) / (N - dofminus)
      }

      # Instrument cross-product at this lag: ZZ_tau = Z_now' * diag(wv) * Z_lag
      # Stata: ZZhat = quadcross(Z[tmatrix[.,1],.], wv, Z[tmatrix[.,2],.])
      if (is.null(weights)) {
        ZZ_tau <- crossprod(basis[i_now, , drop = FALSE],
                            basis[i_lag, , drop = FALSE])
      } else {
        wv2 <- weights[i_now] * weights[i_lag]
        ZZ_tau <- crossprod(basis[i_now, , drop = FALSE],
                            wv2 * basis[i_lag, , drop = FALSE])
      }

      # Kronecker product reduces to scalar * matrix for single equation
      ghat <- sigma_tau * ZZ_tau
      shat <- shat + kw * (ghat + t(ghat))
    }
  }

  shat / N
}


# --------------------------------------------------------------------------
# .compute_hac_vcov
# --------------------------------------------------------------------------
#' Compute HAC variance-covariance matrix
#'
#' Sandwich estimator: `V = bread * meat * bread * N/(N - dofminus)`.
#' This matches the HC scaling pattern used in `.compute_hc_vcov()`.
#' Note: Stata's m_omega divides by (N-dofminus) in Z-space, but our bread
#' is `(X'P_Z X)^{-1}` (unscaled), so we use the HC-style `N/(N-dofminus)`
#' multiplier instead.
#'
#' @param bread K x K bread matrix.
#' @param X_hat N x K projected regressor matrix (sorted by time).
#' @param resid N-vector of residuals (sorted by time).
#' @param time_index List from `.build_time_index()`.
#' @param kernel Canonical kernel name.
#' @param bw Numeric bandwidth.
#' @param N Integer: number of observations.
#' @param K Integer: number of regressors.
#' @param dofminus Integer: large-sample DoF adjustment.
#' @param sdofminus Integer: small-sample DoF adjustment.
#' @param weights Normalized weights (sorted) or NULL.
#' @param weight_type Character: weight type.
#' @return K x K variance-covariance matrix.
#' @keywords internal
.compute_hac_vcov <- function(bread, X_hat, resid, time_index, kernel, bw,
                              N, K, dofminus = 0L, sdofminus = 0L,
                              weights = NULL, weight_type = "aweight") {
  meat <- .hac_meat(X_hat, resid, time_index, kernel, bw,
                    weights, weight_type)
  V <- bread %*% meat %*% bread
  V <- (V + t(V)) / 2
  V <- V * (N / (N - dofminus))
  colnames(V) <- rownames(V) <- colnames(bread)
  V
}


# --------------------------------------------------------------------------
# .compute_ac_vcov
# --------------------------------------------------------------------------
#' Compute AC variance-covariance matrix
#'
#' Sandwich estimator: `V = N * bread * (shat/N) * bread` using the AC
#' (autocorrelation consistent) meat. The `.ac_meat()` returns
#' `shat/N` (per-observation omega), so we multiply by N to get the
#' full VCV. This is analogous to the HAC/HC pattern where
#' `V = bread * meat * bread * N/(N-dofminus)`.
#'
#' @param bread K x K bread matrix.
#' @param X_hat N x K projected regressor matrix (sorted by time).
#' @param resid N-vector of residuals (sorted by time).
#' @param time_index List from `.build_time_index()`.
#' @param kernel Canonical kernel name.
#' @param bw Numeric bandwidth.
#' @param N Integer: number of observations.
#' @param K Integer: number of regressors.
#' @param dofminus Integer: large-sample DoF adjustment.
#' @param sdofminus Integer: small-sample DoF adjustment.
#' @param weights Normalized weights (sorted) or NULL.
#' @param weight_type Character: weight type.
#' @return K x K variance-covariance matrix.
#' @keywords internal
.compute_ac_vcov <- function(bread, X_hat, resid, time_index, kernel, bw,
                             N, K, dofminus = 0L, sdofminus = 0L,
                             small = FALSE,
                             weights = NULL, weight_type = "aweight") {
  # Precompute Z'WZ (using X_hat as "basis" for the VCV computation)
  if (!is.null(weights)) {
    ZwZ <- crossprod(X_hat, weights * X_hat)
  } else {
    ZwZ <- crossprod(X_hat)
  }

  omega <- .ac_meat(X_hat, resid, time_index, kernel, bw,
                     N, dofminus, weights, weight_type, ZwZ)

  V <- bread %*% omega %*% bread
  V <- (V + t(V)) / 2
  V <- V * N
  # Small-sample correction (Stata ivreg2.ado line 1184):
  # V *= (N - dofminus) / (N - K - dofminus - sdofminus)
  if (small) {
    V <- V * (N - dofminus) / (N - K - dofminus - sdofminus)
  }
  colnames(V) <- rownames(V) <- colnames(bread)
  V
}


# --------------------------------------------------------------------------
# .cluster_kernel_meat
# --------------------------------------------------------------------------
#' Compute cluster+kernel meat matrix (Driscoll-Kraay / Thompson time dimension)
#'
#' Aggregates observation-level scores by time period, then accumulates
#' kernel-weighted cross-lag products of time-clustered scores.
#' This is the time-dimension component used by both Driscoll-Kraay (DK)
#' and Thompson (two-way cluster+kernel) VCE.
#'
#' Algorithm matches livreg2.do lines 444–535.
#'
#' @param basis N x P matrix (X_hat for VCV, Z for diagnostics). Must be
#'   in sorted time order (matching `time_index`).
#' @param resid N-vector of residuals, in sorted time order.
#' @param time_index List from `.build_time_index()`.
#' @param kernel Canonical kernel name.
#' @param bw Numeric bandwidth.
#' @param weights Normalized weights (sorted) or NULL.
#' @param weight_type Character: `"aweight"` or `"pweight"`.
#' @return P x P symmetric meat matrix (unscaled).
#' @keywords internal
.cluster_kernel_meat <- function(basis, resid, time_index, kernel, bw,
                                  weights = NULL, weight_type = "aweight") {
  # 1. Observation-level scores
  scores <- .cl_scores(basis, resid, weights)

  # 2. Aggregate scores by time period
  tvar_sorted <- time_index$tvar_sorted
  time_scores <- rowsum(scores, tvar_sorted, reorder = FALSE)

  # 3. Diagonal (tau = 0): crossprod of time-clustered scores
  shat <- crossprod(time_scores)

  # 4. Determine max lag TAU
  ktype <- .kernel_type(kernel)
  if (ktype == "spectral") {
    TAU <- as.integer(time_index$T_span / time_index$tdelta - 1)
  } else {
    TAU <- as.integer(floor(bw))
  }

  if (TAU < 1L) return((shat + t(shat)) / 2)

  # 5. Off-diagonal loop (tau = 1..TAU)
  # Get unique sorted time values (in the order rowsum used)
  unique_times <- as.numeric(rownames(time_scores))
  tdelta <- time_index$tdelta
  n_times <- length(unique_times)

  for (tau in seq_len(TAU)) {
    kw <- .kernel_weights(tau, bw, kernel)
    if (kw == 0) next

    lag_gap <- tau * tdelta

    # Find time periods that have both t and t - lag_gap present
    target_times <- unique_times - lag_gap
    lag_idx <- match(target_times, unique_times)
    valid <- which(!is.na(lag_idx))

    if (length(valid) == 0L) next

    # Cross-product of time-clustered scores at lags
    ghat <- crossprod(time_scores[valid, , drop = FALSE],
                      time_scores[lag_idx[valid], , drop = FALSE])

    shat <- shat + kw * (ghat + t(ghat))
  }

  (shat + t(shat)) / 2
}


# --------------------------------------------------------------------------
# .cluster_kernel_scores_meat
# --------------------------------------------------------------------------
#' Compute cluster+kernel meat from pre-formed score vectors
#'
#' Like `.cluster_kernel_meat()` but takes pre-formed N x P scores directly.
#' Used by the Kleibergen-Paap path where scores are Kronecker products.
#'
#' @param scores N x P score matrix, in sorted time order.
#' @param time_index List from `.build_time_index()`.
#' @param kernel Canonical kernel name.
#' @param bw Numeric bandwidth.
#' @return P x P symmetric meat matrix (unscaled).
#' @keywords internal
.cluster_kernel_scores_meat <- function(scores, time_index, kernel, bw) {
  tvar_sorted <- time_index$tvar_sorted
  time_scores <- rowsum(scores, tvar_sorted, reorder = FALSE)

  shat <- crossprod(time_scores)

  ktype <- .kernel_type(kernel)
  if (ktype == "spectral") {
    TAU <- as.integer(time_index$T_span / time_index$tdelta - 1)
  } else {
    TAU <- as.integer(floor(bw))
  }

  if (TAU < 1L) return((shat + t(shat)) / 2)

  unique_times <- as.numeric(rownames(time_scores))
  tdelta <- time_index$tdelta

  for (tau in seq_len(TAU)) {
    kw <- .kernel_weights(tau, bw, kernel)
    if (kw == 0) next

    lag_gap <- tau * tdelta
    target_times <- unique_times - lag_gap
    lag_idx <- match(target_times, unique_times)
    valid <- which(!is.na(lag_idx))

    if (length(valid) == 0L) next

    ghat <- crossprod(time_scores[valid, , drop = FALSE],
                      time_scores[lag_idx[valid], , drop = FALSE])

    shat <- shat + kw * (ghat + t(ghat))
  }

  (shat + t(shat)) / 2
}


# --------------------------------------------------------------------------
# .compute_cluster_kernel_vcov
# --------------------------------------------------------------------------
#' Compute cluster+kernel VCV (Driscoll-Kraay or Thompson)
#'
#' Sandwich estimator with three internal paths:
#' - DK (one-way cluster+kernel on tvar): meat = cluster_kernel_meat
#' - Thompson (two-way): meat = cluster_meat(ivar) + cluster_kernel_meat(tvar) - hc_meat
#' Small-sample correction: `(N-1)/(N-K-sdofminus) * M/(M-1)` when `small=TRUE`.
#'
#' @param bread K x K bread matrix.
#' @param X_hat N x K projected regressor matrix (sorted by time).
#' @param resid N-vector of residuals (sorted by time).
#' @param cluster_vec Cluster vector (DK: tvar vector; Thompson: list of 2).
#' @param time_index List from `.build_time_index()`.
#' @param kernel Canonical kernel name.
#' @param bw Numeric bandwidth.
#' @param N Integer: number of observations.
#' @param K Integer: number of regressors.
#' @param M Integer: effective cluster count (min(M1, M2) for Thompson).
#' @param small Logical: apply small-sample correction.
#' @param dofminus Integer: large-sample DoF adjustment.
#' @param sdofminus Integer: small-sample DoF adjustment.
#' @param weights Normalized weights (sorted) or NULL.
#' @param weight_type Character: weight type.
#' @param is_twoway Logical: TRUE for Thompson, FALSE for DK.
#' @return K x K variance-covariance matrix.
#' @keywords internal
.compute_cluster_kernel_vcov <- function(bread, X_hat, resid,
                                          cluster_vec, time_index,
                                          kernel, bw,
                                          N, K, M, small,
                                          dofminus = 0L, sdofminus = 0L,
                                          weights = NULL,
                                          weight_type = "aweight",
                                          is_twoway = FALSE) {
  scores <- .cl_scores(X_hat, resid, weights)

  if (is_twoway) {
    # Thompson: CGM decomposition with kernel-smoothed time dimension
    # meat = cluster_meat(ivar) + cluster_kernel_meat(tvar) - hac_meat
    # The third term is the HAC meat (observation-level with kernel lags),
    # NOT just HC meat. This matches Stata livreg2.do line 336 where
    # shat3 = shat * (N-dofminus) — the full HAC meat from the robust block.
    meat_ivar <- crossprod(rowsum(scores, cluster_vec[[1L]], reorder = FALSE))
    meat_ivar <- (meat_ivar + t(meat_ivar)) / 2
    meat_tvar <- .cluster_kernel_meat(X_hat, resid, time_index, kernel, bw,
                                      weights, weight_type)
    meat_hac <- .hac_scores_meat(scores, time_index, kernel, bw)
    meat <- meat_ivar + meat_tvar - meat_hac
  } else {
    # DK: one-way cluster+kernel on tvar
    meat <- .cluster_kernel_meat(X_hat, resid, time_index, kernel, bw,
                                 weights, weight_type)
  }

  V <- bread %*% meat %*% bread
  V <- (V + t(V)) / 2

  if (small) {
    V <- V * ((N - 1) / (N - K - sdofminus)) * (M / (M - 1))
  }

  colnames(V) <- rownames(V) <- colnames(bread)
  V
}
