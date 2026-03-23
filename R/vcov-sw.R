# --------------------------------------------------------------------------
# Stock-Watson (2008) panel-robust VCE (Ticket Q3)
# --------------------------------------------------------------------------
# Heteroskedasticity and autocorrelation robust VCE for fixed-T panels.
# Ported from Stata's livreg2.do lines 548-599.


# --------------------------------------------------------------------------
# .sw_meat
# --------------------------------------------------------------------------
#' Compute Stock-Watson (2008) panel-robust meat matrix
#'
#' Implements the bias-corrected per-panel estimator from Stock & Watson
#' (2008, Econometrica, eq. 6). Panels with T <= 2 are skipped (SW is
#' undefined). The meat has its own internal normalization: divided by
#' `(N - N_panels)`, NOT by `N` or `(N - dofminus)`.
#'
#' @param basis N x P matrix (Z for diagnostics, X_hat for VCV).
#' @param resid N-vector of residuals.
#' @param ivar_vec N-vector of panel identifiers.
#' @param N Integer: number of observations.
#' @param weights Normalized weights or NULL.
#' @param weight_type Character: `"aweight"` or `"pweight"`.
#' @param center Logical: if `TRUE`, subtract the global mean of moment
#'   conditions before computing the outer product.
#' @param eZmean P-vector: pre-computed global mean of scores (used when
#'   `center = TRUE`), or NULL.
#' @return P x P symmetric meat matrix (already normalized by
#'   `1/(N - N_panels)`), or a zero matrix with a warning if all panels
#'   have T <= 2.
#' @keywords internal
.sw_meat <- function(basis, resid, ivar_vec, N, weights = NULL,
                     weight_type = "aweight", center = FALSE,
                     eZmean = NULL) {
  P <- ncol(basis)
  panels <- split(seq_len(N), ivar_vec)

  shat <- matrix(0, P, P)
  bhat <- matrix(0, P, P)
  N_panels <- 0L

  for (panel_rows in panels) {
    Ti <- length(panel_rows)
    if (Ti <= 2L) next
    N_panels <- N_panels + 1L

    ei <- resid[panel_rows]
    Zi <- basis[panel_rows, , drop = FALSE]

    if (is.null(weights)) {
      # Unweighted
      eZ <- ei * Zi  # Ti x P score matrix
      if (center && !is.null(eZmean)) {
        eZ <- sweep(eZ, 2L, eZmean)
      }
      shatsub <- crossprod(eZ)                     # P x P
      sigmahatsub <- sum(ei^2)                      # scalar
      ZZsub <- crossprod(Zi)                        # P x P
    } else {
      wi <- weights[panel_rows]
      # aweight/pweight: quadratic in weights (Stata: (w*wf)^2)
      eZ <- wi * ei * Zi  # Ti x P weighted score matrix
      if (center && !is.null(eZmean)) {
        eZ <- sweep(eZ, 2L, eZmean)
      }
      shatsub <- crossprod(eZ)                      # sum w^2 (eZ)(eZ)'
      sigmahatsub <- sum(wi^2 * ei^2)               # sum w^2 e^2
      ZZsub <- crossprod(wi * Zi)                   # sum w^2 Z'Z
    }

    # Bias correction (SW eq. 6)
    shat <- shat + shatsub * (Ti - 1) / (Ti - 2)
    bhat <- bhat + (ZZsub / Ti) * sigmahatsub / (Ti - 1) / (Ti - 2)
  }

  if (N_panels == 0L) {
    warning("Stock-Watson VCE: all panels have T <= 2; ",
            "SW adjustment is undefined.", call. = FALSE)
    return(matrix(0, P, P))
  }

  # SW normalization (livreg2.do lines 590-598)
  # "We use only an N-n degrees of freedom correction, ignoring k regressors"
  shat <- shat / (N - N_panels)
  bhat <- bhat / N_panels
  shat <- shat - bhat

  # Force symmetry
  (shat + t(shat)) / 2
}


# --------------------------------------------------------------------------
# .sw_cross_meat
# --------------------------------------------------------------------------
#' Compute cross-equation Stock-Watson meat block
#'
#' Used by the reduced-form system VCV to compute the (i,j) block of the
#' stacked meat matrix. Same panel-loop algorithm as `.sw_meat()` but with
#' cross-equation residuals e_i, e_j.
#'
#' @param basis N x P matrix.
#' @param resid_i,resid_j N-vectors of residuals from equations i and j.
#' @param panels List of index vectors (output of `split(seq_len(N), ivar_vec)`).
#' @param N Integer.
#' @param weights Normalized weights or NULL.
#' @param weight_type Character.
#' @param center Logical.
#' @param eZmean_i,eZmean_j P-vectors or NULL.
#' @return P x P matrix (already normalized by `1/(N - N_panels)`).
#' @keywords internal
.sw_cross_meat <- function(basis, resid_i, resid_j, panels, N,
                            weights = NULL, weight_type = "aweight",
                            center = FALSE, eZmean_i = NULL,
                            eZmean_j = NULL) {
  P <- ncol(basis)
  shat <- matrix(0, P, P)
  bhat <- matrix(0, P, P)
  N_panels <- 0L

  for (panel_rows in panels) {
    Ti <- length(panel_rows)
    if (Ti <= 2L) next
    N_panels <- N_panels + 1L

    ei <- resid_i[panel_rows]
    ej <- resid_j[panel_rows]
    Zi <- basis[panel_rows, , drop = FALSE]

    if (is.null(weights)) {
      eZi <- ei * Zi
      eZj <- ej * Zi
      if (center) {
        if (!is.null(eZmean_i)) eZi <- sweep(eZi, 2L, eZmean_i)
        if (!is.null(eZmean_j)) eZj <- sweep(eZj, 2L, eZmean_j)
      }
      shatsub <- crossprod(eZi, eZj)
      sigmahatsub <- sum(ei * ej)
      ZZsub <- crossprod(Zi)
    } else {
      wi <- weights[panel_rows]
      eZi <- wi * ei * Zi
      eZj <- wi * ej * Zi
      if (center) {
        if (!is.null(eZmean_i)) eZi <- sweep(eZi, 2L, eZmean_i)
        if (!is.null(eZmean_j)) eZj <- sweep(eZj, 2L, eZmean_j)
      }
      shatsub <- crossprod(eZi, eZj)
      sigmahatsub <- sum(wi^2 * ei * ej)
      ZZsub <- crossprod(wi * Zi)
    }

    shat <- shat + shatsub * (Ti - 1) / (Ti - 2)
    bhat <- bhat + (ZZsub / Ti) * sigmahatsub / (Ti - 1) / (Ti - 2)
  }

  if (N_panels == 0L) return(matrix(0, P, P))

  shat <- shat / (N - N_panels)
  bhat <- bhat / N_panels
  shat <- shat - bhat
  (shat + t(shat)) / 2
}


# --------------------------------------------------------------------------
# .sw_eZmean
# --------------------------------------------------------------------------
#' Compute global mean of scores for SW centering
#'
#' Helper that computes the centering term for SW meat when `center = TRUE`.
#' Matches the centering logic in `.hc_meat()` for the corresponding
#' weight type.
#'
#' @param basis N x P matrix.
#' @param resid N-vector.
#' @param N Integer.
#' @param weights Normalized weights or NULL.
#' @param weight_type Character.
#' @return P-vector of score means.
#' @keywords internal
.sw_eZmean <- function(basis, resid, N, weights = NULL,
                       weight_type = "aweight") {
  g <- basis * resid  # N x P unweighted scores
  if (is.null(weights)) {
    colMeans(g)
  } else {
    # aweight/pweight: Stata uses wv = (w*wf)^2
    colSums(weights^2 * g) / N
  }
}


# --------------------------------------------------------------------------
# .compute_sw_vcov
# --------------------------------------------------------------------------
#' Compute Stock-Watson (2008) panel-robust VCV
#'
#' Sandwich estimator using `.sw_meat()`. The SW meat has its own internal
#' normalization (`/(N - N_panels)`), so the sandwich multiplies by `N`
#' (not `N/(N-dofminus)` as for HC). Small-sample correction is the
#' standard robust correction `(N-dofminus)/(N-K-dofminus-sdofminus)`.
#'
#' @param bread K x K bread matrix.
#' @param X_hat N x K projected regressors (or X for OLS).
#' @param resid N-vector of residuals.
#' @param ivar_vec N-vector of panel identifiers.
#' @param N,K Integer dimensions.
#' @param small Logical: apply finite-sample correction.
#' @param dofminus,sdofminus Integer DoF adjustments.
#' @param weights Normalized weights or NULL.
#' @param weight_type Character.
#' @param center Logical.
#' @param psd Character or NULL: PSD correction mode.
#' @return K x K variance-covariance matrix.
#' @keywords internal
.compute_sw_vcov <- function(bread, X_hat, resid, ivar_vec, N, K,
                              small = FALSE, dofminus = 0L, sdofminus = 0L,
                              weights = NULL, weight_type = "aweight",
                              center = FALSE, psd = NULL) {
  eZmean <- if (center) {
    .sw_eZmean(X_hat, resid, N, weights, weight_type)
  } else {
    NULL
  }

  sw_shat <- .sw_meat(X_hat, resid, ivar_vec, N, weights, weight_type,
                       center = center, eZmean = eZmean)

  # PSD correction on the meat, matching Stata (livreg2.do:607-617).
  # This must happen here (not on the final VCV) because
  # B * psd_correct(S) * B != psd_correct(B * S * B) in general.
  # The caller (.compute_vcov) skips its own final-VCV PSD pass for SW.
  sw_shat <- .psd_correct(sw_shat, psd)

  # Sandwich: V = bread * sw_shat * bread * N
  # SW meat already normalized by /(N - N_panels), and Stata's assembly
  # multiplies by N (livreg2.do line 162: shat goes into N * bread * shat * bread)
  V <- bread %*% sw_shat %*% bread * N
  V <- (V + t(V)) / 2  # enforce symmetry

  # Small-sample correction (standard robust, ivreg2.ado line 1184)
  if (small) {
    V <- V * ((N - dofminus) / (N - K - dofminus - sdofminus))
  }

  colnames(V) <- rownames(V) <- colnames(bread)
  V
}
