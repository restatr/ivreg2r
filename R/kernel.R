# --------------------------------------------------------------------------
# Kernel weight functions for HAC/AC VCE
# --------------------------------------------------------------------------
# Implements the 8 kernel weight functions used by Stata's ivreg2 for
# HAC (heteroskedasticity and autocorrelation consistent) and AC
# (autocorrelation consistent) variance estimation.
#
# Reference: livreg2.do — s_vkernel() (lines 65–135), m_omega() (lines
# 161–167), m_calckw() (lines 627–664).
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# .validate_kernel
# --------------------------------------------------------------------------
#' Validate and canonicalize a kernel name
#'
#' Matches Stata's abbreviation table from `s_vkernel()` (livreg2.do:96–115),
#' case-insensitive. Accepts Stata's "Danielle" spelling but returns the
#' standard mathematical name "Daniell".
#'
#' @param kernel Character string: kernel name or abbreviation.
#' @return Canonical kernel name string.
#' @keywords internal
.validate_kernel <- function(kernel) {
  if (!is.character(kernel) || length(kernel) != 1L || is.na(kernel)) {
    stop("invalid kernel: must be a single character string", call. = FALSE)
  }
  # Abbreviation lookup table matching Stata's vklist (livreg2.do:96–115)
  kname <- trimws(tolower(kernel))

  # Empty string defaults to Bartlett (Stata's default kernel)
  if (kname == "") return("Bartlett")

  lookup <- c(
    "bar"                = "Bartlett",
    "bartlett"           = "Bartlett",
    "par"                = "Parzen",
    "parzen"             = "Parzen",
    "tru"                = "Truncated",
    "truncated"          = "Truncated",
    "thann"              = "Tukey-Hanning",
    "tukey-hanning"      = "Tukey-Hanning",
    "thamm"              = "Tukey-Hamming",
    "tukey-hamming"      = "Tukey-Hamming",
    "qua"                = "Quadratic Spectral",
    "qs"                 = "Quadratic Spectral",
    "quadratic-spectral" = "Quadratic Spectral",
    "quadratic spectral" = "Quadratic Spectral",
    "dan"                = "Daniell",
    "daniell"            = "Daniell",
    "danielle"           = "Daniell",
    "ten"                = "Tent",
    "tent"               = "Tent"
  )
  canonical <- lookup[kname]
  if (is.na(canonical)) {
    stop("invalid kernel: '", kernel, "'", call. = FALSE)
  }
  unname(canonical)
}


# --------------------------------------------------------------------------
# .kernel_type
# --------------------------------------------------------------------------
#' Classify kernel as lag-window or spectral-window
#'
#' Matches `m_omega()` classification at livreg2.do:161–167.
#'
#' @param kernel Canonical kernel name (from `.validate_kernel()`).
#' @return `"lag"` or `"spectral"`.
#' @keywords internal
.kernel_type <- function(kernel) {
  switch(kernel,
    "Bartlett"           = "lag",
    "Parzen"             = "lag",
    "Truncated"          = "lag",
    "Tukey-Hanning"      = "lag",
    "Tukey-Hamming"      = "lag",
    "Quadratic Spectral" = "spectral",
    "Daniell"            = "spectral",
    "Tent"               = "spectral",
    stop("unknown kernel: '", kernel, "'", call. = FALSE)
  )
}


# --------------------------------------------------------------------------
# .kernel_supports_auto_bw
# --------------------------------------------------------------------------
#' Check whether a kernel supports automatic bandwidth selection
#'
#' Only Bartlett, Parzen, and Quadratic Spectral support `bw = "auto"`.
#' Matches ivreg2.ado:4746.
#'
#' @param kernel Canonical kernel name.
#' @return Logical scalar.
#' @keywords internal
.kernel_supports_auto_bw <- function(kernel) {
  kernel %in% c("Bartlett", "Parzen", "Quadratic Spectral")
}


# --------------------------------------------------------------------------
# .validate_bandwidth
# --------------------------------------------------------------------------
#' Validate bandwidth for HAC estimation
#'
#' Rules from `s_vkernel()` at livreg2.do:71–91 and ivreg2.ado:4746.
#'
#' @param bw Numeric > 0, or the string `"auto"`.
#' @param kernel Canonical kernel name (from `.validate_kernel()`).
#' @return `bw` (invisibly), after validation. Issues warning for bw=1 with
#'   kernels where this means zero lags.
#' @keywords internal
.validate_bandwidth <- function(bw, kernel) {
  if (is.character(bw)) {
    if (length(bw) != 1L || tolower(bw) != "auto") {
      stop("bandwidth must be numeric > 0 or \"auto\"", call. = FALSE)
    }
    if (!.kernel_supports_auto_bw(kernel)) {
      stop("kernel '", kernel, "' is not compatible with bw(auto)", call. = FALSE)
    }
    return(invisible(bw))
  }
  if (!is.numeric(bw) || length(bw) != 1L || is.na(bw) || !is.finite(bw) || bw <= 0) {
    stop("bandwidth must be numeric > 0 or \"auto\"", call. = FALSE)
  }
  # Warn if bw=1 means zero lags for this kernel (flag "0" in Stata's vklist)
  if (bw == 1 && kernel %in% c("Bartlett", "Parzen",
                                 "Tukey-Hanning", "Tukey-Hamming")) {
    warning("kernel='", kernel, "' and bw=1 implies zero lags used. ",
            "Standard errors and test statistics are not ",
            "autocorrelation-consistent.", call. = FALSE)
  }
  invisible(bw)
}


# --------------------------------------------------------------------------
# .kernel_weights
# --------------------------------------------------------------------------
#' Compute kernel weights for HAC estimation
#'
#' Takes raw lag `tau` and bandwidth `bw`, normalizes internally as
#' `x = tau / bw`. Vectorized over `tau`.
#'
#' Formulas match Stata's `m_calckw()` (livreg2.do:627–664) with added
#' singularity guards for numerical stability:
#' - QS near zero: exponential interpolation (from sandwich::kweights)
#' - Daniell at zero: limit = 1
#' - Tent at zero: limit = tau^2 (and 0 when tau = 0)
#' - Lag-window kernels: clipped to 0 for x > 1
#'
#' @param tau Numeric vector of lags (non-negative integers in practice).
#' @param bw Numeric scalar bandwidth (> 0).
#' @param kernel Canonical kernel name (from `.validate_kernel()`).
#' @return Numeric vector of kernel weights, same length as `tau`.
#' @keywords internal
.kernel_weights <- function(tau, bw, kernel) {
  x <- tau / bw

  switch(kernel,
    "Truncated" = {
      # Truncated: kw = 1 for all x, clipped at x > 1
      w <- rep(1, length(x))
      w[x > 1] <- 0
      w
    },

    "Bartlett" = {
      w <- 1 - x
      w[x > 1] <- 0
      w
    },

    "Parzen" = {
      w <- ifelse(x <= 0.5,
                   1 - 6 * x^2 + 6 * x^3,
                   2 * (1 - x)^3)
      w[x > 1] <- 0
      w
    },

    "Tukey-Hanning" = {
      w <- 0.5 + 0.5 * cos(pi * x)
      w[x > 1] <- 0
      w
    },

    "Tukey-Hamming" = {
      w <- 0.54 + 0.46 * cos(pi * x)
      w[x > 1] <- 0
      w
    },

    "Quadratic Spectral" = {
      # QS: 25/(12*pi^2*x^2) * (sin(6*pi*x/5)/(6*pi*x/5) - cos(6*pi*x/5))
      # Singularity guard from sandwich::kweights(): exponential interpolation
      # near zero to avoid 0/0.
      y <- 6 * pi * x / 5
      w <- ifelse(x == 0, 1,
                   (25 / (12 * pi^2 * x^2)) * (sin(y) / y - cos(y)))
      # Exponential interpolation for |x| < 1e-3
      ix <- which(abs(x) > 0 & abs(x) < 1e-3)
      if (length(ix) > 0L) {
        # Compute QS at 1e-3 as anchor point
        y0 <- 6 * pi * 1e-3 / 5
        qs0 <- (25 / (12 * pi^2 * 1e-6)) * (sin(y0) / y0 - cos(y0))
        cf <- 1e6 * log(qs0)
        w[ix] <- exp(cf * x[ix]^2)
      }
      w
    },

    "Daniell" = {
      # sinc function: sin(pi*x) / (pi*x), limit = 1 at x = 0
      ifelse(abs(x) < 1e-10, 1, sin(pi * x) / (pi * x))
    },

    "Tent" = {
      # Stata: 2*(1 - cos(tau * karg)) / karg^2 where karg = tau/bw
      # Note: tau appears explicitly in the cosine argument, not just via x.
      karg <- x  # = tau / bw
      w <- ifelse(abs(karg) < 1e-10,
                   tau^2,
                   2 * (1 - cos(tau * karg)) / karg^2)
      # When tau = 0 and karg ~ 0, the limit is tau^2 = 0
      w[tau == 0] <- 0
      w
    },

    stop("unknown kernel: '", kernel, "'", call. = FALSE)
  )
}


# --------------------------------------------------------------------------
# .auto_bandwidth
# --------------------------------------------------------------------------
#' Automatic bandwidth selection via Newey-West (1994)
#'
#' Implements the Newey-West (1994, REStud 61(4):631-653) plug-in bandwidth
#' selector, matching Stata's `s_abw()` (ivreg2.ado:4796-4911). Supported
#' for Bartlett, Parzen, and Quadratic Spectral kernels only.
#'
#' @param resid N-vector of residuals (sorted by time).
#' @param Z N x L instrument matrix (sorted by time), including exogenous
#'   regressors and (if present) the intercept column.
#' @param time_index List from `.build_time_index()`.
#' @param kernel Canonical kernel name (already validated as auto-compatible).
#' @param has_intercept Logical: whether the last column of Z is an intercept.
#'   If TRUE, the intercept column is excluded from the score process
#'   (Stata zeros it out via the `h` vector).
#' @param N Integer: number of observations (used as autocovariance
#'   denominator, matching Stata's `nobs`).
#' @return Numeric scalar: selected bandwidth (>= 1). Integer for
#'   Bartlett/Parzen, possibly fractional for Quadratic Spectral.
#' @keywords internal
.auto_bandwidth <- function(resid, Z, time_index, kernel, has_intercept, N) {

  # --- Kernel-specific constants (Newey-West 1994, Table II) ---
  # cgamma for Bartlett corrected per Alistair Hall (see Stata source)
  params <- switch(kernel,
    "Bartlett"           = list(expo = 2 / 9,  q = 1L, cgamma = 1.1447),
    "Parzen"             = list(expo = 4 / 25, q = 2L, cgamma = 2.6614),
    "Quadratic Spectral" = list(expo = 2 / 25, q = 2L, cgamma = 1.3221),
    stop("kernel '", kernel, "' not compatible with bw(auto)", call. = FALSE)
  )
  expo   <- params$expo
  q      <- params$q
  cgamma <- params$cgamma

  # --- tobs: time span in index units (Stata's tobs = T_span / tdelta) ---
  tobs <- time_index$T_span / time_index$tdelta

  # --- mstar: pilot lag count ---
  mstar <- trunc(20 * (tobs / 100)^expo)
  if (mstar == 0L) {
    warning("Time span too short for automatic bandwidth selection; ",
            "returning bw = 1.", call. = FALSE)
    return(1)
  }

  # --- Build score process f = rowSums(u * Z_no_intercept) ---
  # Stata: h = J(nrows1,1,1) \ J(nrows2,1,0) — zeros out intercept column
  # Then f = (u :* Z) * h
  # In R, the intercept is typically column 1 ("(Intercept)"), not the last.
  if (has_intercept) {
    intercept_col <- which(colnames(Z) == "(Intercept)")
    Z_score <- Z[, -intercept_col, drop = FALSE]
  } else {
    Z_score <- Z
  }
  # f = (resid * Z_score) %*% 1 — collapse to scalar per obs
  f <- drop((resid * Z_score) %*% rep(1, ncol(Z_score)))

  nobs <- N

  # --- Autocovariances via .lag_pairs() ---
  sigmahat <- numeric(mstar + 1L)

  # j = 0: variance
  sigmahat[1L] <- sum(f * f) / nobs

  # j = 1..mstar
  for (j in seq_len(mstar)) {
    pairs <- .lag_pairs(time_index, j)
    if (nrow(pairs) == 0L) next
    sigmahat[j + 1L] <- sum(f[pairs[, 1L]] * f[pairs[, 2L]]) / nobs
  }

  # --- Weighted sums: shat_q and shat_0 ---
  shat_0 <- sigmahat[1L]
  shat_q <- 0
  for (j in seq_len(mstar)) {
    shat_q <- shat_q + 2 * sigmahat[j + 1L] * j^q
    shat_0 <- shat_0 + 2 * sigmahat[j + 1L]
  }

  # Guard: no autocorrelation detected
  if (abs(shat_0) < .Machine$double.eps) {
    warning("No autocorrelation detected in score process; ",
            "returning bw = 1.", call. = FALSE)
    return(1)
  }

  # --- Optimal bandwidth ---
  expon <- 1 / (2 * q + 1)
  gammahat <- cgamma * ((shat_q / shat_0)^2)^expon
  m <- gammahat * tobs^expon

  # Bartlett/Parzen: truncate to integer; QS: keep fractional
  if (kernel %in% c("Bartlett", "Parzen")) {
    optlag <- min(trunc(m), mstar)
  } else {
    # Quadratic Spectral
    optlag <- min(m, mstar)
  }

  # Return bw = optlag + 1 (Stata convention: bw = max_lag + 1)
  optlag + 1
}
