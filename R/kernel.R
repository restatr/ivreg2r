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
