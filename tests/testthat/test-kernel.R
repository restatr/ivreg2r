# ===========================================================================
# Tests for kernel weight functions (kernel.R)
# ===========================================================================

# ---------------------------------------------------------------------------
# 1. Kernel name validation (.validate_kernel)
# ---------------------------------------------------------------------------

test_that(".validate_kernel maps all abbreviations correctly", {
  # Bartlett

expect_equal(.validate_kernel(""), "Bartlett")
  expect_equal(.validate_kernel("bar"), "Bartlett")
  expect_equal(.validate_kernel("bartlett"), "Bartlett")

  # Parzen
  expect_equal(.validate_kernel("par"), "Parzen")
  expect_equal(.validate_kernel("parzen"), "Parzen")

  # Truncated
  expect_equal(.validate_kernel("tru"), "Truncated")
  expect_equal(.validate_kernel("truncated"), "Truncated")

  # Tukey-Hanning
  expect_equal(.validate_kernel("thann"), "Tukey-Hanning")
  expect_equal(.validate_kernel("tukey-hanning"), "Tukey-Hanning")

  # Tukey-Hamming
  expect_equal(.validate_kernel("thamm"), "Tukey-Hamming")
  expect_equal(.validate_kernel("tukey-hamming"), "Tukey-Hamming")

  # Quadratic Spectral
  expect_equal(.validate_kernel("qua"), "Quadratic Spectral")
  expect_equal(.validate_kernel("qs"), "Quadratic Spectral")
  expect_equal(.validate_kernel("quadratic-spectral"), "Quadratic Spectral")
  expect_equal(.validate_kernel("quadratic spectral"), "Quadratic Spectral")

  # Daniell (accepts Stata's "Danielle" spelling)
  expect_equal(.validate_kernel("dan"), "Daniell")
  expect_equal(.validate_kernel("daniell"), "Daniell")
  expect_equal(.validate_kernel("danielle"), "Daniell")

  # Tent
  expect_equal(.validate_kernel("ten"), "Tent")
  expect_equal(.validate_kernel("tent"), "Tent")
})

test_that(".validate_kernel is case-insensitive", {
  expect_equal(.validate_kernel("BARTLETT"), "Bartlett")
  expect_equal(.validate_kernel("Parzen"), "Parzen")
  expect_equal(.validate_kernel("TUKEY-HANNING"), "Tukey-Hanning")
  expect_equal(.validate_kernel("QS"), "Quadratic Spectral")
  expect_equal(.validate_kernel("DAN"), "Daniell")
  expect_equal(.validate_kernel("Tent"), "Tent")
})

test_that(".validate_kernel handles leading/trailing whitespace", {
  expect_equal(.validate_kernel("  bartlett  "), "Bartlett")
  expect_equal(.validate_kernel(" qs "), "Quadratic Spectral")
})

test_that(".validate_kernel errors on invalid kernel names", {
  expect_error(.validate_kernel("invalid"), "invalid kernel")
  expect_error(.validate_kernel("gaussian"), "invalid kernel")
  expect_error(.validate_kernel("hamming"), "invalid kernel")
  expect_error(.validate_kernel("hanning"), "invalid kernel")
})


# ---------------------------------------------------------------------------
# 2. Kernel type classification (.kernel_type)
# ---------------------------------------------------------------------------

test_that(".kernel_type classifies lag-window kernels correctly", {
  expect_equal(.kernel_type("Bartlett"), "lag")
  expect_equal(.kernel_type("Parzen"), "lag")
  expect_equal(.kernel_type("Truncated"), "lag")
  expect_equal(.kernel_type("Tukey-Hanning"), "lag")
  expect_equal(.kernel_type("Tukey-Hamming"), "lag")
})

test_that(".kernel_type classifies spectral-window kernels correctly", {
  expect_equal(.kernel_type("Quadratic Spectral"), "spectral")
  expect_equal(.kernel_type("Daniell"), "spectral")
  expect_equal(.kernel_type("Tent"), "spectral")
})

test_that(".kernel_type errors on unknown kernel", {
  expect_error(.kernel_type("Unknown"), "unknown kernel")
})


# ---------------------------------------------------------------------------
# 3. .kernel_supports_auto_bw
# ---------------------------------------------------------------------------

test_that(".kernel_supports_auto_bw returns TRUE for supported kernels", {
  expect_true(.kernel_supports_auto_bw("Bartlett"))
  expect_true(.kernel_supports_auto_bw("Parzen"))
  expect_true(.kernel_supports_auto_bw("Quadratic Spectral"))
})

test_that(".kernel_supports_auto_bw returns FALSE for unsupported kernels", {
  expect_false(.kernel_supports_auto_bw("Truncated"))
  expect_false(.kernel_supports_auto_bw("Tukey-Hanning"))
  expect_false(.kernel_supports_auto_bw("Tukey-Hamming"))
  expect_false(.kernel_supports_auto_bw("Daniell"))
  expect_false(.kernel_supports_auto_bw("Tent"))
})


# ---------------------------------------------------------------------------
# 4. Bandwidth validation (.validate_bandwidth)
# ---------------------------------------------------------------------------

test_that(".validate_bandwidth accepts valid numeric bw", {
  expect_silent(.validate_bandwidth(5, "Bartlett"))
  expect_silent(.validate_bandwidth(0.5, "Parzen"))
  expect_silent(.validate_bandwidth(100, "Quadratic Spectral"))
  expect_silent(.validate_bandwidth(2.7, "Daniell"))
})

test_that(".validate_bandwidth errors on non-positive bw", {
  expect_error(.validate_bandwidth(0, "Bartlett"), "numeric > 0")
  expect_error(.validate_bandwidth(-1, "Bartlett"), "numeric > 0")
  expect_error(.validate_bandwidth(-0.5, "Parzen"), "numeric > 0")
})

test_that(".validate_bandwidth errors on non-numeric bw", {
  expect_error(.validate_bandwidth(NA, "Bartlett"), "numeric > 0")
  expect_error(.validate_bandwidth(NaN, "Bartlett"), "numeric > 0")
  expect_error(.validate_bandwidth("foo", "Bartlett"), "numeric > 0")
})

test_that(".validate_bandwidth accepts bw='auto' for compatible kernels", {
  expect_silent(.validate_bandwidth("auto", "Bartlett"))
  expect_silent(.validate_bandwidth("auto", "Parzen"))
  expect_silent(.validate_bandwidth("auto", "Quadratic Spectral"))
})

test_that(".validate_bandwidth errors on bw='auto' for incompatible kernels", {
  expect_error(.validate_bandwidth("auto", "Truncated"), "not compatible")
  expect_error(.validate_bandwidth("auto", "Tukey-Hanning"), "not compatible")
  expect_error(.validate_bandwidth("auto", "Tukey-Hamming"), "not compatible")
  expect_error(.validate_bandwidth("auto", "Daniell"), "not compatible")
  expect_error(.validate_bandwidth("auto", "Tent"), "not compatible")
})

test_that(".validate_bandwidth warns when bw=1 means zero lags", {
  expect_warning(.validate_bandwidth(1, "Bartlett"), "zero lags")
  expect_warning(.validate_bandwidth(1, "Parzen"), "zero lags")
  expect_warning(.validate_bandwidth(1, "Tukey-Hanning"), "zero lags")
  expect_warning(.validate_bandwidth(1, "Tukey-Hamming"), "zero lags")
})

test_that(".validate_bandwidth does NOT warn for bw=1 with non-flagged kernels", {
  expect_silent(.validate_bandwidth(1, "Truncated"))
  expect_silent(.validate_bandwidth(1, "Quadratic Spectral"))
  expect_silent(.validate_bandwidth(1, "Daniell"))
  expect_silent(.validate_bandwidth(1, "Tent"))
})


# ---------------------------------------------------------------------------
# 5. Kernel weight values (.kernel_weights)
# ---------------------------------------------------------------------------

# --- 5a. Boundary values at tau = 0 (all kernels should return 1 or 0) ---

test_that("lag-window kernels return 1 at tau = 0", {
  bw <- 5
  for (k in c("Bartlett", "Parzen", "Truncated",
              "Tukey-Hanning", "Tukey-Hamming")) {
    expect_equal(.kernel_weights(0, bw, k), 1,
                 label = paste(k, "at tau=0"))
  }
})

test_that("spectral kernels return correct limit at tau = 0", {
  bw <- 5
  expect_equal(.kernel_weights(0, bw, "Quadratic Spectral"), 1)
  expect_equal(.kernel_weights(0, bw, "Daniell"), 1)
  expect_equal(.kernel_weights(0, bw, "Tent"), 0)
})

# --- 5b. Lag-window kernels return 0 at x = 1 or beyond ---

test_that("Bartlett returns 0 at tau = bw and beyond", {
  expect_equal(.kernel_weights(5, 5, "Bartlett"), 0)
  expect_equal(.kernel_weights(6, 5, "Bartlett"), 0)
  expect_equal(.kernel_weights(10, 5, "Bartlett"), 0)
})

test_that("Parzen returns 0 at tau = bw and beyond", {
  expect_equal(.kernel_weights(5, 5, "Parzen"), 0)
  expect_equal(.kernel_weights(6, 5, "Parzen"), 0)
})

test_that("Tukey-Hanning returns 0 at tau = bw and beyond", {
  expect_equal(.kernel_weights(5, 5, "Tukey-Hanning"), 0,
               tolerance = 1e-15)
  expect_equal(.kernel_weights(6, 5, "Tukey-Hanning"), 0)
})

test_that("Tukey-Hamming does NOT return 0 at x = 1", {
  # Tukey-Hamming: 0.54 + 0.46*cos(pi) = 0.54 - 0.46 = 0.08
  # But we clip at x > 1, so at exactly x = 1 it's 0.08
  expect_equal(.kernel_weights(5, 5, "Tukey-Hamming"), 0.08,
               tolerance = 1e-14)
  # Beyond x = 1, clipped to 0
  expect_equal(.kernel_weights(6, 5, "Tukey-Hamming"), 0)
})

test_that("Truncated returns 1 for tau < bw and 0 for tau > bw", {
  expect_equal(.kernel_weights(3, 5, "Truncated"), 1)
  expect_equal(.kernel_weights(5, 5, "Truncated"), 1)  # x = 1 exactly
  expect_equal(.kernel_weights(6, 5, "Truncated"), 0)
})

# --- 5c. Known interior values ---

test_that("Bartlett at x = 0.5 gives 0.5", {
  expect_equal(.kernel_weights(2.5, 5, "Bartlett"), 0.5)
})

test_that("Parzen at x = 0.5 gives 0.25 (both branches agree)", {
  expect_equal(.kernel_weights(2.5, 5, "Parzen"), 0.25)
})

test_that("Tukey-Hanning at x = 0.5 gives 0.5", {
  expect_equal(.kernel_weights(2.5, 5, "Tukey-Hanning"), 0.5)
})

test_that("Tukey-Hamming at x = 0.5 gives 0.54", {
  expect_equal(.kernel_weights(2.5, 5, "Tukey-Hamming"), 0.54,
               tolerance = 1e-14)
})

# --- 5d. Cross-check with sandwich::kweights (5 shared kernels) ---

test_that("kernel weights match sandwich::kweights for shared kernels", {
  skip_if_not_installed("sandwich")

  shared_kernels <- c("Truncated", "Bartlett", "Parzen",
                      "Tukey-Hanning", "Quadratic Spectral")
  bw <- 10
  tau_vals <- c(1, 2, 3, 4, 5, 7, 9, 10, 12, 15, 20)

  for (k in shared_kernels) {
    x_vals <- tau_vals / bw
    expected <- sandwich::kweights(x_vals, kernel = k)
    actual <- .kernel_weights(tau_vals, bw, k)
    expect_equal(actual, expected, tolerance = 1e-12,
                 label = paste(k, "vs sandwich"))
  }
})

test_that("QS matches sandwich at many x values", {
  skip_if_not_installed("sandwich")
  bw <- 10
  tau_vals <- seq(0, 30, by = 0.5)
  x_vals <- tau_vals / bw
  expected <- sandwich::kweights(x_vals, kernel = "Quadratic Spectral")
  actual <- .kernel_weights(tau_vals, bw, "Quadratic Spectral")
  expect_equal(actual, expected, tolerance = 1e-10,
               label = "QS fine grid vs sandwich")
})

# --- 5e. Hand-computed reference values for Tukey-Hamming, Daniell, Tent ---

test_that("Tukey-Hamming matches hand-computed values", {
  bw <- 10
  tau <- c(0, 1, 2.5, 5, 7.5, 10)
  x <- tau / bw
  expected <- 0.54 + 0.46 * cos(pi * x)
  expected[x > 1] <- 0
  actual <- .kernel_weights(tau, bw, "Tukey-Hamming")
  expect_equal(actual, expected, tolerance = 1e-14)
})

test_that("Daniell matches hand-computed sinc values", {
  bw <- 10
  tau <- c(1, 2.5, 5, 7.5, 10, 15, 20)
  x <- tau / bw
  expected <- sin(pi * x) / (pi * x)
  actual <- .kernel_weights(tau, bw, "Daniell")
  expect_equal(actual, expected, tolerance = 1e-14)
})

test_that("Daniell produces negative weights where expected", {
  # sinc(x) < 0 for x in (1, 2), so tau/bw in (1, 2) should be negative
  bw <- 10
  tau_neg <- c(11, 12, 15, 18, 19)
  w <- .kernel_weights(tau_neg, bw, "Daniell")
  expect_true(all(w < 0), label = "Daniell negative for x in (1,2)")
})

test_that("Tent matches hand-computed values", {
  # Stata: kw = 2*(1-cos(tau*karg))/karg^2, karg = tau/bw
  tau <- c(2, 3, 1)
  bw <- c(5, 4, 10)
  expected <- c(
    2 * (1 - cos(2 * (2 / 5))) / (2 / 5)^2,  # tau=2, bw=5
    2 * (1 - cos(3 * (3 / 4))) / (3 / 4)^2,  # tau=3, bw=4
    2 * (1 - cos(1 * (1 / 10))) / (1 / 10)^2  # tau=1, bw=10
  )
  for (i in seq_along(tau)) {
    actual <- .kernel_weights(tau[i], bw[i], "Tent")
    expect_equal(actual, expected[i], tolerance = 1e-12,
                 label = sprintf("Tent tau=%d, bw=%d", tau[i], bw[i]))
  }
})


# ---------------------------------------------------------------------------
# 6. Singularity handling
# ---------------------------------------------------------------------------

test_that("QS at x = 0 returns 1, no NaN", {
  w <- .kernel_weights(0, 10, "Quadratic Spectral")
  expect_false(is.nan(w))
  expect_equal(w, 1)
})

test_that("QS at very small x returns value near 1, no NaN", {
  # x = 1e-6 (tau = 1e-5, bw = 10)
  w <- .kernel_weights(1e-5, 10, "Quadratic Spectral")
  expect_false(is.nan(w))
  expect_true(w > 0.999)
  expect_true(w <= 1)
})

test_that("Daniell at x = 0 returns 1, no NaN", {
  w <- .kernel_weights(0, 10, "Daniell")
  expect_false(is.nan(w))
  expect_equal(w, 1)
})

test_that("Tent at tau = 0, any bw returns 0, no NaN", {
  w <- .kernel_weights(0, 10, "Tent")
  expect_false(is.nan(w))
  expect_equal(w, 0)

  w2 <- .kernel_weights(0, 0.5, "Tent")
  expect_false(is.nan(w2))
  expect_equal(w2, 0)
})


# ---------------------------------------------------------------------------
# 7. Vectorized evaluation
# ---------------------------------------------------------------------------

test_that("kernel_weights is vectorized over tau", {
  bw <- 10
  tau_vec <- 0:15

  # Bartlett
  w_bart <- .kernel_weights(tau_vec, bw, "Bartlett")
  expect_length(w_bart, length(tau_vec))
  expect_equal(w_bart[1], 1)           # tau = 0
  expect_equal(w_bart[6], 0.5)         # tau = 5, x = 0.5
  expect_equal(w_bart[11], 0)          # tau = 10, x = 1
  expect_true(all(w_bart[12:16] == 0)) # tau > bw

  # QS
  w_qs <- .kernel_weights(tau_vec, bw, "Quadratic Spectral")
  expect_length(w_qs, length(tau_vec))
  expect_false(any(is.nan(w_qs)))
  expect_equal(w_qs[1], 1)  # tau = 0
})

test_that("mixed vector with zeros and non-zeros works for all kernels", {
  bw <- 5
  tau_vec <- c(0, 1, 2, 3, 4, 5, 6, 7)
  all_kernels <- c("Bartlett", "Parzen", "Truncated",
                   "Tukey-Hanning", "Tukey-Hamming",
                   "Quadratic Spectral", "Daniell", "Tent")

  for (k in all_kernels) {
    w <- .kernel_weights(tau_vec, bw, k)
    expect_length(w, length(tau_vec))
    expect_false(any(is.nan(w)), label = paste(k, "no NaN"))
    expect_false(any(is.infinite(w)), label = paste(k, "no Inf"))
  }
})


# ---------------------------------------------------------------------------
# 8. Lag-window clipping
# ---------------------------------------------------------------------------

test_that("all lag-window kernels return 0 for tau > bw", {
  bw <- 5
  tau_beyond <- c(6, 7, 10, 100)
  lag_kernels <- c("Bartlett", "Parzen", "Truncated",
                   "Tukey-Hanning", "Tukey-Hamming")

  for (k in lag_kernels) {
    w <- .kernel_weights(tau_beyond, bw, k)
    expect_true(all(w == 0), label = paste(k, "clipped beyond bw"))
  }
})

test_that("spectral kernels are NOT clipped beyond x = 1", {
  bw <- 5
  # QS, Daniell, Tent should still return non-zero values beyond x = 1
  w_qs <- .kernel_weights(6, bw, "Quadratic Spectral")
  expect_true(w_qs != 0, label = "QS not clipped")

  # Daniell at x = 1.2 => sin(1.2*pi)/(1.2*pi) which is non-zero
  w_dan <- .kernel_weights(6, bw, "Daniell")
  expect_true(w_dan != 0, label = "Daniell not clipped")

  # Tent: 2*(1-cos(6 * 6/5)) / (6/5)^2 which is non-zero
  w_tent <- .kernel_weights(6, bw, "Tent")
  expect_true(w_tent != 0, label = "Tent not clipped")
})
