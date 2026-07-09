# ==========================================================================
# Tests for M-28: Stored Scalar Results (absorbed into hf suite)
# ==========================================================================
# Verifies yy, yyc, rankxx, rankzz, condxx, condzz, ll, ccev, cdev against
# hf-suite fixtures on the canonical mroz bases H31 (2SLS, help.txt:1274) and
# H61 (OLS via ivreg2, help.txt:1416). Absorbed into the hf suite 2026-07-04
# per planning/22-spec-matrix.md; these quantities are VCE-invariant, so two
# cells (IV overid + OLS) suffice to cover the family. Originally Ticket R2.

# Helper: read a Stata fixture CSV into a named list (skip if absent)
read_sr_fixture <- function(name) {
  path <- fixture_path(paste0("hf_mroz_stored_", name, ".csv"))
  if (!file.exists(path)) skip(paste("Fixture not found:", path))
  read_scalar_fixture(path)
}

# === Data: mroz estimation sample (Stata's e(sample) = lwage non-missing) ===
data(mroz, package = "ivreg2r")
mroz_lw <- mroz[!is.na(mroz$lwage), ]


# ==========================================================================
# Stata parity: mroz H31 — 2SLS, overidentified (help.txt:1274; D5a)
# ==========================================================================
test_that("stored results match Stata: mroz H31 2SLS", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw)
  fx <- read_sr_fixture("H31")

  expect_equal(fit$yy, fx$yy, tolerance = stata_tol$coef)
  expect_equal(fit$yyc, fx$yyc, tolerance = stata_tol$coef)
  expect_identical(as.integer(fit$rank), as.integer(fx$rankxx))
  expect_identical(as.integer(fit$rankzz), as.integer(fx$rankzz))
  expect_equal(fit$condxx, fx$condxx, tolerance = stata_tol$coef)
  expect_equal(fit$condzz, fx$condzz, tolerance = stata_tol$coef)
  expect_equal(fit$ll, fx$ll, tolerance = stata_tol$coef)
  expect_identical(nobs(fit), as.integer(fx$N))
  expect_equal(fit$rss, fx$rss, tolerance = stata_tol$coef)

  # ccev/cdev (single endogenous → length-1 vectors)
  expect_equal(fit$diagnostics$ccev[1], fx$ccev1, tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$cdev[1], fx$cdev1, tolerance = stata_tol$stat)
})


# ==========================================================================
# Stata parity: mroz H61 — OLS via ivreg2, robust (help.txt:1416; D5a)
# ==========================================================================
test_that("stored results match Stata: mroz H61 OLS", {
  # `robust` is part of the verbatim H61 command; these scalars are
  # VCE-invariant so the choice of vcov does not affect this comparison.
  fit <- ivreg2(lwage ~ exper + expersq + age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust")
  fx <- read_sr_fixture("H61")

  expect_equal(fit$yy, fx$yy, tolerance = stata_tol$coef)
  expect_equal(fit$yyc, fx$yyc, tolerance = stata_tol$coef)
  expect_identical(as.integer(fit$rank), as.integer(fx$rankxx))
  expect_identical(as.integer(fit$rankzz), as.integer(fx$rankzz))
  expect_equal(fit$condxx, fx$condxx, tolerance = stata_tol$coef)
  expect_equal(fit$condzz, fx$condzz, tolerance = stata_tol$coef)
  expect_equal(fit$ll, fx$ll, tolerance = stata_tol$coef)
  expect_identical(nobs(fit), as.integer(fx$N))
  expect_equal(fit$rss, fx$rss, tolerance = stata_tol$coef)

  # OLS convention: Z = X, so rankzz == rankxx and condzz == condxx exactly.
  expect_identical(as.integer(fx$rankzz), as.integer(fx$rankxx))
  expect_equal(fx$condzz, fx$condxx)

  # OLS: no ccev/cdev in this fixture or in the fit
  expect_null(fit$diagnostics$ccev)
  expect_null(fit$diagnostics$cdev)
})


# ==========================================================================
# Cross-copy identity: yy/yyc against the direct weighted TSS formulas
# ==========================================================================
test_that("yy/yyc match the direct weighted TSS formulas across estimator paths", {
  # fit$yy (tss_u) and fit$yyc (tss_c) are computed by six copy-pasted
  # blocks: fit-ols.R:64-74, fit-2sls.R:94-104, fit-liml.R:182-192, and
  # fit-gmm.R ~140-145 (gmm2s), ~335-340 (wmatrix), ~619-625 (cue). The
  # Stata-parity fixtures above only exercise the fit-2sls (H31) and
  # fit-ols (H61) copies, so this identity test pins the tss_u/tss_c lines
  # in the copies the Stata fixtures do not reach (M-28 review follow-up).
  form <- lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6
  y <- mroz_lw$lwage
  n <- nrow(mroz_lw)

  fits <- list(
    "2SLS"        = ivreg2(form, data = mroz_lw),
    "LIML"        = ivreg2(form, data = mroz_lw, method = "liml"),
    "GMM2S"       = ivreg2(form, data = mroz_lw, method = "gmm2s"),
    "CUE"         = ivreg2(form, data = mroz_lw, method = "cue"),
    "H31 aweight" = ivreg2(form, data = mroz_lw, weights = age),
    "H31 fweight" = ivreg2(form, data = mroz_lw, weights = age,
                            weight_type = "fweight")
  )

  for (nm in names(fits)) {
    fit <- fits[[nm]]

    # Reconstruct the exact weight vector the tss_u/tss_c blocks consume
    # (`w <- parsed$weights`): .prepare_model() normalizes aweights to sum
    # to N and leaves fweights as-is (redefining N as sum(w)) -- see
    # ivreg2.R, "Validate and normalize weights". fit$weights stores the
    # RAW user weights, so aweight must be re-normalized here to match.
    if (is.null(fit$weights)) {
      w <- rep(1, n)
    } else if (fit$weight_type == "fweight") {
      w <- fit$weights
    } else {
      w <- fit$weights * (n / sum(fit$weights))
    }

    ref_yy <- sum(w * y^2)
    m <- sum(w * y) / sum(w)
    ref_yyc <- sum(w * (y - m)^2)

    expect_equal(fit$yy, ref_yy, tolerance = 1e-10, label = paste0(nm, " yy"))
    expect_equal(fit$yyc, ref_yyc, tolerance = 1e-10, label = paste0(nm, " yyc"))
  }
})


# ==========================================================================
# Edge cases
# ==========================================================================
test_that("ccev/cdev are NULL for OLS", {
  data(card)
  fit <- ivreg2(lwage ~ educ + exper, data = card)
  expect_null(fit$diagnostics$ccev)
  expect_null(fit$diagnostics$cdev)
})

test_that("ccev/cdev are NULL when noid = TRUE", {
  data(card)
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card, noid = TRUE)
  expect_null(fit$diagnostics$ccev)
  expect_null(fit$diagnostics$cdev)
})

test_that("new stored results columns are stored on the fitted object", {
  data(card)
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card)
  expect_false(is.null(fit$yy))
  expect_false(is.null(fit$yyc))
  expect_false(is.null(fit$rank))
  expect_false(is.null(fit$rankzz))
  expect_false(is.null(fit$condxx))
  expect_false(is.null(fit$condzz))
  expect_false(is.null(fit$ll))
  # Assert non-NULL explicitly: min(NULL) returns Inf (not NA), so an
  # is.na() check alone would pass silently if these ever regressed to NULL.
  expect_false(is.null(fit$diagnostics$ccev))
  expect_false(is.null(fit$diagnostics$cdev))
  expect_true(all(is.finite(fit$diagnostics$ccev)))
  expect_true(all(is.finite(fit$diagnostics$cdev)))
})

test_that("OLS has condzz = condxx and ccev/cdev are NULL", {
  data(card)
  fit <- ivreg2(lwage ~ educ + exper, data = card)
  expect_equal(fit$condzz, fit$condxx)
  expect_equal(fit$rankzz, fit$rank)
  expect_null(fit$diagnostics$ccev)
  expect_null(fit$diagnostics$cdev)
})
