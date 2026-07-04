# ==========================================================================
# Tests for M-28: Stored Scalar Results (absorbed into hf suite)
# ==========================================================================
# Verifies yy, yyc, rankxx, rankzz, condxx, condzz, ll, ccev, cdev against
# hf-suite fixtures on the canonical mroz bases H31 (2SLS, help.txt:1274) and
# H61 (OLS via ivreg2, help.txt:1416). Absorbed into the hf suite 2026-07-04
# per planning/22-spec-matrix.md; these quantities are VCE-invariant, so two
# cells (IV overid + OLS) suffice to cover the family. Originally Ticket R2.

# Helper: read a Stata fixture CSV into a named list
read_sr_fixture <- function(name) {
  path <- fixture_path(paste0("hf_mroz_stored_", name, ".csv"))
  if (!file.exists(path)) skip(paste("Fixture not found:", path))
  df <- read.csv(path, stringsAsFactors = FALSE)
  vals <- as.list(df$value)
  names(vals) <- df$quantity
  vals
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

test_that("glance includes new stored results columns", {
  data(card)
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card)
  gl <- glance(fit, diagnostics = TRUE)
  expect_true("yy" %in% names(gl))
  expect_true("yyc" %in% names(gl))
  expect_true("rankxx" %in% names(gl))
  expect_true("rankzz" %in% names(gl))
  expect_true("condxx" %in% names(gl))
  expect_true("condzz" %in% names(gl))
  expect_true("ll" %in% names(gl))
  expect_true("ccev_min" %in% names(gl))
  expect_true("cdev_min" %in% names(gl))
  expect_false(is.na(gl$ccev_min))
  expect_false(is.na(gl$cdev_min))
})

test_that("glance OLS has condzz = condxx and ccev/cdev are NA", {
  data(card)
  fit <- ivreg2(lwage ~ educ + exper, data = card)
  gl <- glance(fit, diagnostics = TRUE)
  expect_equal(gl$condzz, gl$condxx)
  expect_equal(gl$rankzz, gl$rankxx)
  expect_true(is.na(gl$ccev_min))
  expect_true(is.na(gl$cdev_min))
})
