# ==========================================================================
# Tests for Ticket R2: Additional Stored Results
# ==========================================================================
# Verifies yy, yyc, rankxx, rankzz, condxx, condzz, ll, ccev, cdev
# against Stata fixtures.

# Helper: read a Stata fixture CSV into a named list
read_sr_fixture <- function(name) {
  path <- fixture_path(paste0("stored_results_", name, ".csv"))
  if (!file.exists(path)) skip(paste("Fixture not found:", path))
  df <- read.csv(path, stringsAsFactors = FALSE)
  vals <- as.list(df$value)
  names(vals) <- df$quantity
  vals
}


# ==========================================================================
# 2SLS just-identified (IID)
# ==========================================================================
test_that("stored results match Stata: 2SLS just-identified", {
  data(card)
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  fx <- read_sr_fixture("justid")

  expect_equal(fit$yy, fx$yy, tolerance = stata_tol$coef)
  expect_equal(fit$yyc, fx$yyc, tolerance = stata_tol$coef)
  expect_identical(fit$rank, as.integer(fx$rankxx))
  expect_identical(fit$rankzz, as.integer(fx$rankzz))
  expect_equal(fit$condxx, fx$condxx, tolerance = stata_tol$coef)
  expect_equal(fit$condzz, fx$condzz, tolerance = stata_tol$coef)
  expect_equal(fit$ll, fx$ll, tolerance = stata_tol$coef)

  # ccev/cdev (single endogenous → length-1 vectors)
  expect_equal(fit$diagnostics$ccev[1], fx$ccev1, tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$cdev[1], fx$cdev1, tolerance = stata_tol$stat)
})


# ==========================================================================
# 2SLS overidentified (robust)
# ==========================================================================
test_that("stored results match Stata: 2SLS overidentified robust", {
  data(card)
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card, vcov = "robust")
  fx <- read_sr_fixture("overid_robust")

  expect_equal(fit$yy, fx$yy, tolerance = stata_tol$coef)
  expect_equal(fit$yyc, fx$yyc, tolerance = stata_tol$coef)
  expect_identical(fit$rank, as.integer(fx$rankxx))
  expect_identical(fit$rankzz, as.integer(fx$rankzz))
  expect_equal(fit$condxx, fx$condxx, tolerance = stata_tol$coef)
  expect_equal(fit$condzz, fx$condzz, tolerance = stata_tol$coef)
  expect_equal(fit$ll, fx$ll, tolerance = stata_tol$coef)
  expect_equal(fit$diagnostics$ccev[1], fx$ccev1, tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$cdev[1], fx$cdev1, tolerance = stata_tol$stat)
})


# ==========================================================================
# 2SLS overidentified with clustering
# ==========================================================================
test_that("stored results match Stata: 2SLS overid cluster", {
  data(card)
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card, clusters = ~smsa)
  fx <- read_sr_fixture("overid_cluster")

  expect_equal(fit$yy, fx$yy, tolerance = stata_tol$coef)
  expect_equal(fit$yyc, fx$yyc, tolerance = stata_tol$coef)
  expect_identical(fit$rank, as.integer(fx$rankxx))
  expect_identical(fit$rankzz, as.integer(fx$rankzz))
  expect_equal(fit$condxx, fx$condxx, tolerance = stata_tol$coef)
  expect_equal(fit$condzz, fx$condzz, tolerance = stata_tol$coef)
  expect_equal(fit$ll, fx$ll, tolerance = stata_tol$coef)
  expect_equal(fit$diagnostics$ccev[1], fx$ccev1, tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$cdev[1], fx$cdev1, tolerance = stata_tol$stat)
})


# ==========================================================================
# Aweight
# ==========================================================================
test_that("stored results match Stata: aweight", {
  data(card)
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card, weights = age)
  fx <- read_sr_fixture("aweight")

  expect_equal(fit$yy, fx$yy, tolerance = stata_tol$coef)
  expect_equal(fit$yyc, fx$yyc, tolerance = stata_tol$coef)
  expect_identical(fit$rank, as.integer(fx$rankxx))
  expect_identical(fit$rankzz, as.integer(fx$rankzz))
  expect_equal(fit$condxx, fx$condxx, tolerance = stata_tol$coef)
  expect_equal(fit$condzz, fx$condzz, tolerance = stata_tol$coef)
  expect_equal(fit$ll, fx$ll, tolerance = stata_tol$coef)
  expect_equal(fit$diagnostics$ccev[1], fx$ccev1, tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$cdev[1], fx$cdev1, tolerance = stata_tol$stat)
})


# ==========================================================================
# LIML
# ==========================================================================
test_that("stored results match Stata: LIML", {
  data(card)
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card, method = "liml")
  fx <- read_sr_fixture("liml")

  expect_equal(fit$yy, fx$yy, tolerance = stata_tol$coef)
  expect_equal(fit$yyc, fx$yyc, tolerance = stata_tol$coef)
  expect_identical(fit$rank, as.integer(fx$rankxx))
  expect_identical(fit$rankzz, as.integer(fx$rankzz))
  expect_equal(fit$condxx, fx$condxx, tolerance = stata_tol$coef)
  expect_equal(fit$condzz, fx$condzz, tolerance = stata_tol$coef)
  expect_equal(fit$ll, fx$ll, tolerance = stata_tol$coef)
  expect_equal(fit$diagnostics$ccev[1], fx$ccev1, tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$cdev[1], fx$cdev1, tolerance = stata_tol$stat)
})


# ==========================================================================
# Fweight
# ==========================================================================
test_that("stored results match Stata: fweight", {
  data(card)
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card, weights = age, weight_type = "fweight")
  fx <- read_sr_fixture("fweight")

  expect_equal(fit$yy, fx$yy, tolerance = stata_tol$coef)
  expect_equal(fit$yyc, fx$yyc, tolerance = stata_tol$coef)
  expect_identical(fit$rank, as.integer(fx$rankxx))
  expect_identical(fit$rankzz, as.integer(fx$rankzz))
  expect_equal(fit$condxx, fx$condxx, tolerance = stata_tol$coef)
  expect_equal(fit$condzz, fx$condzz, tolerance = stata_tol$coef)
  expect_equal(fit$ll, fx$ll, tolerance = stata_tol$coef)
  expect_equal(fit$diagnostics$ccev[1], fx$ccev1, tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$cdev[1], fx$cdev1, tolerance = stata_tol$stat)
})


# ==========================================================================
# OLS (no instruments)
# ==========================================================================
test_that("stored results match Stata: OLS", {
  data(card)
  fit <- ivreg2(lwage ~ exper + expersq + black + south + educ, data = card)
  fx <- read_sr_fixture("ols")

  expect_equal(fit$yy, fx$yy, tolerance = stata_tol$coef)
  expect_equal(fit$yyc, fx$yyc, tolerance = stata_tol$coef)
  expect_identical(fit$rank, as.integer(fx$rankxx))
  expect_null(fit$rankzz)
  expect_equal(fit$condxx, fx$condxx, tolerance = stata_tol$coef)
  expect_null(fit$condzz)
  expect_equal(fit$ll, fx$ll, tolerance = stata_tol$coef)
  # OLS: no ccev/cdev
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

test_that("glance condzz is NA for OLS", {
  data(card)
  fit <- ivreg2(lwage ~ educ + exper, data = card)
  gl <- glance(fit, diagnostics = TRUE)
  expect_true(is.na(gl$condzz))
  expect_true(is.na(gl$rankzz))
  expect_true(is.na(gl$ccev_min))
  expect_true(is.na(gl$cdev_min))
})
