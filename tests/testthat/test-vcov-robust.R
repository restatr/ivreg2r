# ============================================================================
# Tests: Robust sandwich VCV (Ticket C1)
#
# M-10 re-base (2026-07-06): the card_just_id/card_overid robust fixture cells
# are re-based onto the mroz justid D5a cell (mroz_justid_formula) and the
# hf H41 overid cell (help.txt:1325); the retired *_small fixture cells have
# no direct replacement fixture — their heir is the "robust+small VCV equals
# robust VCV scaled by N/(N-K)" identity test below, extended to cover both
# the justid and overid models.
#
# Fixture naming convention: Stata fixtures named "hc1" were generated with
# ivreg2's `, robust` option (no finite-sample correction). Fixtures named
# "hc1_small" were generated with `, robust small` (with N/(N-K) correction).

data(mroz, package = "ivreg2r")

mroz_justid_formula <- lwage ~ exper + expersq | educ | age
mroz_overid_formula <- lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6

# ============================================================================
# mroz_justid: robust matches Stata `, robust` (no correction)
# ============================================================================

test_that("2SLS robust VCV matches Stata mroz_justid robust fixture", {
  fit <- ivreg2(mroz_justid_formula, data = mroz, vcov = "robust")
  expect_coef_fixture(fit, "mroz_justid_coef_robust.csv")
  expect_vcov_fixture(fit, "mroz_justid_vcov_robust.csv")
})

# ============================================================================
# hf overid (H41): robust matches Stata `, robust`
# ============================================================================

test_that("2SLS robust VCV matches Stata hf_mroz H41 robust fixture", {
  fit <- ivreg2(mroz_overid_formula, data = mroz, vcov = "robust")
  expect_vcov_fixture(fit, "hf_mroz_vcov_H41.csv")
})

# ============================================================================
# small=TRUE VCV = small=FALSE VCV * N/(N-K) (the correction is the difference)
#
# This identity is the fixture-free retirement of the card_just_id/
# card_overid *_small VCV fixture cells (2026-07-06 byte-diff: robust VCVs
# under `small` scale by exactly N/(N-K)); it is checked on both the justid
# and overid models so both retired cells have an identity heir.
# ============================================================================

test_that("robust+small VCV equals robust VCV scaled by N/(N-K): justid", {
  fit_base <- ivreg2(mroz_justid_formula, data = mroz, vcov = "robust")
  fit_small <- ivreg2(mroz_justid_formula, data = mroz, vcov = "robust",
                      small = TRUE)
  N <- fit_base$nobs
  K <- length(coef(fit_base))
  expect_equal(fit_small$vcov, fit_base$vcov * (N / (N - K)),
               tolerance = .Machine$double.eps^0.5)
})

test_that("robust+small VCV equals robust VCV scaled by N/(N-K): overid", {
  fit_base <- ivreg2(mroz_overid_formula, data = mroz, vcov = "robust")
  fit_small <- ivreg2(mroz_overid_formula, data = mroz, vcov = "robust",
                      small = TRUE)
  N <- fit_base$nobs
  K <- length(coef(fit_base))
  expect_equal(fit_small$vcov, fit_base$vcov * (N / (N - K)),
               tolerance = .Machine$double.eps^0.5)
})

# ============================================================================
# Coefficients unchanged by VCV type and small flag
# ============================================================================

test_that("coefficients are identical for iid, robust, robust+small", {
  fit_iid <- ivreg2(mroz_justid_formula, data = mroz, vcov = "iid")
  fit_robust <- ivreg2(mroz_justid_formula, data = mroz, vcov = "robust")
  fit_robust_small <- ivreg2(mroz_justid_formula, data = mroz,
                             vcov = "robust", small = TRUE)
  expect_equal(coef(fit_robust), coef(fit_iid),
               tolerance = .Machine$double.eps^0.5)
  expect_equal(coef(fit_robust_small), coef(fit_iid),
               tolerance = .Machine$double.eps^0.5)
})

# ============================================================================
# OLS: cross-validate against sandwich::vcovHC
# ============================================================================

test_that("OLS robust+small matches sandwich::vcovHC(type='HC1')", {
  skip_if_not_installed("sandwich")

  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "robust", small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  V_sandwich <- sandwich::vcovHC(lm_fit, type = "HC1")

  expect_equal(fit$vcov, V_sandwich, tolerance = 1e-10)
})

test_that("OLS robust matches sandwich::vcovHC(type='HC0')", {
  skip_if_not_installed("sandwich")

  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "robust")
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  V_sandwich <- sandwich::vcovHC(lm_fit, type = "HC0")

  expect_equal(fit$vcov, V_sandwich, tolerance = 1e-10)
})

# ============================================================================
# vcov_type field is stored correctly
# ============================================================================

test_that("vcov_type reports robust", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "robust")
  expect_identical(fit$vcov_type, "robust")
})

test_that("vcov_type reports robust with small", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "robust", small = TRUE)
  expect_identical(fit$vcov_type, "robust")
})

# ============================================================================
# HC0/HC1 are rejected with informative error
# ============================================================================

test_that("vcov = 'HC0' produces informative error", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC0"),
    'vcov = "HC0" is no longer supported.*Use vcov = "robust"'
  )
})

test_that("vcov = 'HC1' produces informative error", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC1"),
    'vcov = "HC1" is no longer supported.*Use vcov = "robust"'
  )
})

# ============================================================================
# VCV matrix is symmetric
# ============================================================================

test_that("robust VCV is symmetric", {
  fit <- ivreg2(mroz_justid_formula, data = mroz, vcov = "robust")
  expect_equal(fit$vcov, t(fit$vcov), tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# sim_no_constant: robust VCV (noconstant model)
# ============================================================================

sim_noconst_path <- fixture_path("sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

test_that("2SLS robust VCV matches Stata sim_no_constant robust fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  vcov_path <- fixture_path("sim_no_constant_vcov_hc1.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, vcov = "robust")
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

test_that("2SLS robust+small VCV matches Stata sim_no_constant robust small fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  vcov_path <- fixture_path("sim_no_constant_vcov_hc1_small.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst,
                vcov = "robust", small = TRUE)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})
