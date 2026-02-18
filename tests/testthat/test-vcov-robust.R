# ============================================================================
# Tests: HC0/HC1 robust sandwich VCV (Ticket C1)
# ============================================================================
#
# Fixture naming convention: Stata fixtures named "hc1" were generated with
# ivreg2's `, robust` option (no finite-sample correction = HC0 in sandwich
# terminology). Fixtures named "hc1_small" were generated with `, robust small`
# (with N/(N-K) correction = HC1 in sandwich terminology).

# --- Helper: load Card data and fixtures ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)
card_path <- file.path(fixture_dir, "card_data.csv")

if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

# Helper: read VCV fixture and convert Stata names to R names
read_vcov_fixture <- function(path) {
  fixture <- read.csv(path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)
  vcov_cols <- grep("^vcov_", names(fixture), value = TRUE)
  V <- as.matrix(fixture[, vcov_cols])
  rownames(V) <- r_names
  col_stata <- sub("^vcov_", "", vcov_cols)
  colnames(V) <- ifelse(col_stata == "_cons", "(Intercept)", col_stata)
  V
}

# Helper: compare VCV matrices element-wise
expect_vcov_equal <- function(V_r, V_stata, tol = stata_tol$vcov) {
  shared <- intersect(rownames(V_r), rownames(V_stata))
  for (rn in shared) {
    for (cn in shared) {
      expect_equal(
        V_r[rn, cn], V_stata[rn, cn],
        tolerance = tol,
        info = paste("VCV mismatch:", rn, cn)
      )
    }
  }
}

# ============================================================================
# card_just_id: HC0 matches Stata `, robust` (no correction)
# ============================================================================

test_that("2SLS HC0 VCV matches Stata card_just_id robust fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  vcov_path <- file.path(fixture_dir, "card_just_id_vcov_hc1.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "HC0")
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

# ============================================================================
# card_just_id: HC1 matches Stata `, robust small` (with N/(N-K) correction)
# ============================================================================

test_that("2SLS HC1 VCV matches Stata card_just_id robust small fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  vcov_path <- file.path(fixture_dir, "card_just_id_vcov_hc1_small.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "HC1")
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

# ============================================================================
# card_overid: HC0 matches Stata `, robust`
# ============================================================================

test_that("2SLS HC0 VCV matches Stata card_overid robust fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  vcov_path <- file.path(fixture_dir, "card_overid_vcov_hc1.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "HC0")
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

# ============================================================================
# card_overid: HC1 matches Stata `, robust small`
# ============================================================================

test_that("2SLS HC1 VCV matches Stata card_overid robust small fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  vcov_path <- file.path(fixture_dir, "card_overid_vcov_hc1_small.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "HC1")
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

# ============================================================================
# HC1 = HC0 * N/(N-K) (the correction is the only difference)
# ============================================================================

test_that("HC1 VCV equals HC0 VCV scaled by N/(N-K)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_hc0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, vcov = "HC0")
  fit_hc1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, vcov = "HC1")
  N <- fit_hc0$nobs
  K <- length(coef(fit_hc0))
  expect_equal(fit_hc1$vcov, fit_hc0$vcov * (N / (N - K)),
               tolerance = .Machine$double.eps^0.5)
})

# ============================================================================
# small flag does not affect robust VCV (only sigma, test distributions)
# ============================================================================

test_that("small flag does not change HC0 VCV", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_a <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                   data = card, vcov = "HC0", small = FALSE)
  fit_b <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                   data = card, vcov = "HC0", small = TRUE)
  expect_equal(fit_a$vcov, fit_b$vcov, tolerance = .Machine$double.eps^0.5)
})

test_that("small flag does not change HC1 VCV", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_a <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                   data = card, vcov = "HC1", small = FALSE)
  fit_b <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                   data = card, vcov = "HC1", small = TRUE)
  expect_equal(fit_a$vcov, fit_b$vcov, tolerance = .Machine$double.eps^0.5)
})

# ============================================================================
# Coefficients unchanged by VCV type
# ============================================================================

test_that("coefficients are identical for iid, HC0, HC1", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_iid <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, vcov = "iid")
  fit_hc0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, vcov = "HC0")
  fit_hc1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, vcov = "HC1")
  expect_equal(coef(fit_hc0), coef(fit_iid),
               tolerance = .Machine$double.eps^0.5)
  expect_equal(coef(fit_hc1), coef(fit_iid),
               tolerance = .Machine$double.eps^0.5)
})

# ============================================================================
# OLS: cross-validate against sandwich::vcovHC
# ============================================================================

test_that("OLS HC1 matches sandwich::vcovHC(type='HC1')", {
  skip_if_not_installed("sandwich")

  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC1")
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  V_sandwich <- sandwich::vcovHC(lm_fit, type = "HC1")

  expect_equal(fit$vcov, V_sandwich, tolerance = 1e-10)
})

test_that("OLS HC0 matches sandwich::vcovHC(type='HC0')", {
  skip_if_not_installed("sandwich")

  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC0")
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  V_sandwich <- sandwich::vcovHC(lm_fit, type = "HC0")

  expect_equal(fit$vcov, V_sandwich, tolerance = 1e-10)
})

# ============================================================================
# vcov_type field is stored correctly
# ============================================================================

test_that("vcov_type reports HC0", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC0")
  expect_identical(fit$vcov_type, "HC0")
})

test_that("vcov_type reports HC1", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC1")
  expect_identical(fit$vcov_type, "HC1")
})

# ============================================================================
# VCV matrix is symmetric
# ============================================================================

test_that("HC1 VCV is symmetric", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "HC1")
  expect_equal(fit$vcov, t(fit$vcov), tolerance = .Machine$double.eps^0.5)
})
