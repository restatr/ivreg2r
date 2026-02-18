# ============================================================================
# Tests: ivreg2() OLS estimation (Ticket A2)
# ============================================================================

# --- Helper: load Card data ---
card_path <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures", "card_data.csv"
)
if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

# ============================================================================
# mtcars: basic OLS with intercept
# ============================================================================

test_that("coefficients match lm()", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(coef(fit), coef(lm_fit), tolerance = .Machine$double.eps^0.5)
})

test_that("residuals match lm()", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$residuals, residuals(lm_fit), tolerance = .Machine$double.eps^0.5)
})

test_that("fitted values match lm()", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$fitted.values, fitted(lm_fit), tolerance = .Machine$double.eps^0.5)
})

test_that("R-squared matches lm()", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$r.squared, summary(lm_fit)$r.squared,
               tolerance = .Machine$double.eps^0.5)
})

test_that("adjusted R-squared matches lm()", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$adj.r.squared, summary(lm_fit)$adj.r.squared,
               tolerance = .Machine$double.eps^0.5)
})

test_that("vcov with small=TRUE matches vcov(lm())", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$vcov, vcov(lm_fit), tolerance = .Machine$double.eps^0.5)
})

test_that("vcov with small=FALSE equals vcov(lm) * (N-K)/N", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = FALSE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  N <- nobs(lm_fit)
  K <- length(coef(lm_fit))
  expected <- vcov(lm_fit) * (N - K) / N
  expect_equal(fit$vcov, expected, tolerance = .Machine$double.eps^0.5)
})

test_that("sigma with small=TRUE matches summary(lm)$sigma", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$sigma, summary(lm_fit)$sigma,
               tolerance = .Machine$double.eps^0.5)
})

test_that("sigma with small=FALSE equals sqrt(RSS/N)", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = FALSE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expected <- sqrt(sum(residuals(lm_fit)^2) / nobs(lm_fit))
  expect_equal(fit$sigma, expected, tolerance = .Machine$double.eps^0.5)
})

test_that("RSS matches lm()", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$rss, sum(residuals(lm_fit)^2),
               tolerance = .Machine$double.eps^0.5)
})

test_that("df.residual matches lm()", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_identical(fit$df.residual, lm_fit$df.residual)
})

test_that("rank matches lm()", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_identical(fit$rank, lm_fit$rank)
})

test_that("nobs is correct", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_identical(fit$nobs, nobs(lm_fit))
})

test_that("class is ivreg2", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_s3_class(fit, "ivreg2")
})

test_that("diagnostics and first_stage are NULL for OLS", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$diagnostics)
  expect_null(fit$first_stage)
})

test_that("endogenous and instruments are empty for OLS", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_identical(fit$endogenous, character(0))
  expect_identical(fit$instruments, character(0))
})

# ============================================================================
# No-intercept model
# ============================================================================

test_that("no-intercept R-squared uses uncentered TSS", {
  fit <- ivreg2(mpg ~ 0 + wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ 0 + wt + hp, data = mtcars)
  # lm() uses uncentered TSS when there's no intercept
  expected_r2 <- summary(lm_fit)$r.squared
  expect_equal(fit$r.squared, expected_r2, tolerance = .Machine$double.eps^0.5)
})

test_that("no-intercept adjusted R-squared formula: 1 - (1 - r2) * N/(N-K)", {
  fit <- ivreg2(mpg ~ 0 + wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ 0 + wt + hp, data = mtcars)
  expect_equal(fit$adj.r.squared, summary(lm_fit)$adj.r.squared,
               tolerance = .Machine$double.eps^0.5)
})

test_that("no-intercept coefficients match lm()", {
  fit <- ivreg2(mpg ~ 0 + wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ 0 + wt + hp, data = mtcars)
  expect_equal(coef(fit), coef(lm_fit), tolerance = .Machine$double.eps^0.5)
})

# ============================================================================
# Card data OLS
# ============================================================================

test_that("Card data OLS coefficients match lm()", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fit <- ivreg2(lwage ~ educ + exper + expersq + black + south, data = card)
  lm_fit <- lm(lwage ~ educ + exper + expersq + black + south, data = card)
  expect_equal(coef(fit), coef(lm_fit), tolerance = .Machine$double.eps^0.5)
})

test_that("Card data OLS vcov (small=TRUE) matches lm()", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fit <- ivreg2(lwage ~ educ + exper + expersq + black + south,
                data = card, small = TRUE)
  lm_fit <- lm(lwage ~ educ + exper + expersq + black + south, data = card)
  expect_equal(fit$vcov, vcov(lm_fit), tolerance = .Machine$double.eps^0.5)
})

# ============================================================================
# NA handling
# ============================================================================

test_that("NAs are dropped and results match lm()", {
  d <- mtcars
  d$mpg[c(1, 5, 10)] <- NA
  d$wt[15] <- NA
  fit <- ivreg2(mpg ~ wt + hp, data = d)
  lm_fit <- lm(mpg ~ wt + hp, data = d)
  expect_equal(coef(fit), coef(lm_fit), tolerance = .Machine$double.eps^0.5)
  expect_identical(fit$nobs, nobs(lm_fit))
  expect_equal(fit$residuals, residuals(lm_fit), tolerance = .Machine$double.eps^0.5)
})

# ============================================================================
# Subset argument
# ============================================================================

test_that("subset argument works correctly", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, subset = cyl == 6)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars, subset = cyl == 6)
  expect_equal(coef(fit), coef(lm_fit), tolerance = .Machine$double.eps^0.5)
  expect_identical(fit$nobs, nobs(lm_fit))
  expect_equal(fit$residuals, residuals(lm_fit), tolerance = .Machine$double.eps^0.5)
})

# ============================================================================
# model/x/y flags
# ============================================================================

test_that("model=TRUE stores model frame", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, model = TRUE)
  expect_false(is.null(fit$model))
  expect_s3_class(fit$model, "data.frame")
})

test_that("model=FALSE omits model frame", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, model = FALSE)
  expect_null(fit$model)
})

test_that("x=TRUE stores model matrices", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, x = TRUE)
  expect_false(is.null(fit$x))
  expect_true(is.matrix(fit$x$X))
  expect_null(fit$x$Z)  # OLS has no Z
})

test_that("x=FALSE omits model matrices", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, x = FALSE)
  expect_null(fit$x)
})

test_that("y=TRUE stores response vector", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, y = TRUE)
  expect_false(is.null(fit$y))
  expect_equal(length(fit$y), nrow(mtcars))
})

test_that("y=FALSE omits response vector", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, y = FALSE)
  expect_null(fit$y)
})

# ============================================================================
# Call stored
# ============================================================================

test_that("call is stored correctly", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_true(is.call(fit$call))
  expect_equal(fit$call[[1L]], quote(ivreg2))
})

# ============================================================================
# IV formula runs (no longer guarded)
# ============================================================================

test_that("IV formula runs without error", {
  set.seed(42)
  d <- data.frame(y = rnorm(20), x = rnorm(20), e = rnorm(20), z = rnorm(20))
  fit <- ivreg2(y ~ x | e | z, data = d)
  expect_s3_class(fit, "ivreg2")
  expect_length(fit$endogenous, 1L)
})

# ============================================================================
# Print method
# ============================================================================

test_that("print method shows OLS Estimation", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_output(print(fit), "OLS Estimation")
})

# ============================================================================
# Argument validation
# ============================================================================

test_that("invalid vcov value errors", {
  expect_error(ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC2"),
               "not yet implemented")
})

test_that("non-formula clusters argument errors", {
  expect_error(ivreg2(mpg ~ wt + hp, data = mtcars, clusters = "cyl"),
               "one-sided formula")
})

test_that("weights argument accepts valid weights", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp)
  expect_s3_class(fit, "ivreg2")
})

test_that("invalid small argument errors", {
  expect_error(ivreg2(mpg ~ wt + hp, data = mtcars, small = "yes"),
               "must be TRUE or FALSE")
})

# ============================================================================
# vcov_type and small are stored
# ============================================================================

test_that("vcov_type and small are stored", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  expect_identical(fit$vcov_type, "iid")
  expect_true(fit$small)
})
