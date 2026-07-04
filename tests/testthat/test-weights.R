# ============================================================================
# Tests: Analytic Weights (Ticket C3)
# ============================================================================

# card_wt (card data + deterministic region/fwt columns) comes from
# helper-fixtures.R; card-based tests use it directly (bundled data).


# ============================================================================
# Weighted OLS vs lm(): exact match (iid VCV)
# ============================================================================

test_that("weighted OLS coefficients match lm() with weights", {
  fit_iv2 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp)
  fit_lm <- lm(mpg ~ wt + hp, data = mtcars, weights = disp)
  expect_equal(coef(fit_iv2), coef(fit_lm),
               tolerance = .Machine$double.eps^0.5)
})

test_that("weighted OLS SEs match lm() with weights (iid, small=TRUE)", {
  fit_iv2 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
                    small = TRUE)
  fit_lm <- lm(mpg ~ wt + hp, data = mtcars, weights = disp)
  se_iv2 <- sqrt(diag(fit_iv2$vcov))
  se_lm <- sqrt(diag(vcov(fit_lm)))
  expect_equal(se_iv2, se_lm, tolerance = .Machine$double.eps^0.5)
})

test_that("weighted OLS sigma is scale-invariant in weights", {
  # Stata aweight convention: normalizing weights to sum to N makes sigma
  # independent of weight scale. This differs from lm() where sigma scales
  # with sqrt(c) when weights are multiplied by c.
  fit1 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp, small = TRUE)
  d2 <- mtcars
  d2$disp10 <- mtcars$disp * 10
  fit2 <- ivreg2(mpg ~ wt + hp, data = d2, weights = disp10, small = TRUE)
  expect_equal(fit1$sigma, fit2$sigma, tolerance = .Machine$double.eps^0.5)
})

test_that("weighted OLS sigma matches Stata convention (differs from lm by sqrt(N/sum(w)))", {
  fit_iv2 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
                    small = TRUE)
  fit_lm <- lm(mpg ~ wt + hp, data = mtcars, weights = disp)
  N <- nrow(mtcars)
  correction <- sqrt(N / sum(mtcars$disp))
  expect_equal(fit_iv2$sigma, summary(fit_lm)$sigma * correction,
               tolerance = .Machine$double.eps^0.5)
})

test_that("weighted OLS R-squared matches lm() summary", {
  fit_iv2 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
                    small = TRUE)
  fit_lm <- lm(mpg ~ wt + hp, data = mtcars, weights = disp)
  expect_equal(fit_iv2$r.squared, summary(fit_lm)$r.squared,
               tolerance = .Machine$double.eps^0.5)
  expect_equal(fit_iv2$adj.r.squared, summary(fit_lm)$adj.r.squared,
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Weighted OLS vs lm(): HC1 VCV
# ============================================================================

test_that("weighted OLS robust+small VCV matches sandwich::vcovHC(type='HC1')", {
  skip_if_not_installed("sandwich")
  fit_iv2 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
                    vcov = "robust", small = TRUE)
  fit_lm <- lm(mpg ~ wt + hp, data = mtcars, weights = disp)
  V_sand <- sandwich::vcovHC(fit_lm, type = "HC1")
  expect_equal(fit_iv2$vcov, V_sand, tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Weighted 2SLS vs Stata: aweight cell-means parity (card_aw_cells_*)
# ============================================================================
# Collapsed data: collapse (mean) lwage (count) n=lwage, by(educ exper
# expersq black south nearc4 nearc2) on the micro card data (1030 cells;
# family M-11, re-based). Fit [aw=n] on the collapsed data — coefficients
# must recover the unweighted micro-data fit (see the identity test below);
# V differs by construction and is only compared against the Stata fixture.

aw_cells_by <- lwage ~ educ + exper + expersq + black + south + nearc4 + nearc2
card_cells <- aggregate(aw_cells_by, data = card_wt, FUN = mean)
card_cells$n <- aggregate(aw_cells_by, data = card_wt, FUN = length)$lwage

test_aw_cells_config <- function(suffix, vcov_arg, small_arg) {
  coef_path <- fixture_path(paste0("card_aw_cells_coef_", suffix, ".csv"))
  vcov_path <- fixture_path(paste0("card_aw_cells_vcov_", suffix, ".csv"))
  diag_path <- fixture_path(paste0("card_aw_cells_diagnostics_", suffix, ".csv"))

  skip_if(!file.exists(coef_path))

  stata_coef <- read_coef_fixture(coef_path)
  stata_vcov <- read_vcov_fixture(vcov_path)
  stata_diag <- read_diagnostics(diag_path)

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card_cells, weights = n,
                vcov = vcov_arg, small = small_arg)

  expect_stata_parity_core(fit, stata_coef, stata_vcov, stata_diag)

  # R-squared
  expect_equal(fit$r.squared, stata_diag$r2, tolerance = stata_tol$coef)

  # Overidentification: Sargan under iid, Hansen J under robust
  if (!is.na(stata_diag$sargan)) {
    expect_equal(fit$diagnostics$overid$stat, stata_diag$sargan,
                 tolerance = stata_tol$stat)
  }
  if (!is.na(stata_diag$j)) {
    expect_equal(fit$diagnostics$overid$stat, stata_diag$j,
                 tolerance = stata_tol$stat)
  }
}

for (cell in list(
  list(suffix = "iid",       vcov_arg = "iid",    small_arg = FALSE),
  list(suffix = "iid_small", vcov_arg = "iid",    small_arg = TRUE),
  list(suffix = "hc1",       vcov_arg = "robust", small_arg = FALSE),
  list(suffix = "hc1_small", vcov_arg = "robust", small_arg = TRUE)
)) {
  test_that(paste("aweight cell-means overid matches Stata:", cell$suffix), {
    test_aw_cells_config(cell$suffix, cell$vcov_arg, cell$small_arg)
  })
}

test_that("aweight cell-means recover the micro-data fit", {
  # Mirrors the .do generator's self-check: coefficients from the collapsed
  # [aw=n] iid fit equal the unweighted micro-data fit (V intentionally
  # differs by construction, so only coefficients are compared here).
  fit_cells <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                      data = card_cells, weights = n)
  fit_micro <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                      data = card_wt)
  expect_equal(coef(fit_cells), coef(fit_micro))
})


# ============================================================================
# Weight validation
# ============================================================================

test_that("negative weights give error", {
  d <- mtcars
  d$w <- c(-1, rep(1, nrow(d) - 1))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w),
    "strictly positive"
  )
})

test_that("zero weights give error", {
  d <- mtcars
  d$w <- c(0, rep(1, nrow(d) - 1))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w),
    "strictly positive"
  )
})

test_that("NA weights give error", {
  d <- mtcars
  d$w <- c(NA, rep(1, nrow(d) - 1))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w, na.action = na.pass),
    "finite and non-missing"
  )
})

test_that("Inf weights give error", {
  d <- mtcars
  d$w <- c(Inf, rep(1, nrow(d) - 1))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w),
    "finite and non-missing"
  )
})

test_that("NULL weights produce unweighted results", {
  fit_null <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit_null$weights)
  expect_equal(coef(fit_null), coef(lm(mpg ~ wt + hp, data = mtcars)),
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Edge cases
# ============================================================================

test_that("equal weights produce same results as unweighted", {
  d <- mtcars
  d$w <- rep(1, nrow(d))
  fit_w <- ivreg2(mpg ~ wt + hp, data = d, weights = w, small = TRUE)
  fit_u <- ivreg2(mpg ~ wt + hp, data = d, small = TRUE)
  expect_equal(coef(fit_w), coef(fit_u),
               tolerance = .Machine$double.eps^0.5)
  expect_equal(fit_w$vcov, fit_u$vcov,
               tolerance = .Machine$double.eps^0.5)
  expect_equal(fit_w$sigma, fit_u$sigma,
               tolerance = .Machine$double.eps^0.5)
  expect_equal(fit_w$r.squared, fit_u$r.squared,
               tolerance = .Machine$double.eps^0.5)
})

test_that("weighted OLS 1-part formula works", {
  fit <- ivreg2(mpg ~ wt + hp + cyl, data = mtcars, weights = disp)
  expect_s3_class(fit, "ivreg2")
  expect_length(coef(fit), 4L)
  expect_equal(fit$weights, mtcars$disp)
})

test_that("weighted IV with clustering works", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card_wt, weights = weight, clusters = ~ region)
  expect_s3_class(fit, "ivreg2")
  expect_identical(fit$vcov_type, "CL")
  expect_false(is.null(fit$weights))
})

test_that("weights are stored in the returned object", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp)
  expect_equal(fit$weights, mtcars$disp)
})
