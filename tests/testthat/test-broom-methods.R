# ============================================================================
# Tests: broom methods — tidy(), glance(), augment() for ivreg2 objects (G2)
# ============================================================================

# --- Helper: load Card data ---
card_path <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures", "card_data.csv"
)
if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)


# ============================================================================
# tidy() — OLS
# ============================================================================

test_that("tidy OLS returns tibble with correct columns (conf.int = TRUE)", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  td <- tidy(fit)
  expect_s3_class(td, "tbl_df")
  expect_named(td, c("term", "estimate", "std.error", "statistic",
                      "p.value", "conf.low", "conf.high"))
  expect_equal(nrow(td), 3L)
})

test_that("tidy OLS conf.int = FALSE omits CI columns", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  td <- tidy(fit, conf.int = FALSE)
  expect_named(td, c("term", "estimate", "std.error", "statistic", "p.value"))
})

test_that("tidy OLS row count = number of coefficients", {
  fit <- ivreg2(mpg ~ wt + hp + disp, data = mtcars, small = TRUE)
  td <- tidy(fit)
  expect_equal(nrow(td), length(coef(fit)))
})

test_that("tidy OLS values match lm (small = TRUE, mtcars)", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  lm_s <- summary(lm_fit)$coefficients
  td <- tidy(fit, conf.int = FALSE)
  expect_equal(td$estimate, unname(lm_s[, "Estimate"]),
               tolerance = .Machine$double.eps^0.5)
  expect_equal(td$std.error, unname(lm_s[, "Std. Error"]),
               tolerance = .Machine$double.eps^0.5)
  expect_equal(td$statistic, unname(lm_s[, "t value"]),
               tolerance = .Machine$double.eps^0.5)
  expect_equal(td$p.value, unname(lm_s[, "Pr(>|t|)"]),
               tolerance = .Machine$double.eps^0.5)
})

test_that("tidy conf.level = 0.99 produces wider CIs than 0.95", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  td_95 <- tidy(fit, conf.level = 0.95)
  td_99 <- tidy(fit, conf.level = 0.99)
  width_95 <- td_95$conf.high - td_95$conf.low
  width_99 <- td_99$conf.high - td_99$conf.low
  expect_true(all(width_99 > width_95))
})


# ============================================================================
# tidy() — IV
# ============================================================================

test_that("tidy IV estimates/SEs match Stata fixtures (card_just_id iid)", {
  skip_if_not(file.exists(card_path), "Card data not available")
  coef_path <- file.path(fixture_dir, "card_just_id_coef_iid.csv")
  skip_if_not(file.exists(coef_path), "Stata fixture not available")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  stata <- read.csv(coef_path, stringsAsFactors = FALSE)
  td <- tidy(fit, conf.int = FALSE)

  # Match by term name — Stata uses _cons for intercept
  stata$term[stata$term == "_cons"] <- "(Intercept)"
  stata <- stata[match(td$term, stata$term), ]

  expect_equal(td$estimate, stata$estimate, tolerance = stata_tol$coef)
  expect_equal(td$std.error, stata$std_error, tolerance = stata_tol$se)
})


# ============================================================================
# glance() — OLS
# ============================================================================

test_that("glance OLS returns 1-row tibble with 17 columns", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  gl <- glance(fit)
  expect_s3_class(gl, "tbl_df")
  expect_equal(nrow(gl), 1L)
  expect_equal(ncol(gl), 30L)
})

test_that("glance OLS has correct column names", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  gl <- glance(fit)
  expected_names <- c("r.squared", "adj.r.squared", "sigma", "statistic",
                      "p.value", "df", "df.residual", "nobs", "vcov_type",
                      "method", "lambda", "kclass_value", "fuller_parameter",
                      "n_clusters1", "n_clusters2",
                      "weak_id_stat", "weak_id_robust_stat",
                      "underid_stat", "underid_p",
                      "overid_stat", "overid_p",
                      "endogeneity_stat", "endogeneity_p",
                      "stock_wright_stat", "stock_wright_p",
                      "stock_wright_df",
                      "orthog_stat", "orthog_p",
                      "rf_f_stat", "rf_f_p")
  expect_named(gl, expected_names)
})

test_that("glance OLS column types are correct", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  gl <- glance(fit)
  expect_type(gl$r.squared, "double")
  expect_type(gl$adj.r.squared, "double")
  expect_type(gl$sigma, "double")
  expect_type(gl$nobs, "integer")
  expect_type(gl$df.residual, "integer")
  expect_type(gl$vcov_type, "character")
})

test_that("glance OLS: all diagnostic columns are NA", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  gl <- glance(fit)
  diag_cols <- c("weak_id_stat", "weak_id_robust_stat",
                 "underid_stat", "underid_p",
                 "overid_stat", "overid_p",
                 "endogeneity_stat", "endogeneity_p",
                 "orthog_stat", "orthog_p")
  for (col in diag_cols) {
    expect_true(is.na(gl[[col]]), info = paste(col, "should be NA for OLS"))
  }
})

test_that("glance OLS r.squared, sigma, nobs match lm", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  lm_s <- summary(lm_fit)
  gl <- glance(fit)
  expect_equal(gl$r.squared, lm_s$r.squared,
               tolerance = .Machine$double.eps^0.5)
  expect_equal(gl$adj.r.squared, lm_s$adj.r.squared,
               tolerance = .Machine$double.eps^0.5)
  expect_equal(gl$sigma, lm_s$sigma,
               tolerance = .Machine$double.eps^0.5)
  expect_equal(gl$nobs, nobs(lm_fit))
})


# ============================================================================
# glance() — IV
# ============================================================================

test_that("glance IV just-identified iid: weak_id_stat present, overid NA", {
  skip_if_not(file.exists(card_path), "Card data not available")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  gl <- glance(fit)
  expect_false(is.na(gl$weak_id_stat))
  # Just-identified: overid should be NA
  expect_true(is.na(gl$overid_stat))
  expect_true(is.na(gl$overid_p))
  # IID: robust stat should be NA
  expect_true(is.na(gl$weak_id_robust_stat))
})

test_that("glance IV just-identified iid: weak_id_stat matches Stata CD F", {
  skip_if_not(file.exists(card_path), "Card data not available")
  diag_path <- file.path(fixture_dir, "card_just_id_diagnostics_iid.csv")
  skip_if_not(file.exists(diag_path), "Stata fixture not available")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  stata <- read.csv(diag_path, stringsAsFactors = FALSE)
  gl <- glance(fit)
  expect_equal(gl$weak_id_stat, stata$cdf, tolerance = stata_tol$stat)
})

test_that("glance IV overidentified HC1: overid and robust stats present", {
  skip_if_not(file.exists(card_path), "Card data not available")
  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
    data = card, vcov = "HC1"
  )
  gl <- glance(fit)
  expect_false(is.na(gl$overid_stat))
  expect_false(is.na(gl$overid_p))
  expect_false(is.na(gl$weak_id_robust_stat))
})


# ============================================================================
# augment() — OLS
# ============================================================================

test_that("augment OLS returns tibble with .fitted and .resid", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  aug <- augment(fit)
  expect_s3_class(aug, "tbl_df")
  expect_true(".fitted" %in% names(aug))
  expect_true(".resid" %in% names(aug))
})

test_that("augment OLS nrow = nobs when data = NULL", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  aug <- augment(fit)
  expect_equal(nrow(aug), nobs(fit))
})

test_that("augment .fitted matches fitted(), .resid matches residuals()", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  aug <- augment(fit)
  expect_equal(aug$.fitted, unname(fitted(fit)))
  expect_equal(aug$.resid, unname(residuals(fit)))
})

test_that("augment with data argument works", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  aug <- augment(fit, data = mtcars)
  expect_s3_class(aug, "tbl_df")
  expect_equal(nrow(aug), nrow(mtcars))
  expect_true(".fitted" %in% names(aug))
  expect_true(".resid" %in% names(aug))
})

test_that("augment errors when model = FALSE and no data supplied", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, model = FALSE)
  expect_error(augment(fit), "model = FALSE")
})

test_that("augment model = FALSE works with data argument", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, model = FALSE)
  aug <- augment(fit, data = mtcars)
  expect_s3_class(aug, "tbl_df")
  expect_equal(nrow(aug), nrow(mtcars))
})


# ============================================================================
# augment() — NA alignment with na.exclude
# ============================================================================

test_that("augment handles na.exclude alignment", {
  dat <- mtcars
  dat$mpg[c(3, 7)] <- NA
  fit <- ivreg2(mpg ~ wt + hp, data = dat, na.action = na.exclude,
                small = TRUE)
  # With data argument: should have NAs at rows 3 and 7
  aug <- augment(fit, data = dat)
  expect_equal(nrow(aug), nrow(dat))
  expect_true(is.na(aug$.fitted[3]))
  expect_true(is.na(aug$.fitted[7]))
  expect_true(is.na(aug$.resid[3]))
  expect_true(is.na(aug$.resid[7]))
  # Non-missing rows should be finite
  expect_true(all(is.finite(aug$.fitted[-c(3, 7)])))
})

test_that("augment with no data and na.exclude: nrow = nobs (no NAs)", {
  dat <- mtcars
  dat$mpg[c(3, 7)] <- NA
  fit <- ivreg2(mpg ~ wt + hp, data = dat, na.action = na.exclude,
                small = TRUE)
  aug <- augment(fit)
  # Model frame only has complete cases
  expect_equal(nrow(aug), nobs(fit))
  expect_equal(nrow(aug), nrow(dat) - 2L)
  # No NAs in model-frame augmentation
  expect_false(anyNA(aug$.fitted))
  expect_false(anyNA(aug$.resid))
})


# ============================================================================
# augment() — IV
# ============================================================================

test_that("augment IV returns expected columns", {
  skip_if_not(file.exists(card_path), "Card data not available")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  aug <- augment(fit)
  expect_s3_class(aug, "tbl_df")
  expect_true(".fitted" %in% names(aug))
  expect_true(".resid" %in% names(aug))
  expect_equal(nrow(aug), nobs(fit))
})
