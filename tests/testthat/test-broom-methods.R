# ============================================================================
# Tests: broom methods — tidy(), glance(), augment() for ivreg2 objects (G2)
#
# M-10 re-base (2026-07-06): card fits are re-based onto the bundled mroz
# justid D5a cell (mroz_justid_formula) and the hf H31 overid cell
# (help.txt:1274), retiring the card_just_id fixtures.
# ============================================================================

data(mroz, package = "ivreg2r")

# mroz_justid_formula and mroz_overid_formula come from helper-fixtures.R.


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

test_that("tidy IV estimates/SEs match Stata fixtures (mroz_justid iid)", {
  coef_path <- fixture_path("mroz_justid_coef_iid.csv")

  fit <- ivreg2(mroz_justid_formula, data = mroz)
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

test_that("glance OLS returns 1-row tibble with expected columns", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  gl <- glance(fit)
  expect_s3_class(gl, "tbl_df")
  expect_equal(nrow(gl), 1L)
  expect_equal(ncol(gl), 15L)
})

test_that("glance OLS has correct column names", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  gl <- glance(fit)
  expected_names <- c("r.squared", "adj.r.squared", "sigma", "statistic",
                      "p.value", "df", "df.residual", "nobs", "vcov_type",
                      "weak_id_stat", "weak_id_robust_stat",
                      "underid_stat", "underid_p",
                      "overid_stat", "overid_p")
  expect_named(gl, expected_names)
})

test_that("glance drops config flags and stored results (reachable on object)", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  gl <- glance(fit)
  dropped <- c("small", "weight_type", "method", "lambda", "kclass_value",
               "coviv", "center", "psd", "kernel", "bw", "kiefer", "sw",
               "cue_convergence", "partial_ct", "yy", "yyc", "condxx", "ll",
               "endogeneity_stat", "orthog_stat", "redundancy_stat",
               "stock_wright_stat", "ar_overid_lr_stat", "rf_f_stat")
  for (col in dropped) {
    expect_false(col %in% names(gl),
                 info = paste(col, "should not be a glance column"))
  }
  # The dropped quantities remain reachable on the fitted object
  expect_type(fit$small, "logical")
  expect_identical(fit$method, "ols")
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
                 "overid_stat", "overid_p")
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
  fit <- ivreg2(mroz_justid_formula, data = mroz)
  gl <- glance(fit)
  expect_false(is.na(gl$weak_id_stat))
  # Just-identified: overid should be NA
  expect_true(is.na(gl$overid_stat))
  expect_true(is.na(gl$overid_p))
  # IID: robust stat should be NA
  expect_true(is.na(gl$weak_id_robust_stat))
})

test_that("glance IV just-identified iid: weak_id_stat matches Stata CD F", {
  diag_path <- fixture_path("mroz_justid_diagnostics_iid.csv")

  fit <- ivreg2(mroz_justid_formula, data = mroz)
  stata <- read.csv(diag_path, stringsAsFactors = FALSE)
  gl <- glance(fit)
  expect_equal(gl$weak_id_stat, stata$cdf, tolerance = stata_tol$stat)
})

test_that("glance IV overidentified HC1: overid and robust stats present", {
  # Overid (H31/H41 arc): sargan/df > 0 here, unlike the justid model above,
  # so overid_stat is meaningfully non-NA.
  fit <- ivreg2(mroz_overid_formula, data = mroz, vcov = "robust")
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
  fit <- ivreg2(mroz_justid_formula, data = mroz)
  aug <- augment(fit)
  expect_s3_class(aug, "tbl_df")
  expect_true(".fitted" %in% names(aug))
  expect_true(".resid" %in% names(aug))
  expect_equal(nrow(aug), nobs(fit))
})


# ============================================================================
# glance() — diagnostics argument
# ============================================================================

test_that("glance diagnostics = FALSE returns compact tibble", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  gl <- glance(fit, diagnostics = FALSE)
  expect_s3_class(gl, "tbl_df")
  expect_equal(nrow(gl), 1L)
  expect_equal(ncol(gl), 9L)
  # Diagnostic test columns should be absent
  diag_cols <- c("weak_id_stat", "weak_id_robust_stat",
                 "underid_stat", "underid_p",
                 "overid_stat", "overid_p")
  for (col in diag_cols) {
    expect_false(col %in% names(gl),
                 info = paste(col, "should not be present with diagnostics = FALSE"))
  }
})

test_that("glance diagnostics = TRUE matches default", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  gl_default <- glance(fit)
  gl_true <- glance(fit, diagnostics = TRUE)
  expect_identical(gl_default, gl_true)
})

test_that("glance diagnostics = FALSE retains model-fit columns", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  gl <- glance(fit, diagnostics = FALSE)
  retained <- c("r.squared", "adj.r.squared", "sigma", "statistic",
                "p.value", "df", "df.residual", "nobs", "vcov_type")
  for (col in retained) {
    expect_true(col %in% names(gl),
                info = paste(col, "should be present with diagnostics = FALSE"))
  }
})


# ============================================================================
# glance() — small column
# ============================================================================

# ============================================================================
# tidy() — exponentiate
# ============================================================================

test_that("tidy exponentiate = TRUE exponentiates estimates", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  td_raw <- tidy(fit, conf.int = FALSE)
  td_exp <- tidy(fit, conf.int = FALSE, exponentiate = TRUE)
  expect_equal(td_exp$estimate, exp(td_raw$estimate))
})

test_that("tidy exponentiate = TRUE exponentiates conf.low and conf.high", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  td_raw <- tidy(fit)
  td_exp <- tidy(fit, exponentiate = TRUE)
  expect_equal(td_exp$conf.low, exp(td_raw$conf.low))
  expect_equal(td_exp$conf.high, exp(td_raw$conf.high))
})

test_that("tidy exponentiate = TRUE does NOT transform std.error", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  td_raw <- tidy(fit, conf.int = FALSE)
  td_exp <- tidy(fit, conf.int = FALSE, exponentiate = TRUE)
  expect_equal(td_exp$std.error, td_raw$std.error)
})

test_that("tidy exponentiate = FALSE (default) is unchanged", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  td_default <- tidy(fit)
  td_false <- tidy(fit, exponentiate = FALSE)
  expect_identical(td_default, td_false)
})

test_that("tidy exponentiate = TRUE works with conf.int = FALSE", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  td_exp <- tidy(fit, conf.int = FALSE, exponentiate = TRUE)
  expect_false("conf.low" %in% names(td_exp))
  expect_false("conf.high" %in% names(td_exp))
  td_raw <- tidy(fit, conf.int = FALSE)
  expect_equal(td_exp$estimate, exp(td_raw$estimate))
})


# ============================================================================
# modelsummary integration
# ============================================================================

test_that("modelsummary works with a single ivreg2 object", {
  skip_if_not_installed("modelsummary")
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  out <- modelsummary::modelsummary(fit, output = "data.frame")
  expect_s3_class(out, "data.frame")
  expect_true(nrow(out) > 0L)
})

test_that("modelsummary works with a named list of models", {
  skip_if_not_installed("modelsummary")
  fit1 <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  fit2 <- ivreg2(mpg ~ wt + hp + disp, data = mtcars, small = TRUE)
  out <- modelsummary::modelsummary(
    list("Model 1" = fit1, "Model 2" = fit2), output = "data.frame"
  )
  expect_s3_class(out, "data.frame")
  expect_true(nrow(out) > 0L)
})

test_that("modelsummary works with IV model", {
  skip_if_not_installed("modelsummary")
  fit <- ivreg2(mroz_justid_formula, data = mroz, vcov = "robust")
  out <- modelsummary::modelsummary(fit, output = "data.frame")
  expect_s3_class(out, "data.frame")
  expect_true(nrow(out) > 0L)
})


# ============================================================================
# small flag — now off glance(), carried on the fitted object
# ============================================================================

test_that("small flag is stored on the fitted object", {
  fit_small <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  fit_large <- ivreg2(mpg ~ wt + hp, data = mtcars, small = FALSE)
  expect_true(fit_small$small)
  expect_false(fit_large$small)
  expect_type(fit_small$small, "logical")
})
