# Tests for predict(..., se.fit = TRUE) on the canonical mroz H31 base.
# Prediction standard errors: SE = sqrt(diag(X %*% V %*% t(X)))

# === Data: mroz estimation sample (Stata's e(sample) = lwage non-missing) ===
data(mroz, package = "ivreg2r")
mroz_lw <- mroz[!is.na(mroz$lwage), ]

# === Helper ===
read_stdp_fixture <- function(prefix, suffix) {
  path <- fixture_path(paste0(prefix, "_stdp_", suffix, ".csv"))
  if (!file.exists(path)) return(NULL)
  read.csv(path)
}

compare_stdp <- function(fit, fixture, tol = stata_tol$se, label = "") {
  stata <- read_stdp_fixture(fixture$prefix, fixture$suffix)
  skip_if(is.null(stata), "stdp fixture not found")
  pred <- predict(fit, se.fit = TRUE)
  expect_equal(pred$se.fit, stata$se_fit, tolerance = tol,
               info = paste(label, "se.fit mismatch"))
  expect_equal(unname(pred$fit), stata$fit, tolerance = stata_tol$coef,
               info = paste(label, "fit mismatch"))
}

compare_stdp_sub <- function(fit, fixture, newdata, tol = stata_tol$se,
                              label = "") {
  stata <- read_stdp_fixture(fixture$prefix, paste0(fixture$suffix, "_sub"))
  skip_if(is.null(stata), "stdp subset fixture not found")
  pred <- predict(fit, newdata = newdata, se.fit = TRUE)
  expect_equal(pred$se.fit, stata$se_fit, tolerance = tol,
               info = paste(label, "se.fit newdata mismatch"))
  expect_equal(unname(pred$fit), stata$fit, tolerance = stata_tol$coef,
               info = paste(label, "fit newdata mismatch"))
}


# ==========================================================================
# Stata parity: mroz — IID
# Canonical base H31 (help.txt:1274); `predict, stdp` per D5a (no worked
# example in the help file).
# ==========================================================================
test_that("stdp matches Stata: mroz IID", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw)
  compare_stdp(fit, list(prefix = "mroz", suffix = "iid"))
})

test_that("stdp newdata matches Stata: mroz IID", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw)
  compare_stdp_sub(fit, list(prefix = "mroz", suffix = "iid"),
                   newdata = mroz_lw[1:10, ])
})


# ==========================================================================
# Stata parity: mroz — robust
# H41/H46 robust spec (help.txt:1325/1350) on the H31 base.
# ==========================================================================
test_that("stdp matches Stata: mroz robust", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust")
  compare_stdp(fit, list(prefix = "mroz", suffix = "robust"))
})

test_that("stdp newdata matches Stata: mroz robust", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust")
  compare_stdp_sub(fit, list(prefix = "mroz", suffix = "robust"),
                   newdata = mroz_lw[1:10, ])
})


# ==========================================================================
# Stata parity: mroz — robust small
# D5a option-variation (small) on the H41 robust spec; verifies the
# small-sample e(V) convention propagates into stdp.
# ==========================================================================
test_that("stdp matches Stata: mroz robust small", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust", small = TRUE)
  compare_stdp(fit, list(prefix = "mroz", suffix = "robust_small"))
})

test_that("stdp newdata matches Stata: mroz robust small", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust", small = TRUE)
  compare_stdp_sub(fit, list(prefix = "mroz", suffix = "robust_small"),
                   newdata = mroz_lw[1:10, ])
})


# ==========================================================================
# Partial blocked
# ==========================================================================
test_that("se.fit errors after partialling", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, partial = "exper")
  expect_error(predict(fit, se.fit = TRUE), "partial")
})


# ==========================================================================
# se.fit = FALSE returns bare vector (default behavior preserved)
# ==========================================================================
test_that("se.fit = FALSE returns numeric vector, not list", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                data = card)
  pred <- predict(fit)
  expect_true(is.numeric(pred))
  expect_false(is.list(pred))
  expect_equal(pred, fitted(fit))
})


# ==========================================================================
# se.fit = TRUE returns list with correct structure
# ==========================================================================
test_that("se.fit = TRUE returns list with fit and se.fit components", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                data = card)
  pred <- predict(fit, se.fit = TRUE)
  expect_true(is.list(pred))
  expect_named(pred, c("fit", "se.fit"))
  expect_true(is.numeric(pred$fit))
  expect_true(is.numeric(pred$se.fit))
  expect_equal(length(pred$fit), nobs(fit))
  expect_equal(length(pred$se.fit), nobs(fit))
})


# ==========================================================================
# Consistency: se.fit equals manual sqrt(diag(X V X'))
# (deliberately the naive n x n formula, independent of .compute_se_fit's
# rowSums algebra; run on mroz, n = 428, so the n x n matrix stays small)
# ==========================================================================
test_that("se.fit matches manual computation", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust", model = TRUE)
  pred <- predict(fit, se.fit = TRUE)

  X <- model.matrix(fit)
  V <- vcov(fit)
  se_manual <- sqrt(diag(X %*% V %*% t(X)))
  expect_equal(pred$se.fit, unname(se_manual), tolerance = 1e-12)
})


# ==========================================================================
# OLS model also works with se.fit
# ==========================================================================
test_that("se.fit works for OLS models", {
  fit <- ivreg2(lwage ~ exper + expersq + educ, data = mroz_lw, model = TRUE)
  pred <- predict(fit, se.fit = TRUE)
  expect_true(is.list(pred))
  expect_equal(length(pred$se.fit), nobs(fit))
  expect_true(all(pred$se.fit >= 0))

  X <- model.matrix(fit)
  V <- vcov(fit)
  se_manual <- sqrt(diag(X %*% V %*% t(X)))
  expect_equal(pred$se.fit, unname(se_manual), tolerance = 1e-12)
})


# ==========================================================================
# newdata se.fit
# ==========================================================================
test_that("se.fit with newdata works", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                data = card)
  sub <- card[1:10, ]
  pred <- predict(fit, newdata = sub, se.fit = TRUE)
  expect_true(is.list(pred))
  expect_equal(length(pred$se.fit), 10L)
  expect_true(all(pred$se.fit >= 0))
})


# ==========================================================================
# All se.fit values are non-negative
# ==========================================================================
test_that("se.fit values are non-negative", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust", small = TRUE)
  pred <- predict(fit, se.fit = TRUE)
  expect_true(all(pred$se.fit >= 0))
})


# ==========================================================================
# Model without stored model/x errors informatively
# ==========================================================================
test_that("se.fit errors when model and x not stored", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                data = card, model = FALSE, x = FALSE)
  expect_error(predict(fit, se.fit = TRUE), "Refit")
})


# ==========================================================================
# na.exclude: NAs reinserted at omitted positions, matching predict.lm
# ==========================================================================
test_that("na.exclude reinserts NAs in predict, fitted, residuals", {
  dat <- mtcars
  dat$mpg[c(3, 7)] <- NA
  fit <- ivreg2(mpg ~ wt + hp, data = dat, na.action = na.exclude,
                small = TRUE, model = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = dat, na.action = na.exclude)

  # All observation-level extractors return length = nrow(original data)
  expect_equal(length(fitted(fit)), nrow(dat))
  expect_equal(length(residuals(fit)), nrow(dat))
  expect_equal(length(predict(fit)), nrow(dat))

  # NAs at the omitted positions
  expect_true(is.na(fitted(fit)[3]))
  expect_true(is.na(fitted(fit)[7]))
  expect_true(is.na(residuals(fit)[3]))
  expect_true(is.na(residuals(fit)[7]))
  expect_true(is.na(predict(fit)[3]))
  expect_true(is.na(predict(fit)[7]))

  # Non-NA values match lm
  expect_equal(fitted(fit), fitted(lm_fit), tolerance = 1e-12)
  expect_equal(residuals(fit), residuals(lm_fit), tolerance = 1e-12)
  expect_equal(predict(fit), predict(lm_fit), tolerance = 1e-12)
})

test_that("se.fit=TRUE reinserts NAs under na.exclude", {
  dat <- mtcars
  dat$mpg[c(3, 7)] <- NA
  fit <- ivreg2(mpg ~ wt + hp, data = dat, na.action = na.exclude,
                small = TRUE, model = TRUE)

  pred_plain <- predict(fit)
  pred_se <- predict(fit, se.fit = TRUE)

  # fit and se.fit both have NAs at positions 3 and 7
  expect_equal(length(pred_se$fit), nrow(dat))
  expect_equal(length(pred_se$se.fit), nrow(dat))
  expect_equal(pred_se$fit, pred_plain)
  expect_true(is.na(pred_se$se.fit[3]))
  expect_true(is.na(pred_se$se.fit[7]))
  expect_true(all(pred_se$se.fit[-c(3, 7)] >= 0))
})

test_that("na.omit (default) returns complete-case vectors", {
  dat <- mtcars
  dat$mpg[c(3, 7)] <- NA
  fit <- ivreg2(mpg ~ wt + hp, data = dat, small = TRUE)

  # Default na.omit: no reinsertion, length = nobs
  expect_equal(length(fitted(fit)), nobs(fit))
  expect_equal(length(residuals(fit)), nobs(fit))
  expect_equal(length(predict(fit)), nobs(fit))
  expect_false(anyNA(fitted(fit)))
})
