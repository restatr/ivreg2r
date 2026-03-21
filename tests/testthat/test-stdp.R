# Tests for predict(..., se.fit = TRUE) — Ticket Q1
# Prediction standard errors: SE = sqrt(diag(X %*% V %*% t(X)))

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
# Stata parity: Card overid — IID
# ==========================================================================
test_that("stdp matches Stata: Card overid IID", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card)
  compare_stdp(fit, list(prefix = "card_overid", suffix = "iid"))
})

test_that("stdp newdata matches Stata: Card overid IID", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card)
  compare_stdp_sub(fit, list(prefix = "card_overid", suffix = "iid"),
                   newdata = card[1:10, ])
})


# ==========================================================================
# Stata parity: Card overid — robust (HC0)
# ==========================================================================
test_that("stdp matches Stata: Card overid HC0", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, vcov = "robust")
  compare_stdp(fit, list(prefix = "card_overid", suffix = "hc0"))
})


# ==========================================================================
# Stata parity: Card overid — robust small (HC1)
# ==========================================================================
test_that("stdp matches Stata: Card overid HC1 small", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, vcov = "robust", small = TRUE)
  compare_stdp(fit, list(prefix = "card_overid", suffix = "hc1_small"))
})


# ==========================================================================
# Stata parity: Card overid — cluster
# ==========================================================================
test_that("stdp matches Stata: Card overid cluster", {
  data(card, package = "ivreg2r")
  fit <- muffle_rank_warnings(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, clusters = ~smsa)
  )
  compare_stdp(fit, list(prefix = "card_overid", suffix = "cluster"))
})


# ==========================================================================
# Stata parity: Card overid — cluster small
# ==========================================================================
test_that("stdp matches Stata: Card overid cluster small", {
  data(card, package = "ivreg2r")
  fit <- muffle_rank_warnings(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, clusters = ~smsa, small = TRUE)
  )
  compare_stdp(fit, list(prefix = "card_overid", suffix = "cluster_small"))
})

test_that("stdp newdata matches Stata: Card overid cluster small", {
  data(card, package = "ivreg2r")
  fit <- muffle_rank_warnings(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, clusters = ~smsa, small = TRUE)
  )
  compare_stdp_sub(fit, list(prefix = "card_overid", suffix = "cluster_small"),
                   newdata = card[1:10, ])
})


# ==========================================================================
# Stata parity: Card overid — aweight
# ==========================================================================
test_that("stdp matches Stata: Card overid aweight", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, weights = weight)
  compare_stdp(fit, list(prefix = "card_overid", suffix = "aw"))
})


# ==========================================================================
# Stata parity: Card overid — LIML
# ==========================================================================
test_that("stdp matches Stata: Card overid LIML", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, method = "liml")
  compare_stdp(fit, list(prefix = "card_overid", suffix = "liml"))
})


# ==========================================================================
# Stata parity: Card overid — Fuller(1)
# ==========================================================================
test_that("stdp matches Stata: Card overid Fuller(1)", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, method = "liml", fuller = 1)
  compare_stdp(fit, list(prefix = "card_overid", suffix = "fuller1"))
})


# ==========================================================================
# Stata parity: Card just-identified — IID
# ==========================================================================
test_that("stdp matches Stata: Card justid IID", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                data = card)
  compare_stdp(fit, list(prefix = "card_justid", suffix = "iid"))
})

test_that("stdp newdata matches Stata: Card justid IID", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                data = card)
  compare_stdp_sub(fit, list(prefix = "card_justid", suffix = "iid"),
                   newdata = card[1:10, ])
})


# ==========================================================================
# Stata parity: Card just-identified — robust (HC0)
# ==========================================================================
test_that("stdp matches Stata: Card justid HC0", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                data = card, vcov = "robust")
  compare_stdp(fit, list(prefix = "card_justid", suffix = "hc0"))
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
# ==========================================================================
test_that("se.fit matches manual computation", {
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, vcov = "robust", model = TRUE)
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
  data(card, package = "ivreg2r")
  fit <- ivreg2(lwage ~ exper + expersq + educ, data = card, model = TRUE)
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
  data(card, package = "ivreg2r")
  fit <- muffle_rank_warnings(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, clusters = ~smsa, small = TRUE)
  )
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
  expect_error(predict(fit, se.fit = TRUE), "refit")
})


# ==========================================================================
# na.exclude consistency: se.fit=TRUE fit matches se.fit=FALSE
# ==========================================================================
test_that("se.fit=TRUE fit matches predict(fit) under na.exclude", {
  dat <- mtcars
  dat$mpg[c(3, 7)] <- NA
  fit <- ivreg2(mpg ~ wt + hp, data = dat, na.action = na.exclude,
                small = TRUE, model = TRUE)
  pred_plain <- predict(fit)
  pred_se <- predict(fit, se.fit = TRUE)
  expect_equal(length(pred_se$fit), length(pred_plain))
  expect_equal(pred_se$fit, pred_plain)
  expect_equal(length(pred_se$se.fit), nobs(fit))
  expect_true(all(pred_se$se.fit >= 0))
})
