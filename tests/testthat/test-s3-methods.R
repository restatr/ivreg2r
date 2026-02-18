# ============================================================================
# Tests: S3 methods for ivreg2 objects (Ticket G1)
# ============================================================================

# --- Helper: load Card data ---
card_path <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures", "card_data.csv"
)
if (file.exists(card_path)) {
  card <- read.csv(card_path)
}


# ============================================================================
# Simple extractors — OLS (vs lm on mtcars)
# ============================================================================

test_that("coef() returns named vector matching lm", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(coef(fit), coef(lm_fit), tolerance = .Machine$double.eps^0.5)
})

test_that("vcov() returns matrix matching lm", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(vcov(fit), vcov(lm_fit), tolerance = .Machine$double.eps^0.5)
})

test_that("residuals() returns vector matching lm", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(residuals(fit), residuals(lm_fit),
               tolerance = .Machine$double.eps^0.5)
})

test_that("fitted() returns vector matching lm", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fitted(fit), fitted(lm_fit),
               tolerance = .Machine$double.eps^0.5)
})

test_that("nobs() returns integer matching lm", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_identical(nobs(fit), nobs(lm_fit))
})

test_that("formula() returns the original formula", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  f <- formula(fit)
  expect_s3_class(f, "formula")
  # Should contain all three variables
  expect_true("mpg" %in% all.vars(f))
  expect_true("wt" %in% all.vars(f))
  expect_true("hp" %in% all.vars(f))
})

test_that("formula() forwards ... to formula.Formula (rhs subsetting)", {
  skip_if_not(file.exists(card_path), "Card data not available")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  # Full formula preserves multi-part structure
  f_full <- formula(fit)
  expect_true(grepl("\\|", deparse(f_full)))
  # rhs = 1 extracts only the exogenous part (no "|")
  f_rhs1 <- formula(fit, rhs = 1L)
  expect_false(grepl("\\|", deparse(f_rhs1)))
  expect_true(all(c("exper", "expersq", "black", "south") %in% all.vars(f_rhs1)))
})


# ============================================================================
# confint — OLS
# ============================================================================

test_that("confint with small=TRUE matches confint(lm())", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(confint(fit), confint(lm_fit),
               tolerance = .Machine$double.eps^0.5)
})

test_that("confint with small=FALSE uses qnorm (narrower than t)", {
  fit_small <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  fit_large <- ivreg2(mpg ~ wt + hp, data = mtcars, small = FALSE)
  ci_small <- confint(fit_small)
  ci_large <- confint(fit_large)
  # Normal CIs are narrower than t-CIs for same level
  width_small <- ci_small[, 2] - ci_small[, 1]
  width_large <- ci_large[, 2] - ci_large[, 1]
  expect_true(all(width_large < width_small))
})

test_that("confint small=FALSE matches manual qnorm computation", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = FALSE)
  cf <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  z <- qnorm(0.025)
  expected <- cbind(cf + z * se, cf - z * se)
  expect_equal(unname(confint(fit)), unname(expected),
               tolerance = .Machine$double.eps^0.5)
})

test_that("confint parm subsetting by name", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  ci <- confint(fit, parm = "wt")
  expect_equal(nrow(ci), 1L)
  expect_equal(rownames(ci), "wt")
})

test_that("confint parm subsetting by index", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  ci <- confint(fit, parm = 2)
  expect_equal(nrow(ci), 1L)
  expect_equal(rownames(ci), "wt")
})

test_that("confint level argument works", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  ci_90 <- confint(fit, level = 0.90)
  ci_99 <- confint(fit, level = 0.99)
  width_90 <- ci_90[, 2] - ci_90[, 1]
  width_99 <- ci_99[, 2] - ci_99[, 1]
  expect_true(all(width_99 > width_90))
})


# ============================================================================
# confint — IV
# ============================================================================

test_that("confint IV small=TRUE uses t distribution", {
  skip_if_not(exists("card"))
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, small = TRUE)
  cf <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  crit <- qt(0.025, df = fit$df.residual)
  expected <- cbind(cf + crit * se, cf - crit * se)
  expect_equal(unname(confint(fit)), unname(expected),
               tolerance = .Machine$double.eps^0.5)
})

test_that("confint IV small=FALSE uses qnorm", {
  skip_if_not(exists("card"))
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, small = FALSE)
  cf <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  crit <- qnorm(0.025)
  expected <- cbind(cf + crit * se, cf - crit * se)
  expect_equal(unname(confint(fit)), unname(expected),
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# predict
# ============================================================================

test_that("predict with no newdata returns fitted.values", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  expect_identical(predict(fit), fitted(fit))
})

test_that("predict OLS matches predict(lm())", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  nd <- data.frame(wt = c(2.5, 3.0, 3.5), hp = c(100, 150, 200))
  expect_equal(predict(fit, newdata = nd), predict(lm_fit, newdata = nd),
               tolerance = .Machine$double.eps^0.5)
})

test_that("predict with newdata matches manual X %*% coef", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  nd <- data.frame(wt = c(2.5, 3.0), hp = c(100, 150))
  X <- model.matrix(~ wt + hp, data = nd)
  expected <- drop(X %*% coef(fit))
  expect_equal(predict(fit, newdata = nd), expected,
               tolerance = .Machine$double.eps^0.5)
})

test_that("predict works with factor variables", {
  dat <- mtcars
  dat$cyl_f <- factor(dat$cyl)
  fit <- ivreg2(mpg ~ wt + cyl_f, data = dat, small = TRUE)
  nd <- data.frame(wt = 3.0, cyl_f = factor(c("4", "6", "8"),
                                              levels = c("4", "6", "8")))
  pred <- predict(fit, newdata = nd)
  expect_length(pred, 3L)
  expect_true(all(is.finite(pred)))
})


# ============================================================================
# summary / print.summary — OLS
# ============================================================================

test_that("summary returns object of class summary.ivreg2", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  s <- summary(fit)
  expect_s3_class(s, "summary.ivreg2")
})

test_that("summary coefficient table has 4 columns", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  s <- summary(fit)
  expect_equal(ncol(s$coef_table), 4L)
  expect_equal(nrow(s$coef_table), 3L)  # intercept, wt, hp
})

test_that("summary small=TRUE labels t value", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  s <- summary(fit)
  expect_true("t value" %in% colnames(s$coef_table))
  expect_true("Pr(>|t|)" %in% colnames(s$coef_table))
})

test_that("summary small=FALSE labels z value", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = FALSE)
  s <- summary(fit)
  expect_true("z value" %in% colnames(s$coef_table))
  expect_true("Pr(>|z|)" %in% colnames(s$coef_table))
})

test_that("print.summary OLS output contains expected strings", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  out <- capture.output(summary(fit))
  out_text <- paste(out, collapse = "\n")
  expect_match(out_text, "OLS Estimation")
  expect_match(out_text, "R-squared:")
  expect_match(out_text, "Root MSE:")
  # No diagnostics section for OLS
  expect_false(grepl("Weak identification", out_text))
  expect_false(grepl("Underidentification", out_text))
})


# ============================================================================
# summary / print.summary — IV
# ============================================================================

test_that("print.summary IV output contains 2SLS and diagnostics", {
  skip_if_not(exists("card"))
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, small = TRUE)
  out <- capture.output(summary(fit))
  out_text <- paste(out, collapse = "\n")
  expect_match(out_text, "2SLS Estimation")
  expect_match(out_text, "Weak identification")
  expect_match(out_text, "Underidentification")
  expect_match(out_text, "exactly identified")
  expect_match(out_text, "Instrumented:")
  expect_match(out_text, "Excluded instruments:")
})

test_that("print.summary IV robust shows KP stats", {
  skip_if_not(exists("card"))
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "HC1", small = TRUE)
  out <- capture.output(summary(fit))
  out_text <- paste(out, collapse = "\n")
  expect_match(out_text, "Kleibergen-Paap")
})

test_that("print.summary clustered shows cluster info", {
  skip_if_not(exists("card"))
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, clusters = ~smsa, small = TRUE)
  out <- capture.output(summary(fit))
  out_text <- paste(out, collapse = "\n")
  expect_match(out_text, "Cluster")
  expect_match(out_text, "smsa")
})

test_that("summary OLS p-values match lm summary", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  s <- summary(fit)
  lm_s <- summary(lm_fit)
  # Compare p-values
  expect_equal(s$coef_table[, 4], lm_s$coefficients[, 4],
               tolerance = .Machine$double.eps^0.5)
})

test_that("summary OLS t-statistics match lm summary", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  s <- summary(fit)
  lm_s <- summary(lm_fit)
  expect_equal(s$coef_table[, 3], lm_s$coefficients[, 3],
               tolerance = .Machine$double.eps^0.5)
})
