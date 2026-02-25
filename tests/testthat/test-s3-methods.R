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


# ============================================================================
# terms
# ============================================================================

test_that("terms() returns regressors terms by default (OLS)", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  tt <- terms(fit)
  expect_s3_class(tt, "terms")
  expect_true(all(c("wt", "hp") %in% attr(tt, "term.labels")))
})

test_that("terms() returns all three components for IV model", {
  skip_if_not(exists("card"))
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  # Regressors
  tt_r <- terms(fit, component = "regressors")
  expect_s3_class(tt_r, "terms")
  expect_true("educ" %in% attr(tt_r, "term.labels"))

  # Instruments (excluded)
  tt_i <- terms(fit, component = "instruments")
  expect_s3_class(tt_i, "terms")
  expect_true("nearc4" %in% attr(tt_i, "term.labels"))

  # Full
  tt_f <- terms(fit, component = "full")
  expect_s3_class(tt_f, "terms")
})

test_that("terms(, component = 'instruments') returns NULL for OLS", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(terms(fit, component = "instruments"))
})


# ============================================================================
# model.matrix
# ============================================================================

test_that("model.matrix regressors dimensions match (x = TRUE)", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, x = TRUE)
  X <- model.matrix(fit, component = "regressors")
  expect_equal(nrow(X), nobs(fit))
  expect_equal(colnames(X), names(coef(fit)))
})

test_that("model.matrix regressors matches stored x$X", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, x = TRUE)
  X <- model.matrix(fit, component = "regressors")
  expect_identical(X, fit$x$X)
})

test_that("model.matrix instruments returns NULL for OLS", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, x = TRUE)
  expect_null(model.matrix(fit, component = "instruments"))
})

test_that("model.matrix projected equals X for OLS", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, x = TRUE)
  X <- model.matrix(fit, component = "regressors")
  Xhat <- model.matrix(fit, component = "projected")
  expect_equal(Xhat, X)
})

test_that("model.matrix IV dimensions (x = TRUE)", {
  skip_if_not(exists("card"))
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, x = TRUE)
  X <- model.matrix(fit, component = "regressors")
  Z <- model.matrix(fit, component = "instruments")
  Xhat <- model.matrix(fit, component = "projected")

  # X has K columns (regressors including intercept)
  expect_equal(nrow(X), nobs(fit))
  expect_equal(ncol(X), length(coef(fit)))

  # Z has L columns (exog + excluded IVs including intercept)
  expect_equal(nrow(Z), nobs(fit))
  expect_true(ncol(Z) >= ncol(X))  # overidentified or just-identified

  # Projected has same dims as X
  expect_equal(dim(Xhat), dim(X))
})

test_that("model.matrix projected differs from X for IV", {
  skip_if_not(exists("card"))
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, x = TRUE)
  X <- model.matrix(fit, component = "regressors")
  Xhat <- model.matrix(fit, component = "projected")
  # Projected endogenous regressors differ from raw
  expect_false(isTRUE(all.equal(X[, "educ"], Xhat[, "educ"])))
})

test_that("model.matrix reconstructs from model frame (x = FALSE)", {
  fit_x <- ivreg2(mpg ~ wt + hp, data = mtcars, x = TRUE, model = TRUE)
  fit_m <- ivreg2(mpg ~ wt + hp, data = mtcars, x = FALSE, model = TRUE)
  X_from_x <- model.matrix(fit_x, component = "regressors")
  X_from_m <- model.matrix(fit_m, component = "regressors")
  expect_equal(colnames(X_from_m), colnames(X_from_x))
  expect_equal(as.numeric(X_from_m), as.numeric(X_from_x))
})

test_that("model.matrix IV reconstructs from model frame (x = FALSE)", {
  skip_if_not(exists("card"))
  fit_x <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card, x = TRUE, model = TRUE)
  fit_m <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card, x = FALSE, model = TRUE)

  X_from_x <- model.matrix(fit_x, component = "regressors")
  X_from_m <- model.matrix(fit_m, component = "regressors")
  expect_equal(colnames(X_from_m), colnames(X_from_x))
  expect_equal(as.numeric(X_from_m), as.numeric(X_from_x))

  Z_from_x <- model.matrix(fit_x, component = "instruments")
  Z_from_m <- model.matrix(fit_m, component = "instruments")
  expect_equal(colnames(Z_from_m), colnames(Z_from_x))
  expect_equal(as.numeric(Z_from_m), as.numeric(Z_from_x))
})

test_that("model.matrix errors with x = FALSE, model = FALSE", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, x = FALSE, model = FALSE)
  expect_error(model.matrix(fit), "not enough information")
})

test_that("model.matrix reconstruction matches stored after collinearity", {
  # Create data with a collinear excluded instrument (z2 = z1)
  set.seed(42)
  n <- 200
  z1 <- rnorm(n)
  z2 <- z1  # perfectly collinear
  x <- z1 + rnorm(n)
  y <- x + rnorm(n)
  dat <- data.frame(y = y, x = x, z1 = z1, z2 = z2)

  # z2 should be dropped for collinearity
  fit_x <- suppressWarnings(suppressMessages(
    ivreg2(y ~ 1 | x | z1 + z2, data = dat, x = TRUE, model = TRUE)
  ))
  fit_m <- suppressWarnings(suppressMessages(
    ivreg2(y ~ 1 | x | z1 + z2, data = dat, x = FALSE, model = TRUE)
  ))

  # Verify z2 was dropped
  expect_true(length(fit_x$dropped_instruments) > 0L)

  # Reconstructed Z should match stored Z
  Z_from_x <- model.matrix(fit_x, component = "instruments")
  Z_from_m <- model.matrix(fit_m, component = "instruments")
  expect_equal(colnames(Z_from_m), colnames(Z_from_x))
  expect_equal(ncol(Z_from_m), ncol(Z_from_x))
  expect_equal(as.numeric(Z_from_m), as.numeric(Z_from_x))
})

test_that("model.matrix with factor variables", {
  dat <- mtcars
  dat$cyl_f <- factor(dat$cyl)
  fit <- ivreg2(mpg ~ wt + cyl_f, data = dat, x = TRUE)
  X <- model.matrix(fit, component = "regressors")
  expect_true("cyl_f6" %in% colnames(X))
  expect_true("cyl_f8" %in% colnames(X))
})

test_that("model.matrix weighted projected uses weights", {
  skip_if_not(exists("card"))
  fit_w <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card, weights = weight, x = TRUE)
  fit_u <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card, x = TRUE)
  Xhat_w <- model.matrix(fit_w, component = "projected")
  Xhat_u <- model.matrix(fit_u, component = "projected")
  # Weighted and unweighted projections should differ
  expect_false(isTRUE(all.equal(Xhat_w, Xhat_u)))
})


# ============================================================================
# update
# ============================================================================

test_that("update changes vcov type", {
  fit_iid <- ivreg2(mpg ~ wt + hp, data = mtcars)
  fit_hc1 <- update(fit_iid, vcov = "HC1")
  expect_equal(fit_hc1$vcov_type, "HC1")
  # Coefficients unchanged
  expect_equal(coef(fit_iid), coef(fit_hc1))
  # SEs differ
  expect_false(isTRUE(all.equal(vcov(fit_iid), vcov(fit_hc1))))
})

test_that("update changes formula (drop regressor)", {
  fit1 <- ivreg2(mpg ~ wt + hp, data = mtcars)
  fit2 <- update(fit1, . ~ . - hp)
  expect_equal(length(coef(fit2)), 2L)  # intercept + wt
  expect_false("hp" %in% names(coef(fit2)))
})

test_that("update evaluate = FALSE returns unevaluated call", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  cl <- update(fit, vcov = "HC0", evaluate = FALSE)
  expect_true(is.call(cl))
  expect_equal(cl$vcov, "HC0")
})

test_that("update changes data", {
  fit1 <- ivreg2(mpg ~ wt + hp, data = mtcars)
  mtcars2 <- mtcars[1:20, ]
  fit2 <- update(fit1, data = mtcars2)
  expect_equal(nobs(fit2), 20L)
})

test_that("update IV model changes vcov", {
  skip_if_not(exists("card"))
  fit_iid <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card)
  fit_hc1 <- update(fit_iid, vcov = "HC1")
  expect_equal(fit_hc1$vcov_type, "HC1")
  expect_equal(coef(fit_iid), coef(fit_hc1))
})
