# ============================================================================
# Tests: Factor variable support (Ticket L1)
# ============================================================================
# Verify that factor variables (which expand to multiple dummy columns) work
# correctly in all formula positions: exogenous, endogenous, and excluded
# instruments. Uses synthetic data (not Stata fixtures) since R and Stata
# handle factor expansion with different naming conventions.


# --- Helper: generate synthetic data with factors ---
.make_factor_data <- function(n = 500, seed = 42) {
  set.seed(seed)
  region <- factor(sample(1:4, n, replace = TRUE))
  x1 <- rnorm(n)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  # Structural equation: y = 1 + x1 + endo + region dummies + error
  # First stage: endo = z1 + z2 + noise
  endo <- 0.5 * z1 + 0.3 * z2 + rnorm(n)
  y <- 1 + 0.5 * x1 + 0.8 * endo +
    0.3 * (region == "2") + 0.6 * (region == "3") +
    0.1 * (region == "4") + rnorm(n)
  data.frame(y = y, x1 = x1, endo = endo, z1 = z1, z2 = z2,
             region = region)
}


# ============================================================================
# Factor as exogenous regressor
# ============================================================================

test_that("factor exogenous regressor: diagnostics produce numeric results", {
  d <- .make_factor_data()

  fit <- ivreg2(y ~ x1 + region | endo | z1 + z2, data = d)

  # Coefficients should include region dummies
  cf <- coef(fit)
  expect_true("region2" %in% names(cf))
  expect_true("region3" %in% names(cf))
  expect_true("region4" %in% names(cf))
  expect_true(all(is.finite(cf)))

  # Diagnostics should all be numeric (no NA from failed match())
  s <- summary(fit)
  expect_true(is.finite(s$diagnostics$weak_id$stat))
  expect_true(is.finite(s$diagnostics$underid$stat))
  expect_true(is.finite(s$diagnostics$anderson_rubin$f_stat))
  expect_true(is.finite(s$diagnostics$endogeneity$stat))
})

test_that("factor exogenous regressor: predict(newdata) works", {
  d <- .make_factor_data()
  fit <- ivreg2(y ~ x1 + region | endo | z1 + z2, data = d)
  preds <- predict(fit, newdata = d[1:5, ])
  expect_length(preds, 5L)
  expect_true(all(is.finite(preds)))
  # Should match fitted values for the same rows
  expect_equal(unname(preds), unname(fitted(fit)[1:5]), tolerance = 1e-10)
})


# ============================================================================
# Factor as endogenous regressor
# ============================================================================

test_that("factor endogenous regressor: diagnostics produce numeric results", {
  d <- .make_factor_data()
  # Treat region as endogenous (needs enough instruments)
  # region expands to 3 dummies → need >= 3 excluded instruments
  d$z3 <- rnorm(nrow(d))
  d$z4 <- rnorm(nrow(d))

  fit <- ivreg2(y ~ x1 | region | z1 + z2 + z3 + z4, data = d)

  # $endogenous should contain the term label, not dummy names

  expect_equal(fit$endogenous, "region")

  # Coefficients should include region dummies
  cf <- coef(fit)
  expect_true(all(c("region2", "region3", "region4") %in% names(cf)))

  # Diagnostics should be numeric
  s <- summary(fit)
  expect_true(is.finite(s$diagnostics$weak_id$stat))
  expect_true(is.finite(s$diagnostics$underid$stat))
  expect_true(is.finite(s$diagnostics$anderson_rubin$f_stat))
  expect_true(is.finite(s$diagnostics$endogeneity$stat))

  # Summary footer should not list factor dummies as included instruments
  footer <- capture.output(print(s))
  incl_line <- grep("Included instruments:", footer, value = TRUE)
  expect_false(grepl("region2", incl_line))
  expect_false(grepl("region3", incl_line))
  expect_false(grepl("region4", incl_line))

  # First-stage: one entry per dummy column
  expect_equal(length(fit$first_stage), 3L)
  expect_true(all(c("region2", "region3", "region4") %in%
                     names(fit$first_stage)))
  for (nm in names(fit$first_stage)) {
    expect_true(is.finite(fit$first_stage[[nm]]$f_stat))
  }
})

test_that("factor endogenous: predict(newdata) works", {
  d <- .make_factor_data()
  d$z3 <- rnorm(nrow(d))
  d$z4 <- rnorm(nrow(d))
  fit <- ivreg2(y ~ x1 | region | z1 + z2 + z3 + z4, data = d)
  preds <- predict(fit, newdata = d[1:5, ])
  expect_length(preds, 5L)
  expect_true(all(is.finite(preds)))
})


# ============================================================================
# Factor as excluded instrument
# ============================================================================

test_that("factor excluded instrument: identification tests work", {
  d <- .make_factor_data()
  # region as excluded instrument (3 dummies → overidentified)
  fit <- ivreg2(y ~ x1 | endo | region, data = d)

  # Diagnostics should be numeric
  s <- summary(fit)
  expect_true(is.finite(s$diagnostics$weak_id$stat))
  expect_true(is.finite(s$diagnostics$underid$stat))
  expect_true(is.finite(s$diagnostics$anderson_rubin$f_stat))
  expect_true(is.finite(s$diagnostics$endogeneity$stat))

  # Overidentified (3 excluded instruments, 1 endogenous) → overid test should work
  expect_true(is.finite(s$diagnostics$overid$stat))
  expect_equal(s$diagnostics$overid$df, 2L)

  # $instruments should contain the term label
  expect_equal(fit$instruments, "region")
})

test_that("factor excluded instrument: predict(newdata) works", {
  d <- .make_factor_data()
  fit <- ivreg2(y ~ x1 | endo | region, data = d)
  preds <- predict(fit, newdata = d[1:5, ])
  expect_length(preds, 5L)
  expect_true(all(is.finite(preds)))
})


# ============================================================================
# Collinear factor dummies
# ============================================================================

test_that("collinear factor dummies are handled correctly", {
  d <- .make_factor_data()
  # Create a variable perfectly collinear with one of the region dummies
  d$region2_dup <- as.numeric(d$region == "2")

  # This should drop the collinear column and still work
  expect_warning(
    fit <- ivreg2(y ~ x1 + region2_dup + region | endo | z1 + z2, data = d),
    "collinear"
  )
  expect_true(all(is.finite(coef(fit))))

  # predict(newdata) should work even with dropped columns
  preds <- predict(fit, newdata = d[1:5, ])
  expect_length(preds, 5L)
  expect_true(all(is.finite(preds)))
  expect_equal(unname(preds), unname(fitted(fit)[1:5]), tolerance = 1e-10)
})


# ============================================================================
# Non-default contrasts on factor endogenous
# ============================================================================

test_that("non-default contrasts on factor endo: predict uses stored contrasts", {
  d <- .make_factor_data()
  d$z3 <- rnorm(nrow(d))
  d$z4 <- rnorm(nrow(d))

  # Use sum contrasts instead of default treatment contrasts
  d$region_sum <- d$region
  contrasts(d$region_sum) <- contr.sum(4)

  fit <- ivreg2(y ~ x1 | region_sum | z1 + z2 + z3 + z4, data = d)

  # Contrasts should be stored
  expect_false(is.null(fit$contrasts))
  expect_true("region_sum" %in% names(fit$contrasts))

  # predict(newdata) should use the stored contrasts
  preds <- predict(fit, newdata = d[1:5, ])
  expect_length(preds, 5L)
  expect_true(all(is.finite(preds)))
  expect_equal(unname(preds), unname(fitted(fit)[1:5]), tolerance = 1e-10)
})


# ============================================================================
# endog/orthog with factor terms
# ============================================================================

test_that("endog with factor term expands to all dummies", {
  d <- .make_factor_data()
  # Two endogenous: endo (numeric) + region (factor, 3 dummies)
  # Need >= 4 excluded instruments
  d$z3 <- rnorm(nrow(d))
  d$z4 <- rnorm(nrow(d))
  d$z5 <- rnorm(nrow(d))

  fit <- ivreg2(y ~ x1 | endo + region | z1 + z2 + z3 + z4 + z5,
                data = d, endog = "region")

  # Endogeneity test should have df = 3 (3 dummies for region)
  expect_equal(fit$diagnostics$endogeneity$df, 3L)
  expect_true(is.finite(fit$diagnostics$endogeneity$stat))
})

test_that("endog with numeric term still works with factors present", {
  d <- .make_factor_data()
  d$z3 <- rnorm(nrow(d))
  d$z4 <- rnorm(nrow(d))
  d$z5 <- rnorm(nrow(d))

  fit <- ivreg2(y ~ x1 | endo + region | z1 + z2 + z3 + z4 + z5,
                data = d, endog = "endo")

  # Endogeneity test should have df = 1 (just endo)
  expect_equal(fit$diagnostics$endogeneity$df, 1L)
  expect_true(is.finite(fit$diagnostics$endogeneity$stat))
})


# ============================================================================
# LIML with factor variables
# ============================================================================

test_that("LIML with factor endogenous works", {
  d <- .make_factor_data()
  d$z3 <- rnorm(nrow(d))
  d$z4 <- rnorm(nrow(d))

  fit <- ivreg2(y ~ x1 | region | z1 + z2 + z3 + z4, data = d,
                method = "liml")

  expect_true(all(is.finite(coef(fit))))
  expect_true(is.finite(fit$lambda))
  expect_true(fit$lambda >= 1)
})

test_that("LIML with factor excluded instrument works", {
  d <- .make_factor_data()

  fit <- ivreg2(y ~ x1 | endo | region, data = d, method = "liml")

  expect_true(all(is.finite(coef(fit))))
  expect_true(is.finite(fit$lambda))
})


# ============================================================================
# Robust VCV with factor variables
# ============================================================================

test_that("HC1 with factor variables works", {
  d <- .make_factor_data()

  fit <- ivreg2(y ~ x1 + region | endo | z1 + z2, data = d,
                vcov = "HC1", small = TRUE)

  expect_true(all(is.finite(coef(fit))))
  s <- summary(fit)
  expect_true(is.finite(s$diagnostics$weak_id$stat))
  expect_true(is.finite(s$diagnostics$weak_id_robust$stat))
})


# ============================================================================
# Reduced form with factor variables
# ============================================================================

test_that("reduced form with factor endogenous works", {
  d <- .make_factor_data()
  d$z3 <- rnorm(nrow(d))
  d$z4 <- rnorm(nrow(d))

  fit <- ivreg2(y ~ x1 | region | z1 + z2 + z3 + z4, data = d,
                reduced_form = "system")

  expect_false(is.null(fit$reduced_form))
  expect_equal(fit$reduced_form$mode, "system")
})
