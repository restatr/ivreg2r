# ==========================================================================
# Tests for IV collinearity metadata (Ticket A3)
#
# Verifies that .parse_formula() and ivreg2() correctly track metadata when
# collinear endogenous regressors or instruments are dropped/reclassified.
# ==========================================================================

make_collin_data <- function(n = 200) {
  set.seed(99)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  x1 <- rnorm(n)
  endo1 <- 0.5 * z1 + rnorm(n)
  endo2 <- 0.3 * z2 + rnorm(n)
  y <- 1 + x1 + endo1 + 0.5 * endo2 + rnorm(n)
  data.frame(y = y, x1 = x1, endo1 = endo1, endo2 = endo2,
             z1 = z1, z2 = z2)
}

# ==========================================================================
# 1. Collinear endogenous regressor dropped
# ==========================================================================
test_that("collinear endogenous regressor: endo_names excludes it, dropped_regressors includes it", {
  d <- make_collin_data()
  d$endo1_dup <- d$endo1
  result <- suppressWarnings(
    .parse_formula(y ~ x1 | endo1 + endo1_dup | z1, data = d)
  )
  # endo1_dup should be dropped (pass 1: intra-endogenous collinearity)
  expect_false("endo1_dup" %in% result$endo_names)
  expect_true("endo1_dup" %in% result$dropped_regressors)
  expect_false("endo1_dup" %in% colnames(result$X))
  # endo1 survives
  expect_true("endo1" %in% result$endo_names)
  expect_true("endo1" %in% colnames(result$X))
})

# ==========================================================================
# 2. Collinear excluded instrument dropped
# ==========================================================================
test_that("collinear excluded instrument: excluded_names excludes it, dropped_instruments includes it", {
  d <- make_collin_data()
  d$z1_dup <- d$z1
  result <- suppressWarnings(
    .parse_formula(y ~ x1 | endo1 | z1 + z1_dup, data = d)
  )
  expect_false("z1_dup" %in% result$excluded_names)
  expect_true("z1_dup" %in% result$dropped_instruments)
  expect_false("z1_dup" %in% colnames(result$Z))
  # z1 survives
  expect_true("z1" %in% result$excluded_names)
  expect_true("z1" %in% colnames(result$Z))
})

# ==========================================================================
# 3. Endogenous variable collinear with instruments → reclassified
# ==========================================================================
test_that("endo collinear with instruments: reclassified as exogenous", {
  d <- make_collin_data()
  # Make endo1 identical to z1 → collinear with instrument set
  d$endo1 <- d$z1
  result <- suppressWarnings(
    .parse_formula(y ~ x1 | endo1 | z1, data = d)
  )
  # endo1 should be reclassified

  expect_true("endo1" %in% result$reclassified_endogenous)
  expect_false("endo1" %in% result$endo_names)
  expect_true("endo1" %in% result$exog_names)
  # K1 reduced, K2 increased
  expect_equal(result$K1, 0L)
  expect_true("endo1" %in% colnames(result$X))
  # Should degenerate to OLS since no endogenous vars remain
  expect_false(result$is_iv)
  expect_null(result$Z)
})

test_that("endo collinear with instruments: warning issued", {
  d <- make_collin_data()
  d$endo1 <- d$z1
  expect_warning(
    .parse_formula(y ~ x1 | endo1 | z1, data = d),
    "Endogenous variable.*collinear with instruments.*Now treated as exogenous: endo1"
  )
})

# ==========================================================================
# 4. Partial reclassification (1 of 2 endo vars reclassified)
# ==========================================================================
test_that("partial reclassification: one endo reclassified, one remains", {
  d <- make_collin_data()
  # Make endo1 collinear with z1, but endo2 stays independent
  d$endo1 <- d$z1
  result <- suppressWarnings(
    .parse_formula(y ~ x1 | endo1 + endo2 | z1 + z2, data = d)
  )
  expect_true("endo1" %in% result$reclassified_endogenous)
  expect_false("endo1" %in% result$endo_names)
  expect_true("endo1" %in% result$exog_names)
  # endo2 remains endogenous
  expect_true("endo2" %in% result$endo_names)
  expect_false("endo2" %in% result$reclassified_endogenous)
  # Model is still IV
  expect_true(result$is_iv)
  expect_equal(result$K1, 1L)
})

# ==========================================================================
# 5. Fields propagate to ivreg2 object
# ==========================================================================
test_that("collinearity fields propagate to ivreg2 object", {
  d <- make_collin_data()
  d$z1_dup <- d$z1
  fit <- suppressWarnings(
    ivreg2(y ~ x1 | endo1 | z1 + z1_dup, data = d)
  )
  expect_true("z1_dup" %in% fit$dropped_instruments)
  expect_false("z1_dup" %in% fit$instruments)
  expect_true("z1" %in% fit$instruments)
  expect_equal(fit$reclassified_endogenous, character(0))
  expect_true("endo1" %in% fit$endogenous)
})

test_that("reclassification fields propagate to ivreg2 object", {
  d <- make_collin_data()
  d$endo1 <- d$z1
  fit <- suppressWarnings(
    ivreg2(y ~ x1 | endo1 + endo2 | z1 + z2, data = d)
  )
  expect_true("endo1" %in% fit$reclassified_endogenous)
  expect_false("endo1" %in% fit$endogenous)
  expect_true("endo2" %in% fit$endogenous)
})

# ==========================================================================
# 6. All endogenous vars reclassified → degenerates to OLS
# ==========================================================================
test_that("all endo reclassified degenerates to OLS", {
  d <- make_collin_data()
  d$endo1 <- d$z1
  result <- suppressWarnings(
    .parse_formula(y ~ x1 | endo1 | z1, data = d)
  )
  expect_false(result$is_iv)
  expect_null(result$Z)
  expect_equal(result$K1, 0L)
  expect_equal(result$endo_names, character(0L))
  expect_equal(result$excluded_names, character(0L))
  # Reclassified var should be in exog and in X
  expect_true("endo1" %in% result$exog_names)
  expect_true("endo1" %in% colnames(result$X))
})

test_that("all endo reclassified: ivreg2 runs OLS", {
  d <- make_collin_data()
  d$endo1 <- d$z1
  fit <- suppressWarnings(
    ivreg2(y ~ x1 | endo1 | z1, data = d)
  )
  # Should be OLS since endo reclassified
  expect_equal(length(fit$endogenous), 0L)
  expect_true("endo1" %in% fit$reclassified_endogenous)
  # Should produce coefficients (i.e., runs successfully)
  expect_true(length(fit$coefficients) > 0L)
})

# ==========================================================================
# 7. Numerical results correct after reclassification
# ==========================================================================
test_that("numerical results match OLS after full reclassification", {
  d <- make_collin_data()
  d$endo1 <- d$z1
  # Full reclassification: endo1 = z1 → becomes exogenous
  fit_iv <- suppressWarnings(
    ivreg2(y ~ x1 | endo1 | z1, data = d)
  )
  # Equivalent OLS
  fit_ols <- ivreg2(y ~ x1 + endo1, data = d)
  expect_equal(coef(fit_iv), coef(fit_ols), tolerance = 1e-10)
})

# ==========================================================================
# 8. return list has reclassified_endogenous field
# ==========================================================================
test_that("parsed_formula return list includes reclassified_endogenous", {
  d <- make_collin_data()
  result <- .parse_formula(y ~ x1 | endo1 | z1, data = d)
  expect_true("reclassified_endogenous" %in% names(result))
  expect_equal(result$reclassified_endogenous, character(0L))
})

# ==========================================================================
# 9. All endo dropped without reclassification clears name vectors
#    (Codex review fix: guard must not require reclassification)
# ==========================================================================
test_that("all endo dropped (zero column) clears endo_names/excluded_names", {
  d <- make_collin_data()
  # Make endo1 a constant zero → collinear with intercept, dropped not reclassified
  d$endo1 <- 0
  result <- suppressWarnings(
    .parse_formula(y ~ x1 | endo1 | z1, data = d)
  )
  # endo1 should be dropped as collinear (with intercept) not reclassified
  # After degeneration to OLS, endo_names must be empty
  expect_equal(result$endo_names, character(0L))
  expect_equal(result$excluded_names, character(0L))
  # ivreg2 object should not claim endogenous variables
  fit <- suppressWarnings(
    ivreg2(y ~ x1 | endo1 | z1, data = d)
  )
  expect_equal(length(fit$endogenous), 0L)
  expect_equal(length(fit$instruments), 0L)
})

# ==========================================================================
# 10. Reclassified-then-dropped var does not leak into exog_names
#     (Codex review fix: filter surviving reclassified by colnames(X))
# ==========================================================================
test_that("reclassified var dropped in pass 3 does not appear in exog_names", {
  d <- make_collin_data()
  # Make endo1 = z1 (triggers reclassification) AND also = x1 (triggers
  # re-drop in pass 3 when rechecking combined [exog, excluded, endo])
  d$endo1 <- d$z1
  d$x1 <- d$z1  # now x1, z1, endo1 are all identical
  result <- suppressWarnings(
    .parse_formula(y ~ x1 | endo1 | z1, data = d)
  )
  # exog_names should not contain variables absent from X
  for (nm in result$exog_names) {
    expect_true(nm %in% colnames(result$X),
                info = paste0("exog_names contains '", nm, "' not in X"))
  }
})
