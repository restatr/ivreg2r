# ============================================================================
# Tests: ivreg2() 2SLS estimation (Ticket B1)
# ============================================================================

# --- Helper: load Card data and fixtures ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)
card_path <- file.path(fixture_dir, "card_data.csv")
coef_path <- file.path(fixture_dir, "card_just_id_coef_iid.csv")
vcov_path <- file.path(fixture_dir, "card_just_id_vcov_iid.csv")

if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

# ============================================================================
# Coefficients match Stata card_just_id fixture (iid, small=FALSE)
# ============================================================================

test_that("2SLS coefficients match Stata card_just_id fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  skip_if(!file.exists(coef_path), "Coefficient fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  fixture <- read.csv(coef_path)

  for (i in seq_len(nrow(fixture))) {
    term <- fixture$term[i]
    r_name <- if (term == "_cons") "(Intercept)" else term
    expect_true(r_name %in% names(coef(fit)),
                info = paste("Missing coefficient:", r_name))
    expect_equal(
      unname(coef(fit)[r_name]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coefficient mismatch:", r_name)
    )
  }
})

# ============================================================================
# VCV matches Stata fixture (iid, small=FALSE)
# ============================================================================

test_that("2SLS vcov matches Stata card_just_id fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  fixture <- read.csv(vcov_path)

  # Map Stata names to R names
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  # Extract fixture VCV matrix
  vcov_cols <- grep("^vcov_", names(fixture), value = TRUE)
  V_stata <- as.matrix(fixture[, vcov_cols])
  rownames(V_stata) <- r_names
  col_stata <- sub("^vcov_", "", vcov_cols)
  colnames(V_stata) <- ifelse(col_stata == "_cons", "(Intercept)", col_stata)

  # Compare each element in shared order
  shared <- intersect(r_names, rownames(fit$vcov))
  for (rn in shared) {
    for (cn in shared) {
      expect_equal(
        fit$vcov[rn, cn], V_stata[rn, cn],
        tolerance = stata_tol$vcov,
        info = paste("VCV mismatch:", rn, cn)
      )
    }
  }
})

# ============================================================================
# Object structure
# ============================================================================

test_that("2SLS object has correct structure", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)

  # Coefficients
  expect_length(coef(fit), 6L)
  expect_named(coef(fit),
               c("(Intercept)", "exper", "expersq", "black", "south", "educ"),
               ignore.order = TRUE)

  # Residuals and fitted values
  expect_length(fit$residuals, 3010L)
  expect_length(fit$fitted.values, 3010L)

  # Rank and nobs
  expect_identical(fit$rank, 6L)
  expect_identical(fit$nobs, 3010L)

  # RSS is positive
  expect_gt(fit$rss, 0)

  # Class
  expect_s3_class(fit, "ivreg2")
})

# ============================================================================
# Residuals use original X (not X_hat)
# ============================================================================

test_that("fitted + residuals == y (original X used)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, y = TRUE)
  expect_equal(fit$fitted.values + fit$residuals, fit$y,
               tolerance = .Machine$double.eps^0.5)
})

# ============================================================================
# Cross-check against ivreg package
# ============================================================================

test_that("2SLS coefficients match ivreg::ivreg()", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  skip_if_not_installed("ivreg")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  iv_fit <- ivreg::ivreg(
    lwage ~ educ + exper + expersq + black + south | nearc4 + exper + expersq + black + south,
    data = card
  )

  # Match coefficients in shared order
  shared <- intersect(names(coef(fit)), names(coef(iv_fit)))
  expect_equal(coef(fit)[shared], coef(iv_fit)[shared],
               tolerance = .Machine$double.eps^0.5)
})

test_that("2SLS residuals match ivreg::ivreg()", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  skip_if_not_installed("ivreg")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  iv_fit <- ivreg::ivreg(
    lwage ~ educ + exper + expersq + black + south | nearc4 + exper + expersq + black + south,
    data = card
  )
  expect_equal(unname(fit$residuals), unname(residuals(iv_fit)),
               tolerance = .Machine$double.eps^0.5)
})

# ============================================================================
# Print method shows "2SLS Estimation"
# ============================================================================

test_that("print method shows 2SLS Estimation for IV model", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  expect_output(print(fit), "2SLS Estimation")
})

# ============================================================================
# endogenous and instruments are populated
# ============================================================================

test_that("endogenous and instruments are populated for IV", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  expect_identical(fit$endogenous, "educ")
  expect_identical(fit$instruments, "nearc4")
})

# ============================================================================
# model/x/y flags work for IV
# ============================================================================

test_that("x=TRUE stores both X and Z for IV", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, x = TRUE)
  expect_false(is.null(fit$x))
  expect_true(is.matrix(fit$x$X))
  expect_true(is.matrix(fit$x$Z))
  expect_equal(ncol(fit$x$X), 6L)  # intercept + exper + expersq + black + south + educ
  expect_equal(ncol(fit$x$Z), 6L)  # intercept + exper + expersq + black + south + nearc4
})

test_that("y=TRUE stores response for IV", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, y = TRUE)
  expect_false(is.null(fit$y))
  expect_length(fit$y, 3010L)
})

test_that("model=TRUE stores model frame for IV", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, model = TRUE)
  expect_false(is.null(fit$model))
  expect_s3_class(fit$model, "data.frame")
})

# ============================================================================
# OLS still works (regression guard)
# ============================================================================

test_that("OLS estimation still works after 2SLS implementation", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(coef(fit), coef(lm_fit), tolerance = .Machine$double.eps^0.5)
  expect_output(print(fit), "OLS Estimation")
  expect_identical(fit$endogenous, character(0))
  expect_identical(fit$instruments, character(0))
})

# ============================================================================
# small=TRUE changes sigma denominator
# ============================================================================

# ============================================================================
# Rank-deficient X_hat errors
# ============================================================================

test_that("rank-deficient X_hat raises informative error", {
  set.seed(123)
  n <- 50
  x <- rnorm(n)
  z <- rnorm(n)
  # Construct e so that its projection onto Z = [1, x, z] lies in span(1, x).
  # This means z has exactly zero predictive power for e after partialling out x.
  Z_mat <- cbind(1, x, z)
  noise <- rnorm(n)
  # Project noise out of Z's column space
  noise_orth <- noise - Z_mat %*% solve(crossprod(Z_mat), crossprod(Z_mat, noise))
  e <- 1 + 0.5 * x + noise_orth  # X is full rank, but X_hat is not
  y <- 1 + x + e + rnorm(n)
  d <- data.frame(y = y, x = x, e = e, z = z)
  expect_error(
    ivreg2(y ~ x | e | z, data = d),
    "rank-deficient"
  )
})

# ============================================================================
# small=TRUE changes sigma denominator
# ============================================================================

test_that("2SLS small=TRUE uses N-K denominator", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_small <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                      data = card, small = TRUE)
  fit_large <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                      data = card, small = FALSE)

  N <- fit_large$nobs
  K <- length(coef(fit_large))
  expect_equal(fit_small$sigma^2 * (N - K), fit_large$sigma^2 * N,
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# sim_no_constant: y ~ 0 + x1 | endo1 | z1+z2 (noconstant)
# ============================================================================

sim_noconst_path <- file.path(fixture_dir, "sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

test_that("2SLS coefficients match Stata sim_no_constant iid fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  noconst_coef_path <- file.path(fixture_dir, "sim_no_constant_coef_iid.csv")
  skip_if(!file.exists(noconst_coef_path), "Coefficient fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst)
  fixture <- read.csv(noconst_coef_path)

  for (i in seq_len(nrow(fixture))) {
    term <- fixture$term[i]
    r_name <- if (term == "_cons") "(Intercept)" else term
    expect_true(r_name %in% names(coef(fit)),
                info = paste("Missing coefficient:", r_name))
    expect_equal(
      unname(coef(fit)[r_name]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coefficient mismatch:", r_name)
    )
  }
})

test_that("no-constant model has no intercept in coef()", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst)
  expect_false("(Intercept)" %in% names(coef(fit)))
})

test_that("2SLS vcov matches Stata sim_no_constant iid fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  noconst_vcov_path <- file.path(fixture_dir, "sim_no_constant_vcov_iid.csv")
  skip_if(!file.exists(noconst_vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst)
  fixture <- read.csv(noconst_vcov_path)

  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)
  vcov_cols <- grep("^vcov_", names(fixture), value = TRUE)
  V_stata <- as.matrix(fixture[, vcov_cols])
  rownames(V_stata) <- r_names
  col_stata <- sub("^vcov_", "", vcov_cols)
  colnames(V_stata) <- ifelse(col_stata == "_cons", "(Intercept)", col_stata)

  shared <- intersect(r_names, rownames(fit$vcov))
  for (rn in shared) {
    for (cn in shared) {
      expect_equal(
        fit$vcov[rn, cn], V_stata[rn, cn],
        tolerance = stata_tol$vcov,
        info = paste("VCV mismatch:", rn, cn)
      )
    }
  }
})
