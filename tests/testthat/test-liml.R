# ============================================================================
# Tests: LIML / Fuller / k-class estimation (Ticket H1)
# ============================================================================

# --- Helper: load Card data and fixtures ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)
card_path <- file.path(fixture_dir, "card_data.csv")
sim_multi_endo_path <- file.path(fixture_dir, "sim_multi_endo_data.csv")

if (file.exists(card_path)) {
  card <- read.csv(card_path)
}
if (file.exists(sim_multi_endo_path)) {
  sim_multi_endo <- read.csv(sim_multi_endo_path)
}

# --- Helper: compare coefficients against Stata fixture ---
compare_coefs <- function(fit, fixture_path, tol = stata_tol$coef) {
  fixture <- read.csv(fixture_path)
  for (i in seq_len(nrow(fixture))) {
    term <- fixture$term[i]
    r_name <- if (term == "_cons") "(Intercept)" else term
    expect_true(r_name %in% names(coef(fit)),
                info = paste("Missing coefficient:", r_name))
    expect_equal(
      unname(coef(fit)[r_name]), fixture$estimate[i],
      tolerance = tol,
      info = paste("Coefficient mismatch:", r_name)
    )
  }
}

# --- Helper: compare VCV against Stata fixture ---
compare_vcov <- function(fit, fixture_path, tol = stata_tol$vcov) {
  fixture <- read.csv(fixture_path)
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
        tolerance = tol,
        info = paste("VCV mismatch:", rn, cn)
      )
    }
  }
}

# --- Helper: compare diagnostics against Stata fixture ---
compare_diagnostics <- function(fit, fixture_path, tol_stat = stata_tol$stat,
                                tol_pval = stata_tol$pval) {
  diag <- read.csv(fixture_path)
  expect_equal(fit$sigma, diag$rmse, tolerance = stata_tol$coef,
               info = "sigma/rmse mismatch")
  expect_equal(fit$rss, diag$rss, tolerance = stata_tol$coef,
               info = "RSS mismatch")
  expect_equal(fit$r.squared, diag$r2, tolerance = stata_tol$coef,
               info = "R-squared mismatch")
  # Model F
  if (!is.na(diag$F_stat)) {
    expect_equal(fit$model_f, diag$F_stat, tolerance = tol_stat,
                 info = "Model F mismatch")
  }
}


# ============================================================================
# 1. LIML overidentified — coefficients, VCV, diagnostics match Stata
# ============================================================================

test_that("LIML overid coefficients match Stata (iid)", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  skip_if(!file.exists(file.path(fixture_dir, "card_liml_overid_coef_iid.csv")),
          "LIML fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")
  compare_coefs(fit, file.path(fixture_dir, "card_liml_overid_coef_iid.csv"))
})

test_that("LIML overid VCV matches Stata (iid)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")
  compare_vcov(fit, file.path(fixture_dir, "card_liml_overid_vcov_iid.csv"))
})

test_that("LIML overid diagnostics match Stata (iid)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")
  diag <- read.csv(file.path(fixture_dir, "card_liml_overid_diagnostics_iid.csv"))

  # Lambda
  expect_equal(fit$lambda, diag$lambda, tolerance = stata_tol$coef,
               info = "lambda mismatch")
  # kclass
  expect_equal(fit$kclass_value, diag$kclass, tolerance = stata_tol$coef,
               info = "kclass mismatch")
  # Sigma, RSS
  compare_diagnostics(fit, file.path(fixture_dir, "card_liml_overid_diagnostics_iid.csv"))
  # Sargan (uses LIML residuals)
  expect_equal(fit$diagnostics$overid$stat, diag$sargan, tolerance = stata_tol$stat,
               info = "Sargan mismatch with LIML residuals")
})

test_that("LIML overid matches Stata (iid, small=TRUE)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml", small = TRUE)
  compare_coefs(fit, file.path(fixture_dir, "card_liml_overid_coef_iid_small.csv"))
  compare_vcov(fit, file.path(fixture_dir, "card_liml_overid_vcov_iid_small.csv"))
  compare_diagnostics(fit, file.path(fixture_dir, "card_liml_overid_diagnostics_iid_small.csv"))
})


# ============================================================================
# 2. LIML exactly-identified — lambda=1, matches 2SLS
# ============================================================================

test_that("LIML exactly-identified: lambda=1, matches 2SLS", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_liml <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, method = "liml")
  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card)

  # lambda must be exactly 1 for just-identified

  expect_equal(fit_liml$lambda, 1.0)
  expect_equal(fit_liml$kclass_value, 1.0)

  # Coefficients must match 2SLS (near machine precision; different
  # computation paths — cross-product vs QR — give ~1e-11 differences)
  expect_equal(coef(fit_liml), coef(fit_2sls), tolerance = 1e-10)
  expect_equal(vcov(fit_liml), vcov(fit_2sls), tolerance = 1e-10)
  expect_equal(fit_liml$sigma, fit_2sls$sigma, tolerance = 1e-10)
})

test_that("LIML justid matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, method = "liml")
  compare_coefs(fit, file.path(fixture_dir, "card_liml_justid_coef_iid.csv"))
  compare_vcov(fit, file.path(fixture_dir, "card_liml_justid_vcov_iid.csv"))

  diag <- read.csv(file.path(fixture_dir, "card_liml_justid_diagnostics_iid.csv"))
  expect_equal(fit$lambda, diag$lambda, tolerance = stata_tol$coef)
  expect_equal(fit$kclass_value, diag$kclass, tolerance = stata_tol$coef)
})


# ============================================================================
# 3. Fuller(1) and Fuller(4)
# ============================================================================

test_that("Fuller(1) matches Stata (iid)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 1)
  compare_coefs(fit, file.path(fixture_dir, "card_fuller1_overid_coef_iid.csv"))
  compare_vcov(fit, file.path(fixture_dir, "card_fuller1_overid_vcov_iid.csv"))

  diag <- read.csv(file.path(fixture_dir, "card_fuller1_overid_diagnostics_iid.csv"))
  expect_equal(fit$lambda, diag$lambda, tolerance = stata_tol$coef,
               info = "lambda mismatch")
  expect_equal(fit$kclass_value, diag$kclass, tolerance = stata_tol$coef,
               info = "kclass mismatch (fuller)")
  expect_equal(fit$fuller_parameter, 1)
  expect_equal(fit$method, "liml")
})

test_that("Fuller(1) matches Stata (iid, small=TRUE)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 1, small = TRUE)
  compare_coefs(fit, file.path(fixture_dir, "card_fuller1_overid_coef_iid_small.csv"))
  compare_vcov(fit, file.path(fixture_dir, "card_fuller1_overid_vcov_iid_small.csv"))
})

test_that("Fuller(4) matches Stata (iid)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 4)
  compare_coefs(fit, file.path(fixture_dir, "card_fuller4_overid_coef_iid.csv"))
  compare_vcov(fit, file.path(fixture_dir, "card_fuller4_overid_vcov_iid.csv"))

  diag <- read.csv(file.path(fixture_dir, "card_fuller4_overid_diagnostics_iid.csv"))
  expect_equal(fit$kclass_value, diag$kclass, tolerance = stata_tol$coef)
  expect_equal(fit$fuller_parameter, 4)
})

test_that("Fuller(4) matches Stata (iid, small=TRUE)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 4, small = TRUE)
  compare_coefs(fit, file.path(fixture_dir, "card_fuller4_overid_coef_iid_small.csv"))
})


# ============================================================================
# 4. kclass(0.5) matches Stata
# ============================================================================

test_that("kclass(0.5) matches Stata (iid)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, kclass = 0.5)
  compare_coefs(fit, file.path(fixture_dir, "card_kclass_half_coef_iid.csv"))
  compare_vcov(fit, file.path(fixture_dir, "card_kclass_half_vcov_iid.csv"))

  expect_equal(fit$method, "kclass")
  expect_equal(fit$kclass_value, 0.5)
  expect_true(is.na(fit$lambda))
})

test_that("kclass(0.5) matches Stata (iid, small=TRUE)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, kclass = 0.5, small = TRUE)
  compare_coefs(fit, file.path(fixture_dir, "card_kclass_half_coef_iid_small.csv"))
  compare_vcov(fit, file.path(fixture_dir, "card_kclass_half_vcov_iid_small.csv"))
})


# ============================================================================
# 5. kclass(1) matches 2SLS exactly
# ============================================================================

test_that("kclass(1) matches 2SLS exactly", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_k1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                   data = card, kclass = 1)
  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card)

  # Near machine-precision agreement (cross-product vs QR paths differ ~1e-11)
  expect_equal(coef(fit_k1), coef(fit_2sls), tolerance = 1e-10)
  expect_equal(vcov(fit_k1), vcov(fit_2sls), tolerance = 1e-10)
  expect_equal(fit_k1$sigma, fit_2sls$sigma, tolerance = 1e-10)
  expect_equal(fit_k1$rss, fit_2sls$rss, tolerance = 1e-10)
})

test_that("kclass(1) matches Stata kclass(1) fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, kclass = 1)
  compare_coefs(fit, file.path(fixture_dir, "card_kclass_1_coef_iid.csv"))
  compare_vcov(fit, file.path(fixture_dir, "card_kclass_1_vcov_iid.csv"))
})


# ============================================================================
# 6. kclass(0) matches OLS exactly
# ============================================================================

test_that("kclass(0) matches OLS exactly", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_k0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                   data = card, kclass = 0)
  fit_ols <- ivreg2(lwage ~ exper + expersq + black + south + educ, data = card)

  # OLS with all regressors (including educ)
  expect_equal(unname(coef(fit_k0)), unname(coef(fit_ols)), tolerance = 1e-10)
  expect_equal(fit_k0$sigma, fit_ols$sigma, tolerance = 1e-10)
})


# ============================================================================
# 7. Multi-endogenous LIML (K1 > 1)
# ============================================================================

test_that("Multi-endogenous LIML matches Stata (iid)", {
  skip_if(!file.exists(sim_multi_endo_path), "Simulated data not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, method = "liml")
  compare_coefs(fit, file.path(fixture_dir, "sim_multi_endo_liml_coef_iid.csv"))
  compare_vcov(fit, file.path(fixture_dir, "sim_multi_endo_liml_vcov_iid.csv"))

  diag <- read.csv(file.path(fixture_dir, "sim_multi_endo_liml_diagnostics_iid.csv"))
  expect_equal(fit$lambda, diag$lambda, tolerance = stata_tol$coef,
               info = "multi-endo lambda mismatch")
  expect_equal(fit$kclass_value, diag$kclass, tolerance = stata_tol$coef,
               info = "multi-endo kclass mismatch")
})

test_that("Multi-endogenous LIML matches Stata (iid, small=TRUE)", {
  skip_if(!file.exists(sim_multi_endo_path), "Simulated data not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, method = "liml", small = TRUE)
  compare_coefs(fit, file.path(fixture_dir, "sim_multi_endo_liml_coef_iid_small.csv"))
  compare_vcov(fit, file.path(fixture_dir, "sim_multi_endo_liml_vcov_iid_small.csv"))
})

test_that("Multi-endogenous Fuller(1) matches Stata (iid)", {
  skip_if(!file.exists(sim_multi_endo_path), "Simulated data not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, fuller = 1)
  compare_coefs(fit, file.path(fixture_dir, "sim_multi_endo_fuller1_coef_iid.csv"))
  compare_vcov(fit, file.path(fixture_dir, "sim_multi_endo_fuller1_vcov_iid.csv"))

  diag <- read.csv(file.path(fixture_dir, "sim_multi_endo_fuller1_diagnostics_iid.csv"))
  expect_equal(fit$kclass_value, diag$kclass, tolerance = stata_tol$coef)
  expect_equal(fit$lambda, diag$lambda, tolerance = stata_tol$coef)
})


# ============================================================================
# 8. Weighted LIML
# ============================================================================

test_that("Weighted LIML matches Stata (iid)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, weights = weight, method = "liml")
  compare_coefs(fit, file.path(fixture_dir, "card_liml_weighted_coef_iid.csv"))
  compare_vcov(fit, file.path(fixture_dir, "card_liml_weighted_vcov_iid.csv"))

  diag <- read.csv(file.path(fixture_dir, "card_liml_weighted_diagnostics_iid.csv"))
  expect_equal(fit$lambda, diag$lambda, tolerance = stata_tol$coef,
               info = "weighted lambda mismatch")
  expect_equal(fit$sigma, diag$rmse, tolerance = stata_tol$coef,
               info = "weighted sigma mismatch")
})

test_that("Weighted LIML matches Stata (iid, small=TRUE)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, weights = weight, method = "liml", small = TRUE)
  compare_coefs(fit, file.path(fixture_dir, "card_liml_weighted_coef_iid_small.csv"))
  compare_vcov(fit, file.path(fixture_dir, "card_liml_weighted_vcov_iid_small.csv"))
})


# ============================================================================
# 9. Return object fields
# ============================================================================

test_that("LIML return object has correct fields", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")

  expect_equal(fit$method, "liml")
  expect_true(is.numeric(fit$lambda))
  expect_true(!is.na(fit$lambda))
  expect_true(fit$lambda > 1)  # overidentified: lambda > 1
  expect_equal(fit$kclass_value, fit$lambda)  # pure LIML: k = lambda
  expect_equal(fit$fuller_parameter, 0)
})

test_that("Fuller return object has correct fields", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 1)

  expect_equal(fit$method, "liml")
  expect_true(!is.na(fit$lambda))
  expect_equal(fit$fuller_parameter, 1)
  # k = lambda - fuller/(N-L)
  expected_k <- fit$lambda - 1 / (nobs(fit) - 7)  # L=7 for this model
  expect_equal(fit$kclass_value, expected_k, tolerance = 1e-12)
})

test_that("kclass return object has correct fields", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, kclass = 0.5)

  expect_equal(fit$method, "kclass")
  expect_true(is.na(fit$lambda))  # no eigenvalue computed
  expect_equal(fit$kclass_value, 0.5)
  expect_equal(fit$fuller_parameter, 0)
})

test_that("2SLS return object has correct method field", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)

  expect_equal(fit$method, "2sls")
  expect_true(is.na(fit$lambda))
  expect_true(is.na(fit$kclass_value))
  expect_equal(fit$fuller_parameter, 0)
})

test_that("OLS return object has correct method field", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)

  expect_equal(fit$method, "ols")
  expect_true(is.na(fit$lambda))
  expect_true(is.na(fit$kclass_value))
  expect_equal(fit$fuller_parameter, 0)
})


# ============================================================================
# 10. Error tests
# ============================================================================

test_that("LIML + robust errors with informative message", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "liml", vcov = "HC1"),
    "not yet implemented"
  )
})

test_that("LIML + cluster errors with informative message", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "liml", clusters = ~smsa66),
    "not yet implemented"
  )
})

test_that("kclass + cluster errors with informative message", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, kclass = 0.5, clusters = ~smsa66),
    "not yet implemented"
  )
})

test_that("LIML + kclass is rejected", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "liml", kclass = 0.5),
    "Cannot specify `kclass` with"
  )
})

test_that("fuller + kclass is rejected", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, fuller = 1, kclass = 0.5),
    "Cannot specify both"
  )
})

test_that("fuller < 0 is rejected", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, fuller = -1),
    "non-negative"
  )
})

test_that("kclass < 0 is rejected", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, kclass = -0.5),
    "non-negative"
  )
})

test_that("LIML without IV is rejected", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, method = "liml"),
    "requires an IV model"
  )
})

test_that("kclass without IV is rejected", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, kclass = 0.5),
    "requires an IV model"
  )
})

test_that("invalid method is rejected", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, method = "gmm"),
    'must be one of'
  )
})

test_that("fuller auto-promotes method to liml", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  # fuller with default method="2sls" should auto-promote
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 1)
  expect_equal(fit$method, "liml")
})

test_that("kclass auto-promotes method to kclass", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  # kclass with default method="2sls" should auto-promote
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, kclass = 0.5)
  expect_equal(fit$method, "kclass")
})


# ============================================================================
# 11. Display and broom methods
# ============================================================================

test_that("print.ivreg2 shows LIML estimation type", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")
  output <- capture.output(print(fit))
  expect_true(any(grepl("LIML Estimation", output)))
})

test_that("print.ivreg2 shows Fuller LIML estimation type", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 1)
  output <- capture.output(print(fit))
  expect_true(any(grepl("Fuller LIML Estimation", output)))
})

test_that("print.ivreg2 shows k-class estimation type", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, kclass = 0.5)
  output <- capture.output(print(fit))
  expect_true(any(grepl("k-class Estimation", output)))
})

test_that("summary.ivreg2 shows lambda and kclass", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")
  output <- capture.output(print(summary(fit)))
  expect_true(any(grepl("lambda:", output)))
  expect_true(any(grepl("kclass:", output)))
})

test_that("glance includes LIML columns", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")
  g <- glance(fit)

  expect_true("method" %in% names(g))
  expect_true("lambda" %in% names(g))
  expect_true("kclass_value" %in% names(g))
  expect_true("fuller_parameter" %in% names(g))
  expect_equal(g$method, "liml")
  expect_false(is.na(g$lambda))
  expect_equal(g$fuller_parameter, 0)
})

test_that("glance includes Fuller columns", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 1)
  g <- glance(fit)

  expect_equal(g$method, "liml")
  expect_equal(g$fuller_parameter, 1)
  expect_false(is.na(g$kclass_value))
})


# ============================================================================
# 12. Existing functionality unchanged
# ============================================================================

test_that("2SLS results unchanged after LIML addition", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card)
  compare_coefs(fit, file.path(fixture_dir, "card_overid_coef_iid.csv"))
  compare_vcov(fit, file.path(fixture_dir, "card_overid_vcov_iid.csv"))
})

test_that("OLS results unchanged after LIML addition", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$method, "ols")
  expect_true(is.na(fit$lambda))
})
