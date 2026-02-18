# ============================================================================
# Tests: Summary statistics — R², adj R², sigma, RSS (Ticket F2)
# ============================================================================
#
# Fixture-based tests verify r2, r2_a, r2u, r2c, rss, rmse, sigmasq against
# Stata's ivreg2 output. Structural tests verify invariance properties.

# --- Helpers ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)

read_diagnostics <- function(path) {
  read.csv(path)
}

# --- Load datasets ---
card_path <- file.path(fixture_dir, "card_data.csv")
if (file.exists(card_path)) {

  card <- read.csv(card_path)
}

sim_multi_path <- file.path(fixture_dir, "sim_multi_endo_data.csv")
if (file.exists(sim_multi_path)) {
  sim_multi <- read.csv(sim_multi_path)
}

sim_noconst_path <- file.path(fixture_dir, "sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

sim_cluster_path <- file.path(fixture_dir, "sim_cluster_data.csv")
if (file.exists(sim_cluster_path)) {
  sim_cluster <- read.csv(sim_cluster_path)
}


# ============================================================================
# Structural tests (no fixtures)
# ============================================================================

test_that("OLS r.squared matches summary(lm(...))$r.squared", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$r.squared, summary(lm_fit)$r.squared, tolerance = 1e-10)
})

test_that("OLS adj.r.squared matches summary(lm(...))$adj.r.squared", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$adj.r.squared, summary(lm_fit)$adj.r.squared,
               tolerance = 1e-10)
})

test_that("R-squared is invariant to small", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = TRUE)
  expect_identical(fit1$r.squared, fit2$r.squared)
  expect_identical(fit1$adj.r.squared, fit2$adj.r.squared)
  expect_identical(fit1$r2u, fit2$r2u)
  expect_identical(fit1$r2c, fit2$r2c)
})

test_that("R-squared is invariant to VCE type", {
  skip_if(!file.exists(card_path), "card data not found")
  fit_iid <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card, vcov = "iid")
  fit_hc  <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card, vcov = "HC1", small = TRUE)
  expect_identical(fit_iid$r.squared, fit_hc$r.squared)
  expect_identical(fit_iid$adj.r.squared, fit_hc$adj.r.squared)
  expect_identical(fit_iid$r2u, fit_hc$r2u)
  expect_identical(fit_iid$r2c, fit_hc$r2c)
})

test_that("Sigma changes with small: sigma_small^2 * (N-K) == sigma_large^2 * N", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = TRUE)
  expect_equal(fit2$sigma^2 * (fit2$nobs - fit2$rank),
               fit1$sigma^2 * fit1$nobs,
               tolerance = 1e-12)
})

test_that("Negative R-squared for IV is allowed (no error/warning)", {
  skip_if(!file.exists(card_path), "card data not found")
  expect_no_warning(
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card)
  )
  expect_true(fit$r.squared < 0)
})

test_that("mss = tss - rss identity holds", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  tss <- fit$rss / (1 - fit$r.squared)
  expect_equal(fit$mss, tss - fit$rss, tolerance = 1e-12)
})

test_that("With intercept: r2 == r2c; without: r2 == r2u", {
  # With intercept
  fit_int <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_identical(fit_int$r.squared, fit_int$r2c)

  # Without intercept
  fit_noint <- ivreg2(mpg ~ 0 + wt + hp, data = mtcars)
  expect_identical(fit_noint$r.squared, fit_noint$r2u)
})


# ============================================================================
# Helper: run summary stats fixture comparison for one spec/VCE combination
# ============================================================================

check_summary_stats <- function(fit, fixture_path, label) {
  fixture <- read_diagnostics(fixture_path)

  test_that(paste("r.squared matches Stata", label), {
    expect_equal(fit$r.squared, fixture$r2,
                 tolerance = stata_tol$coef)
  })

  test_that(paste("adj.r.squared matches Stata", label), {
    expect_equal(fit$adj.r.squared, fixture$r2_a,
                 tolerance = stata_tol$coef)
  })

  test_that(paste("r2u matches Stata", label), {
    expect_equal(fit$r2u, fixture$r2u,
                 tolerance = stata_tol$coef)
  })

  test_that(paste("r2c matches Stata", label), {
    expect_equal(fit$r2c, fixture$r2c,
                 tolerance = stata_tol$coef)
  })

  test_that(paste("rss matches Stata", label), {
    expect_equal(fit$rss, fixture$rss,
                 tolerance = stata_tol$coef)
  })

  test_that(paste("sigma matches Stata rmse", label), {
    expect_equal(fit$sigma, fixture$rmse,
                 tolerance = stata_tol$coef)
  })

  test_that(paste("sigma^2 matches Stata sigmasq", label), {
    expect_equal(fit$sigma^2, fixture$sigmasq,
                 tolerance = stata_tol$coef)
  })
}


# ============================================================================
# card_just_id: lwage ~ exper+expersq+black+south | educ | nearc4
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1"),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("card_just_id_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_just_id", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card, vcov = vce_combo$vcov, small = vce_combo$small)
    check_summary_stats(fit, fixture_file, label)
  }
}


# ============================================================================
# card_overid: lwage ~ exper+expersq+black+south | educ | nearc2+nearc4
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1"),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("card_overid_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_overid", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                  data = card, vcov = vce_combo$vcov, small = vce_combo$small)
    check_summary_stats(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_multi_endo: y ~ x1+x2 | endo1+endo2 | z1+z2+z3+z4
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1"),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_multi_endo_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_multi_endo", vce_combo$suffix)

  if (file.exists(sim_multi_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                  data = sim_multi, vcov = vce_combo$vcov,
                  small = vce_combo$small)
    check_summary_stats(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_no_constant: y ~ 0+x1 | endo1 | z1+z2
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1"),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_no_constant_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_no_constant", vce_combo$suffix)

  if (file.exists(sim_noconst_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2,
                  data = sim_noconst, vcov = vce_combo$vcov,
                  small = vce_combo$small)
    check_summary_stats(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_cluster: y ~ x1 | endo1 | z1+z2, clusters=~cluster_id
# Includes IID, HC, and CL variants
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, clusters = NULL,          suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  clusters = NULL,          suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, clusters = NULL,          suffix = "hc1"),
  list(vcov = "HC1",  small = TRUE,  clusters = NULL,          suffix = "hc1_small"),
  list(vcov = "iid",  small = FALSE, clusters = ~cluster_id,   suffix = "cl"),
  list(vcov = "iid",  small = TRUE,  clusters = ~cluster_id,   suffix = "cl_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_cluster_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_cluster", vce_combo$suffix)

  if (file.exists(sim_cluster_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                  data = sim_cluster, vcov = vce_combo$vcov,
                  small = vce_combo$small, clusters = vce_combo$clusters)
    check_summary_stats(fit, fixture_file, label)
  }
}


# ============================================================================
# card_just_id_weighted: lwage ~ exper+expersq+black+south | educ | nearc4
# with weights=weight; cluster on smsa66 (M=2)
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, clusters = NULL,       suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  clusters = NULL,       suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, clusters = NULL,       suffix = "hc1"),
  list(vcov = "HC1",  small = TRUE,  clusters = NULL,       suffix = "hc1_small"),
  list(vcov = "iid",  small = FALSE, clusters = ~smsa66,    suffix = "cl"),
  list(vcov = "iid",  small = TRUE,  clusters = ~smsa66,    suffix = "cl_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("card_just_id_weighted_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_just_id_weighted", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card, weights = weight,
                  vcov = vce_combo$vcov, small = vce_combo$small,
                  clusters = vce_combo$clusters)
    check_summary_stats(fit, fixture_file, label)
  }
}
