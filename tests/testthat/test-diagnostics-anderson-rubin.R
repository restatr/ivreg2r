# ============================================================================
# Tests: Anderson-Rubin Test (Ticket E3)
# ============================================================================
#
# Weak-instrument-robust test of H0: all endogenous coefficients = 0.
# Verifies against Stata ivreg2 fixtures.

# --- Helpers ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)

read_diagnostics <- function(path) {
  read.csv(path, check.names = FALSE)
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

sim_cluster_path <- file.path(fixture_dir, "sim_cluster_data.csv")
if (file.exists(sim_cluster_path)) {
  sim_cluster <- read.csv(sim_cluster_path)
}


# ============================================================================
# Structural tests (no fixtures)
# ============================================================================

test_that("anderson_rubin is NULL for OLS models", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$diagnostics$anderson_rubin)
})

test_that("anderson_rubin has all expected fields", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  ar <- fit$diagnostics$anderson_rubin
  expect_type(ar, "list")
  expect_named(ar, c("f_stat", "f_p", "f_df1", "f_df2",
                      "chi2_stat", "chi2_p", "chi2_df"),
               ignore.order = TRUE)
})

test_that("small=TRUE and small=FALSE produce identical AR stats (IID)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = TRUE)
  ar1 <- fit1$diagnostics$anderson_rubin
  ar2 <- fit2$diagnostics$anderson_rubin
  expect_equal(ar1$f_stat, ar2$f_stat)
  expect_equal(ar1$chi2_stat, ar2$chi2_stat)
  expect_equal(ar1$f_p, ar2$f_p)
  expect_equal(ar1$chi2_p, ar2$chi2_p)
})

test_that("small=TRUE and small=FALSE produce identical AR stats (HC1)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1", small = TRUE)
  ar1 <- fit1$diagnostics$anderson_rubin
  ar2 <- fit2$diagnostics$anderson_rubin
  expect_equal(ar1$f_stat, ar2$f_stat)
  expect_equal(ar1$chi2_stat, ar2$chi2_stat)
})

# HC0 and HC1 are numerically identical for AR (the HC1 N/(N-K) correction
# doesn't apply to the AR meat). This test guards against someone accidentally
# introducing a correction inside the AR function.
test_that("HC0 and HC1 produce identical AR stats", {
  skip_if(!file.exists(card_path), "card data not found")
  fit0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC0")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1")
  ar0 <- fit0$diagnostics$anderson_rubin
  ar1 <- fit1$diagnostics$anderson_rubin
  expect_equal(ar0$f_stat, ar1$f_stat)
  expect_equal(ar0$chi2_stat, ar1$chi2_stat)
  expect_equal(ar0$f_p, ar1$f_p)
  expect_equal(ar0$chi2_p, ar1$chi2_p)
})


# ============================================================================
# Helper: run fixture comparison for one spec/VCE combination
# ============================================================================

check_anderson_rubin <- function(fit, fixture_path, label) {
  fixture <- read_diagnostics(fixture_path)
  ar <- fit$diagnostics$anderson_rubin

  test_that(paste("AR F-stat matches Stata", label), {
    expect_equal(ar$f_stat, fixture$arf,
                 tolerance = stata_tol$stat)
  })

  test_that(paste("AR F p-value matches Stata", label), {
    expect_equal(ar$f_p, fixture$arfp,
                 tolerance = stata_tol$pval)
  })

  test_that(paste("AR chi2 matches Stata", label), {
    expect_equal(ar$chi2_stat, fixture$archi2,
                 tolerance = stata_tol$stat)
  })

  test_that(paste("AR chi2 p-value matches Stata", label), {
    expect_equal(ar$chi2_p, fixture$archi2p,
                 tolerance = stata_tol$pval)
  })

  test_that(paste("AR F df1 matches Stata", label), {
    expect_identical(ar$f_df1, as.integer(fixture$ardf))
  })

  test_that(paste("AR F df2 matches Stata", label), {
    expect_identical(ar$f_df2, as.integer(fixture$ardf_r))
  })

  test_that(paste("AR chi2 df matches Stata", label), {
    expect_identical(ar$chi2_df, as.integer(fixture$ardf))
  })
}


# ============================================================================
# card_just_id: lwage ~ exper+expersq+black+south | educ | nearc4
# K1=1, L1=1, just-identified
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC1",  small = FALSE, suffix = "hc1"),
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
    check_anderson_rubin(fit, fixture_file, label)
  }
}


# ============================================================================
# card_overid: lwage ~ exper+expersq+black+south | educ | nearc2+nearc4
# K1=1, L1=2, overidentified
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC1",  small = FALSE, suffix = "hc1"),
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
    check_anderson_rubin(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_cluster: y ~ x1 | endo1 | z1+z2, clusters=~cluster_id
# M=50, includes IID/HC1/CL fixtures
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, clusters = NULL,          suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  clusters = NULL,          suffix = "iid_small"),
  list(vcov = "HC1",  small = FALSE, clusters = NULL,          suffix = "hc1"),
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
    check_anderson_rubin(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_multi_endo: y ~ x1+x2 | endo1+endo2 | z1+z2+z3+z4
# K1=2, L1=4
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC1",  small = FALSE, suffix = "hc1"),
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
    check_anderson_rubin(fit, fixture_file, label)
  }
}


# ============================================================================
# card_just_id_weighted: lwage ~ exper+expersq+black+south | educ | nearc4
# with weights=weight; cluster on smsa66 (M=2)
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, clusters = NULL,       suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  clusters = NULL,       suffix = "iid_small"),
  list(vcov = "HC1",  small = FALSE, clusters = NULL,       suffix = "hc1"),
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
    check_anderson_rubin(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_no_constant: y ~ 0 + x1 | endo1 | z1+z2 (noconstant)
# ============================================================================

sim_noconst_path <- file.path(fixture_dir, "sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC1",  small = FALSE, suffix = "hc1"),
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
    check_anderson_rubin(fit, fixture_file, label)
  }
}
