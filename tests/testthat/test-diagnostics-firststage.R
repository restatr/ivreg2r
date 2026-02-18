# ============================================================================
# Tests: First-Stage Diagnostics (Ticket E1)
# ============================================================================
#
# First-stage F-test, partial R², Shea partial R², and RMSE.
# Verifies against Stata ivreg2 fixtures.

# --- Helpers ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)

read_firststage <- function(path) {
  read.csv(path, check.names = FALSE)
}

get_fs_value <- function(fixture, stat, endo_name) {
  as.numeric(fixture[fixture$statistic == stat, endo_name])
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

test_that("first_stage is NULL for OLS models", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$first_stage)
})

test_that("first_stage is named list with names == endogenous", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  expect_type(fit$first_stage, "list")
  expect_named(fit$first_stage, "educ")
})

test_that("first_stage has all expected fields", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  expected_fields <- c("f_stat", "f_p", "f_df1", "f_df2",
                       "partial_r2", "shea_partial_r2", "rmse",
                       "coefficients", "residuals", "fitted.values",
                       "sw_f", "sw_f_p", "sw_f_df1", "sw_f_df2",
                       "sw_chi2", "sw_chi2_p", "sw_partial_r2",
                       "ap_f", "ap_f_p", "ap_f_df1", "ap_f_df2",
                       "ap_chi2", "ap_chi2_p", "ap_partial_r2")
  expect_named(fit$first_stage$educ, expected_fields, ignore.order = TRUE)
})

test_that("first_stage coefficients/residuals/fitted have correct lengths", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  fs <- fit$first_stage$educ
  expect_length(fs$coefficients, 6L)   # L = 6 instruments
  expect_length(fs$residuals, fit$nobs)
  expect_length(fs$fitted.values, fit$nobs)
})

test_that("multi-endo first_stage has correct names", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi)
  expect_named(fit$first_stage, c("endo1", "endo2"))
})

test_that("K1=1: shea_partial_r2 approx equals partial_r2", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  fs <- fit$first_stage$educ
  expect_equal(fs$shea_partial_r2, fs$partial_r2, tolerance = 1e-4)
})

test_that("small does not change first-stage F-stat (IID)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = TRUE)
  expect_equal(fit1$first_stage$educ$f_stat,
               fit2$first_stage$educ$f_stat)
})

test_that("small does not change first-stage F-stat (HC1)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1", small = TRUE)
  expect_equal(fit1$first_stage$educ$f_stat,
               fit2$first_stage$educ$f_stat)
})

test_that("HC0 and HC1 produce same first-stage F-stat", {
  skip_if(!file.exists(card_path), "card data not found")
  fit0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC0")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1")
  expect_equal(fit0$first_stage$educ$f_stat,
               fit1$first_stage$educ$f_stat)
})

test_that("K1=1: SW/AP F equals standard F", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  fs <- fit$first_stage$educ
  expect_equal(fs$sw_f, fs$f_stat)
  expect_equal(fs$ap_f, fs$f_stat)
  expect_equal(fs$sw_f_p, fs$f_p)
  expect_equal(fs$ap_f_p, fs$f_p)
})

test_that("K1=1: SW/AP partial R2 equals partial R2", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  fs <- fit$first_stage$educ
  expect_equal(fs$sw_partial_r2, fs$partial_r2)
  expect_equal(fs$ap_partial_r2, fs$partial_r2)
})

test_that("K1=1: SW/AP df1 equals L1", {
  skip_if(!file.exists(card_path), "card data not found")
  # L1=1 (just identified)
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card)
  expect_identical(fit1$first_stage$educ$sw_f_df1, 1L)
  expect_identical(fit1$first_stage$educ$ap_f_df1, 1L)
  # L1=2 (overidentified)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                 data = card)
  expect_identical(fit2$first_stage$educ$sw_f_df1, 2L)
  expect_identical(fit2$first_stage$educ$ap_f_df1, 2L)
})

test_that("K1>1: SW/AP df1 equals L1 - K1 + 1", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi)
  # L1=4, K1=2 => Fdf1 = 3
  expect_identical(fit$first_stage$endo1$sw_f_df1, 3L)
  expect_identical(fit$first_stage$endo1$ap_f_df1, 3L)
  expect_identical(fit$first_stage$endo2$sw_f_df1, 3L)
  expect_identical(fit$first_stage$endo2$ap_f_df1, 3L)
})

test_that("K1=1: SW/AP identity holds with cluster VCE", {
  skip_if(!file.exists(sim_cluster_path), "sim_cluster data not found")
  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id)
  fs <- fit$first_stage$endo1
  expect_equal(fs$sw_f, fs$f_stat)
  expect_equal(fs$ap_f, fs$f_stat)
  expect_equal(fs$sw_partial_r2, fs$partial_r2)
  expect_equal(fs$ap_partial_r2, fs$partial_r2)
})

test_that("K1=1: SW/AP identity holds with weights", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight)
  fs <- fit$first_stage$educ
  expect_equal(fs$sw_f, fs$f_stat)
  expect_equal(fs$ap_f, fs$f_stat)
  expect_equal(fs$sw_partial_r2, fs$partial_r2)
  expect_equal(fs$ap_partial_r2, fs$partial_r2)
})

test_that("SW/AP partial R2 invariant to VCE type (multi-endo)", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit_iid <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                    data = sim_multi, vcov = "iid")
  fit_hc1 <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                    data = sim_multi, vcov = "HC1")
  for (endo in c("endo1", "endo2")) {
    expect_equal(fit_iid$first_stage[[endo]]$sw_partial_r2,
                 fit_hc1$first_stage[[endo]]$sw_partial_r2)
    expect_equal(fit_iid$first_stage[[endo]]$ap_partial_r2,
                 fit_hc1$first_stage[[endo]]$ap_partial_r2)
  }
})


# ============================================================================
# Helper: run fixture comparison for one spec/VCE combination
# ============================================================================

check_firststage <- function(fit, fixture_path, endo_name, label) {
  fixture <- read_firststage(fixture_path)
  fs <- fit$first_stage[[endo_name]]

  test_that(paste("RMSE matches Stata", label), {
    expect_equal(fs$rmse,
                 get_fs_value(fixture, "rmse", endo_name),
                 tolerance = stata_tol$coef)
  })

  test_that(paste("partial R2 matches Stata", label), {
    expect_equal(fs$partial_r2,
                 get_fs_value(fixture, "pr2", endo_name),
                 tolerance = stata_tol$stat)
  })

  test_that(paste("Shea partial R2 matches Stata", label), {
    expect_equal(fs$shea_partial_r2,
                 get_fs_value(fixture, "sheapr2", endo_name),
                 tolerance = stata_tol$stat)
  })

  test_that(paste("F-stat matches Stata", label), {
    expect_equal(fs$f_stat,
                 get_fs_value(fixture, "F", endo_name),
                 tolerance = stata_tol$stat)
  })

  test_that(paste("F df matches Stata", label), {
    expect_identical(fs$f_df1,
                     as.integer(get_fs_value(fixture, "df", endo_name)))
    expect_identical(fs$f_df2,
                     as.integer(get_fs_value(fixture, "df_r", endo_name)))
  })

  test_that(paste("F p-value matches Stata", label), {
    expect_equal(fs$f_p,
                 get_fs_value(fixture, "pvalue", endo_name),
                 tolerance = stata_tol$pval)
  })

  # --- SW diagnostics ---
  test_that(paste("SW F matches Stata", label), {
    expect_equal(fs$sw_f,
                 get_fs_value(fixture, "SWF", endo_name),
                 tolerance = stata_tol$stat)
  })

  test_that(paste("SW F df matches Stata", label), {
    expect_identical(fs$sw_f_df1,
                     as.integer(get_fs_value(fixture, "SWFdf1", endo_name)))
    expect_identical(fs$sw_f_df2,
                     as.integer(get_fs_value(fixture, "SWFdf2", endo_name)))
  })

  test_that(paste("SW F p-value matches Stata", label), {
    expect_equal(fs$sw_f_p,
                 get_fs_value(fixture, "SWFp", endo_name),
                 tolerance = stata_tol$pval)
  })

  test_that(paste("SW chi2 matches Stata", label), {
    expect_equal(fs$sw_chi2,
                 get_fs_value(fixture, "SWchi2", endo_name),
                 tolerance = stata_tol$stat)
  })

  test_that(paste("SW chi2 p-value matches Stata", label), {
    expect_equal(fs$sw_chi2_p,
                 get_fs_value(fixture, "SWchi2p", endo_name),
                 tolerance = stata_tol$pval)
  })

  test_that(paste("SW partial R2 matches Stata", label), {
    expect_equal(fs$sw_partial_r2,
                 get_fs_value(fixture, "SWr2", endo_name),
                 tolerance = stata_tol$stat)
  })

  # --- AP diagnostics ---
  test_that(paste("AP F matches Stata", label), {
    expect_equal(fs$ap_f,
                 get_fs_value(fixture, "APF", endo_name),
                 tolerance = stata_tol$stat)
  })

  test_that(paste("AP F df matches Stata", label), {
    expect_identical(fs$ap_f_df1,
                     as.integer(get_fs_value(fixture, "APFdf1", endo_name)))
    expect_identical(fs$ap_f_df2,
                     as.integer(get_fs_value(fixture, "APFdf2", endo_name)))
  })

  test_that(paste("AP F p-value matches Stata", label), {
    expect_equal(fs$ap_f_p,
                 get_fs_value(fixture, "APFp", endo_name),
                 tolerance = stata_tol$pval)
  })

  test_that(paste("AP chi2 matches Stata", label), {
    expect_equal(fs$ap_chi2,
                 get_fs_value(fixture, "APchi2", endo_name),
                 tolerance = stata_tol$stat)
  })

  test_that(paste("AP chi2 p-value matches Stata", label), {
    expect_equal(fs$ap_chi2_p,
                 get_fs_value(fixture, "APchi2p", endo_name),
                 tolerance = stata_tol$pval)
  })

  test_that(paste("AP partial R2 matches Stata", label), {
    expect_equal(fs$ap_partial_r2,
                 get_fs_value(fixture, "APr2", endo_name),
                 tolerance = stata_tol$stat)
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
    paste0("card_just_id_firststage_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_just_id", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card, vcov = vce_combo$vcov, small = vce_combo$small)
    check_firststage(fit, fixture_file, "educ", label)
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
    paste0("card_overid_firststage_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_overid", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                  data = card, vcov = vce_combo$vcov, small = vce_combo$small)
    check_firststage(fit, fixture_file, "educ", label)
  }
}


# ============================================================================
# sim_multi_endo: y ~ x1+x2 | endo1+endo2 | z1+z2+z3+z4
# K1=2, L1=4, shea != pr2
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC1",  small = FALSE, suffix = "hc1"),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_multi_endo_firststage_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_multi_endo", vce_combo$suffix)

  if (file.exists(sim_multi_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                  data = sim_multi, vcov = vce_combo$vcov,
                  small = vce_combo$small)
    for (endo in c("endo1", "endo2")) {
      check_firststage(fit, fixture_file, endo,
                       paste(label, endo))
    }
  }
}

# Multi-endo: Shea != partial R² (structural)
test_that("multi-endo: Shea partial R2 differs from partial R2", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi)
  # Shea should be strictly less than partial R² for multi-endo
  expect_true(fit$first_stage$endo1$shea_partial_r2 <
              fit$first_stage$endo1$partial_r2)
  expect_true(fit$first_stage$endo2$shea_partial_r2 <
              fit$first_stage$endo2$partial_r2)
})


# ============================================================================
# sim_no_constant: y ~ 0+x1 | endo1 | z1+z2
# No intercept
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC1",  small = FALSE, suffix = "hc1"),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_no_constant_firststage_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_no_constant", vce_combo$suffix)

  if (file.exists(sim_noconst_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2,
                  data = sim_noconst, vcov = vce_combo$vcov,
                  small = vce_combo$small)
    check_firststage(fit, fixture_file, "endo1", label)
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
    paste0("sim_cluster_firststage_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_cluster", vce_combo$suffix)

  if (file.exists(sim_cluster_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                  data = sim_cluster, vcov = vce_combo$vcov,
                  small = vce_combo$small, clusters = vce_combo$clusters)
    check_firststage(fit, fixture_file, "endo1", label)
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
    paste0("card_just_id_weighted_firststage_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_just_id_weighted", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card, weights = weight,
                  vcov = vce_combo$vcov, small = vce_combo$small,
                  clusters = vce_combo$clusters)
    check_firststage(fit, fixture_file, "educ", label)
  }
}
