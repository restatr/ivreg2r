# ============================================================================
# Tests: Instrument Redundancy Test (Ticket P1)
# ============================================================================
#
# Tests H0: specified excluded instruments are redundant (zero explanatory
# power for endogenous regressors) conditional on maintained instruments.
# KP rk LM test of H0: rank=0. Verifies against Stata ivreg2 fixtures
# (e(redstat), e(redp), e(reddf)).

# --- Helpers ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)

read_redundancy_fixture <- function(path) {
  read.csv(path, check.names = FALSE)
}

# --- Load datasets ---
card_path <- file.path(fixture_dir, "card_data.csv")
if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

sim_multi_endo_path <- file.path(fixture_dir, "sim_multi_endo_data.csv")
if (file.exists(sim_multi_endo_path)) {
  sim_multi_endo <- read.csv(sim_multi_endo_path)
}

ts_hac_path <- file.path(fixture_dir, "ts_hac_data.csv")
if (file.exists(ts_hac_path)) {
  ts_hac <- read.csv(ts_hac_path)
}


# ============================================================================
# Structural tests (no fixtures)
# ============================================================================

test_that("redundancy is NULL when not requested", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4, data = card)
  expect_null(fit$diagnostics$redundancy)
})

test_that("redundant = character(0) is silently skipped", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, redundant = character(0))
  expect_null(fit$diagnostics$redundancy)
})

test_that("redundancy is NULL for OLS models", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, redundant = "wt")
  expect_null(fit$diagnostics)
})

test_that("redundancy has all expected fields", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, redundant = "nearc2")
  red <- fit$diagnostics$redundancy
  expect_type(red, "list")
  expect_named(red, c("stat", "p", "df", "test_name", "tested_vars"),
               ignore.order = TRUE)
  expect_equal(red$tested_vars, "nearc2")
  expect_identical(red$df, 1L)
  expect_equal(red$test_name, "Redundancy test (LM)")
})

test_that("invalid redundant variable name produces error", {
  skip_if(!file.exists(card_path), "card data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, redundant = "not_a_var"),
    "not in the excluded instrument list"
  )
})

test_that("redundant rejects exogenous variable", {
  skip_if(!file.exists(card_path), "card data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, redundant = "exper"),
    "not in the excluded instrument list"
  )
})

test_that("redundant rejects endogenous variable", {
  skip_if(!file.exists(card_path), "card data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, redundant = "educ"),
    "not in the excluded instrument list"
  )
})

test_that("redundant must be character or NULL", {
  skip_if(!file.exists(card_path), "card data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, redundant = 1),
    "character vector or NULL"
  )
})

test_that("duplicate redundant entries are deduplicated with warning", {
  skip_if(!file.exists(card_path), "card data not found")
  fit_single <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                       data = card, redundant = "nearc2")
  expect_warning(
    fit_dup <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                      data = card, redundant = c("nearc2", "nearc2")),
    "duplicate entries"
  )
  expect_equal(fit_dup$diagnostics$redundancy$stat,
               fit_single$diagnostics$redundancy$stat)
  expect_equal(fit_dup$diagnostics$redundancy$df,
               fit_single$diagnostics$redundancy$df)
})

test_that("small=TRUE and small=FALSE produce identical redundancy stats (IID)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, redundant = "nearc2", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, redundant = "nearc2", small = TRUE)
  expect_equal(fit1$diagnostics$redundancy$stat,
               fit2$diagnostics$redundancy$stat)
})

test_that("small=TRUE and small=FALSE produce identical redundancy stats (HC1)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, vcov = "HC1", redundant = "nearc2", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, vcov = "HC1", redundant = "nearc2", small = TRUE)
  expect_equal(fit1$diagnostics$redundancy$stat,
               fit2$diagnostics$redundancy$stat)
})

test_that("HC0 and HC1 produce identical redundancy stat", {
  skip_if(!file.exists(card_path), "card data not found")
  fit0 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, vcov = "HC0", redundant = "nearc2")
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, vcov = "HC1", redundant = "nearc2")
  expect_equal(fit0$diagnostics$redundancy$stat,
               fit1$diagnostics$redundancy$stat)
})

test_that("redundancy appears in glance output", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, redundant = "nearc2")
  g <- glance(fit)
  expect_true("redundancy_stat" %in% names(g))
  expect_true("redundancy_p" %in% names(g))
  expect_false(is.na(g$redundancy_stat))
  expect_false(is.na(g$redundancy_p))
})

test_that("redundancy is NA in glance when not requested", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card)
  g <- glance(fit)
  expect_true(is.na(g$redundancy_stat))
  expect_true(is.na(g$redundancy_p))
})


# ============================================================================
# Property tests
# ============================================================================

test_that("redundancy stat is non-negative", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, redundant = "nearc2")
  expect_gte(fit$diagnostics$redundancy$stat, 0)
})

test_that("df = K1 * L_tested", {
  skip_if(!file.exists(card_path), "card data not found")
  # K1=1, L_tested=1 → df=1
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, redundant = "nearc2")
  expect_identical(fit1$diagnostics$redundancy$df, 1L)

  # K1=1, L_tested=2 → df=2
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, redundant = c("nearc2", "nearc4"))
  expect_identical(fit2$diagnostics$redundancy$df, 2L)
})

test_that("df = K1 * L_tested for multi-endo models", {
  skip_if(!file.exists(sim_multi_endo_path), "sim_multi_endo data not found")
  # K1=2, L_tested=1 → df=2
  fit1 <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                 data = sim_multi_endo, redundant = "z1")
  expect_identical(fit1$diagnostics$redundancy$df, 2L)

  # K1=2, L_tested=2 → df=4
  fit2 <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                 data = sim_multi_endo, redundant = c("z1", "z2"))
  expect_identical(fit2$diagnostics$redundancy$df, 4L)
})

test_that("redundancy suppressed by b0 (no diagnostics)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, redundant = "nearc2",
                b0 = c("(Intercept)" = 0, exper = 0, expersq = 0, educ = 0))
  # b0 suppresses all id diagnostics
  expect_null(fit$diagnostics$redundancy)
})


# ============================================================================
# Helper: run fixture comparison for one redundancy configuration
# ============================================================================

check_redundancy <- function(fit, fixture_path, label) {
  fixture <- read_redundancy_fixture(fixture_path)
  red <- fit$diagnostics$redundancy

  test_that(paste("Redundancy stat matches Stata", label), {
    expect_false(is.null(red))
    expect_equal(red$stat, fixture$redstat,
                 tolerance = stata_tol$stat)
  })

  test_that(paste("Redundancy p-value matches Stata", label), {
    expect_equal(red$p, fixture$redp,
                 tolerance = stata_tol$pval)
  })

  test_that(paste("Redundancy df matches Stata", label), {
    expect_identical(red$df, as.integer(fixture$reddf))
  })
}


# ============================================================================
# Card overid: redundant(nearc2), IID / HC0 / HC1 small / cluster
# Model: lwage ~ exper + expersq | educ | nearc2 + nearc4
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, clusters = NULL, suffix = "nearc2_iid"),
  list(vcov = "HC0",  small = FALSE, clusters = NULL, suffix = "nearc2_hc0"),
  list(vcov = "HC1",  small = TRUE,  clusters = NULL, suffix = "nearc2_hc1_small"),
  list(vcov = "HC1",  small = FALSE, clusters = ~smsa, suffix = "nearc2_cluster")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("card_overid_redundancy_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_overid", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                  data = card, vcov = vce_combo$vcov, small = vce_combo$small,
                  clusters = vce_combo$clusters, redundant = "nearc2")
    check_redundancy(fit, fixture_file, label)
  }
}


# ============================================================================
# Card overid: redundant(nearc4), IID — different instrument
# ============================================================================

test_that("redundant(nearc4) IID matches Stata", {
  fixture_file <- file.path(fixture_dir, "card_overid_redundancy_nearc4_iid.csv")
  skip_if(!file.exists(card_path), "card data not found")
  skip_if(!file.exists(fixture_file), "fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, redundant = "nearc4")
  fixture <- read_redundancy_fixture(fixture_file)
  expect_equal(fit$diagnostics$redundancy$stat, fixture$redstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$redundancy$p, fixture$redp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$redundancy$df, as.integer(fixture$reddf))
})


# ============================================================================
# Card overid: redundant(nearc2 nearc4), IID — test all excluded
# ============================================================================

test_that("redundant(nearc2 nearc4) IID matches Stata", {
  fixture_file <- file.path(fixture_dir, "card_overid_redundancy_both_iid.csv")
  skip_if(!file.exists(card_path), "card data not found")
  skip_if(!file.exists(fixture_file), "fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, redundant = c("nearc2", "nearc4"))
  fixture <- read_redundancy_fixture(fixture_file)
  expect_equal(fit$diagnostics$redundancy$stat, fixture$redstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$redundancy$p, fixture$redp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$redundancy$df, as.integer(fixture$reddf))
})


# ============================================================================
# Card overid + dofminus: redundant(nearc2), IID and cluster
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  clusters = NULL,   suffix = "nearc2_iid"),
  list(vcov = "HC1",  clusters = ~smsa,  suffix = "nearc2_cluster")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("card_overid_dof_redundancy_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_overid_dof", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                  data = card, vcov = vce_combo$vcov,
                  clusters = vce_combo$clusters,
                  redundant = "nearc2", dofminus = 1L)
    check_redundancy(fit, fixture_file, label)
  }
}


# ============================================================================
# Card overid + aweight: redundant(nearc2), IID
# ============================================================================

test_that("redundant(nearc2) aweight IID matches Stata", {
  fixture_file <- file.path(fixture_dir, "card_overid_aw_redundancy_nearc2_iid.csv")
  skip_if(!file.exists(card_path), "card data not found")
  skip_if(!file.exists(fixture_file), "fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, redundant = "nearc2", weights = weight)
  fixture <- read_redundancy_fixture(fixture_file)
  expect_equal(fit$diagnostics$redundancy$stat, fixture$redstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$redundancy$p, fixture$redp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$redundancy$df, as.integer(fixture$reddf))
})


# ============================================================================
# Card overid + gmm2s: redundant(nearc2), robust
# ============================================================================

test_that("redundant(nearc2) gmm2s robust matches Stata", {
  fixture_file <- file.path(fixture_dir, "card_overid_gmm_redundancy_nearc2_robust.csv")
  skip_if(!file.exists(card_path), "card data not found")
  skip_if(!file.exists(fixture_file), "fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, redundant = "nearc2",
                method = "gmm2s", vcov = "HC1")
  fixture <- read_redundancy_fixture(fixture_file)
  expect_equal(fit$diagnostics$redundancy$stat, fixture$redstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$redundancy$p, fixture$redp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$redundancy$df, as.integer(fixture$reddf))
})


# ============================================================================
# sim_multi_endo: redundant(z1), IID and HC1 — K1=2 path
# Model: y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid", suffix = "z1_iid"),
  list(vcov = "HC1", suffix = "z1_hc1")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_multi_endo_redundancy_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_multi_endo", vce_combo$suffix)

  if (file.exists(sim_multi_endo_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                  data = sim_multi_endo, vcov = vce_combo$vcov,
                  redundant = "z1")
    check_redundancy(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_multi_endo: redundant(z1 z2), IID — multi-instrument multi-endo
# ============================================================================

test_that("redundant(z1 z2) multi-endo IID matches Stata", {
  fixture_file <- file.path(fixture_dir, "sim_multi_endo_redundancy_z1z2_iid.csv")
  skip_if(!file.exists(sim_multi_endo_path), "sim_multi_endo data not found")
  skip_if(!file.exists(fixture_file), "fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, redundant = c("z1", "z2"))
  fixture <- read_redundancy_fixture(fixture_file)
  expect_equal(fit$diagnostics$redundancy$stat, fixture$redstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$redundancy$p, fixture$redp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$redundancy$df, as.integer(fixture$reddf))
})


# ============================================================================
# ts_hac: redundant(z1), Bartlett bw=3 — HAC path
# Model: y ~ w | x | z1 + z2
# ============================================================================

test_that("redundant(z1) HAC Bartlett bw=3 matches Stata", {
  fixture_file <- file.path(fixture_dir, "ts_hac_redundancy_z1_bartlett_bw3.csv")
  skip_if(!file.exists(ts_hac_path), "ts_hac data not found")
  skip_if(!file.exists(fixture_file), "fixture not found")

  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_hac,
                redundant = "z1", vcov = "HAC",
                kernel = "Bartlett", bw = 3, tvar = "t")
  fixture <- read_redundancy_fixture(fixture_file)
  expect_equal(fit$diagnostics$redundancy$stat, fixture$redstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$redundancy$p, fixture$redp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$redundancy$df, as.integer(fixture$reddf))
})
