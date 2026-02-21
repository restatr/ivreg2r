# ============================================================================
# Tests: Overidentification Tests (Ticket D1)
# ============================================================================
#
# Sargan (IID) and Hansen J (robust/cluster) overidentification tests.
# Verifies against Stata ivreg2 fixtures.

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

sim_path <- file.path(fixture_dir, "sim_cluster_data.csv")
if (file.exists(sim_path)) {
  sim_cluster <- read.csv(sim_path)
}


# ============================================================================
# Sargan test: IID, card_overid
# ============================================================================

test_that("Sargan stat matches Stata card_overid iid fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_overid_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Sargan")
  expect_equal(fit$diagnostics$overid$stat, fixture$sargan,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$sarganp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$sargandf))
})


# ============================================================================
# Sargan test: IID small — same stat as non-small
# ============================================================================

test_that("Sargan stat matches Stata card_overid iid_small fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_overid_diagnostics_iid_small.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, small = TRUE)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, fixture$sargan,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$sarganp,
               tolerance = stata_tol$pval)
})

test_that("small does not change Sargan statistic", {
  skip_if(!file.exists(card_path), "card data not found")

  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                 data = card, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                 data = card, small = TRUE)
  expect_equal(fit1$diagnostics$overid$stat, fit2$diagnostics$overid$stat)
})


# ============================================================================
# Hansen J test: HC, card_overid
# ============================================================================

test_that("Hansen J stat matches Stata card_overid hc1 fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_overid_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$jdf))
})

test_that("Hansen J stat matches Stata card_overid hc1_small fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_overid_diagnostics_hc1_small.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "HC1", small = TRUE)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
})

test_that("small does not change Hansen J statistic (HC)", {
  skip_if(!file.exists(card_path), "card data not found")

  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                 data = card, vcov = "HC1", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                 data = card, vcov = "HC1", small = TRUE)
  expect_equal(fit1$diagnostics$overid$stat, fit2$diagnostics$overid$stat)
})


# ============================================================================
# Hansen J test: Cluster, sim_cluster
# ============================================================================

test_that("Hansen J stat matches Stata sim_cluster cl fixture", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")
  diag_path <- file.path(fixture_dir, "sim_cluster_diagnostics_cl.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$jdf))
})

test_that("Hansen J stat matches Stata sim_cluster cl_small fixture", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")
  diag_path <- file.path(fixture_dir, "sim_cluster_diagnostics_cl_small.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id, small = TRUE)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
})

test_that("small does not change Hansen J statistic (cluster)", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")

  fit1 <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                 data = sim_cluster, clusters = ~cluster_id, small = FALSE)
  fit2 <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                 data = sim_cluster, clusters = ~cluster_id, small = TRUE)
  expect_equal(fit1$diagnostics$overid$stat, fit2$diagnostics$overid$stat)
})


# ============================================================================
# Exactly identified: stat = 0, p = NA
# ============================================================================

test_that("exactly identified IID model has stat=0, df=0, p=NA", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_just_id_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, 0)
  expect_true(is.na(fit$diagnostics$overid$p))
  expect_identical(fit$diagnostics$overid$df, 0L)
  expect_equal(fit$diagnostics$overid$test_name, "Sargan")
})

test_that("exactly identified HC model has stat=0, df=0, p=NA", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_just_id_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, 0)
  expect_true(is.na(fit$diagnostics$overid$p))
  expect_identical(fit$diagnostics$overid$df, 0L)
  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
})


# ============================================================================
# OLS: diagnostics is NULL
# ============================================================================

test_that("OLS model has NULL diagnostics", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$diagnostics)
})


# ============================================================================
# HC0 produces same J as HC1 (no small-sample correction on Omega)
# ============================================================================

test_that("HC0 and HC1 produce identical Hansen J statistic", {
  skip_if(!file.exists(card_path), "card data not found")

  fit_hc0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "HC0")
  fit_hc1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "HC1")
  expect_equal(fit_hc0$diagnostics$overid$stat,
               fit_hc1$diagnostics$overid$stat)
})


# ============================================================================
# sim_no_constant: y ~ 0 + x1 | endo1 | z1+z2 (noconstant, overid_df=1)
# ============================================================================

sim_noconst_path <- file.path(fixture_dir, "sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

test_that("Sargan stat matches Stata sim_no_constant iid fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- file.path(fixture_dir, "sim_no_constant_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Sargan")
  expect_equal(fit$diagnostics$overid$stat, fixture$sargan,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$sarganp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$sargandf))
})

test_that("Sargan stat matches Stata sim_no_constant iid_small fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- file.path(fixture_dir, "sim_no_constant_diagnostics_iid_small.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, small = TRUE)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, fixture$sargan,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$sarganp,
               tolerance = stata_tol$pval)
})

test_that("Hansen J stat matches Stata sim_no_constant hc1 fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- file.path(fixture_dir, "sim_no_constant_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$jdf))
})

test_that("Hansen J stat matches Stata sim_no_constant hc1_small fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- file.path(fixture_dir, "sim_no_constant_diagnostics_hc1_small.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst,
                vcov = "HC1", small = TRUE)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
})

test_that("small does not change overid statistic (sim_no_constant)", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")

  fit1 <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, small = FALSE)
  fit2 <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, small = TRUE)
  expect_equal(fit1$diagnostics$overid$stat, fit2$diagnostics$overid$stat)

  fit3 <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst,
                 vcov = "HC1", small = FALSE)
  fit4 <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst,
                 vcov = "HC1", small = TRUE)
  expect_equal(fit3$diagnostics$overid$stat, fit4$diagnostics$overid$stat)
})


# ============================================================================
# sim_multi_endo: y ~ x1+x2 | endo1+endo2 | z1+z2+z3+z4 (overid_df=2)
# ============================================================================

sim_multi_path <- file.path(fixture_dir, "sim_multi_endo_data.csv")
if (file.exists(sim_multi_path)) {
  sim_multi <- read.csv(sim_multi_path)
}

test_that("Sargan stat matches Stata sim_multi_endo iid fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- file.path(fixture_dir, "sim_multi_endo_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Sargan")
  expect_equal(fit$diagnostics$overid$stat, fixture$sargan,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$sarganp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$sargandf))
})

test_that("Sargan stat matches Stata sim_multi_endo iid_small fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- file.path(fixture_dir, "sim_multi_endo_diagnostics_iid_small.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi, small = TRUE)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, fixture$sargan,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$sarganp,
               tolerance = stata_tol$pval)
})

test_that("Hansen J stat matches Stata sim_multi_endo hc1 fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- file.path(fixture_dir, "sim_multi_endo_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$jdf))
})

test_that("Hansen J stat matches Stata sim_multi_endo hc1_small fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- file.path(fixture_dir, "sim_multi_endo_diagnostics_hc1_small.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi, vcov = "HC1", small = TRUE)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
})

test_that("small does not change overid statistic (sim_multi_endo)", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")

  fit1 <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                 data = sim_multi, small = FALSE)
  fit2 <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                 data = sim_multi, small = TRUE)
  expect_equal(fit1$diagnostics$overid$stat, fit2$diagnostics$overid$stat)

  fit3 <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                 data = sim_multi, vcov = "HC1", small = FALSE)
  fit4 <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                 data = sim_multi, vcov = "HC1", small = TRUE)
  expect_equal(fit3$diagnostics$overid$stat, fit4$diagnostics$overid$stat)
})
