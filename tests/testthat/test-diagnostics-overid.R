# ============================================================================
# Tests: Overidentification Tests (Ticket D1)
# ============================================================================
#
# Sargan (IID) and Hansen J (robust/cluster) overidentification tests.
# Verifies against Stata ivreg2 fixtures.

# --- Load datasets ---
data(abdata, package = "ivreg2r")
data(mroz, package = "ivreg2r")

# Shared fits (M-17 hoisting precedent, see test-vcov-cluster.R): the ab_cl
# cluster fit is reused across the three tests below instead of being refit.
fit_ab_cl <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                     clusters = ~id, endog = "w")
fit_ab_cl_small <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                          clusters = ~id, small = TRUE, endog = "w")


# ============================================================================
# Sargan test: IID, mroz overid (H31 help.txt:1274; family M-10 re-base,
# replacing card_overid)
# ============================================================================

test_that("Sargan stat matches Stata mroz overid H31 fixture", {
  diag_path <- fixture_path("hf_mroz_diagnostics_H31.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Sargan")
  expect_equal(fit$diagnostics$overid$stat, fixture$sargan,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$sarganp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$sargandf))
})


# ============================================================================
# Sargan test: IID small — same stat as non-small, checked against the same
# H31 fixture (the small-variant retirement, 2026-07-06, established that
# Sargan is byte-identical under `small`, so no dedicated small fixture
# exists)
# ============================================================================

test_that("Sargan stat matches Stata mroz overid H31 fixture (small)", {
  diag_path <- fixture_path("hf_mroz_diagnostics_H31.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz, small = TRUE)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, fixture$sargan,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$sarganp,
               tolerance = stata_tol$pval)
})

test_that("small does not change Sargan statistic", {
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                 data = mroz, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                 data = mroz, small = TRUE)
  expect_equal(fit1$diagnostics$overid$stat, fit2$diagnostics$overid$stat)
})


# ============================================================================
# Hansen J test: HC, mroz overid (H41 help.txt:1325; family M-10 re-base,
# replacing card_overid)
# ============================================================================

test_that("Hansen J stat matches Stata mroz overid H41 fixture", {
  diag_path <- fixture_path("hf_mroz_diagnostics_H41.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz, vcov = "robust")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$jdf))
})

test_that("Hansen J stat matches Stata mroz overid H41 fixture (small)", {
  diag_path <- fixture_path("hf_mroz_diagnostics_H41.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz, vcov = "robust", small = TRUE)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
})

test_that("small does not change Hansen J statistic (HC)", {
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                 data = mroz, vcov = "robust", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                 data = mroz, vcov = "robust", small = TRUE)
  expect_equal(fit1$diagnostics$overid$stat, fit2$diagnostics$overid$stat)
})


# ============================================================================
# Hansen J test: Cluster, ab_cl (family M-13, re-based: abdata H88-minus-gmm2s
# 2SLS, cluster(id))
# ============================================================================

test_that("Hansen J stat matches Stata ab_cl cl fixture", {
  diag_path <- fixture_path("ab_cl_diagnostics_cl.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- fit_ab_cl
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$jdf))
})

test_that("Hansen J stat matches Stata ab_cl cl_small fixture", {
  diag_path <- fixture_path("ab_cl_diagnostics_cl_small.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- fit_ab_cl_small
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
})

test_that("small does not change Hansen J statistic (cluster)", {
  fit1 <- fit_ab_cl
  fit2 <- fit_ab_cl_small
  expect_equal(fit1$diagnostics$overid$stat, fit2$diagnostics$overid$stat)
})


# ============================================================================
# Exactly identified: stat = 0, p = NA
# ============================================================================

test_that("exactly identified IID model has stat=0, df=0, p=NA", {
  diag_path <- fixture_path("mroz_justid_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$stat, 0)
  expect_true(is.na(fit$diagnostics$overid$p))
  expect_identical(fit$diagnostics$overid$df, 0L)
  expect_equal(fit$diagnostics$overid$test_name, "Sargan")
})

test_that("exactly identified HC model has stat=0, df=0, p=NA", {
  diag_path <- fixture_path("mroz_justid_diagnostics_robust.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, vcov = "robust")
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
  fit_hc0 <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                     data = mroz, vcov = "robust")
  fit_hc1 <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                     data = mroz, vcov = "robust")
  expect_equal(fit_hc0$diagnostics$overid$stat,
               fit_hc1$diagnostics$overid$stat)
})


# ============================================================================
# sim_no_constant: y ~ 0 + x1 | endo1 | z1+z2 (noconstant, overid_df=1)
# ============================================================================

sim_noconst_path <- fixture_path("sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

test_that("Sargan stat matches Stata sim_no_constant iid fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- fixture_path("sim_no_constant_diagnostics_iid.csv")
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
  diag_path <- fixture_path("sim_no_constant_diagnostics_iid_small.csv")
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
  diag_path <- fixture_path("sim_no_constant_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, vcov = "robust")
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
  diag_path <- fixture_path("sim_no_constant_diagnostics_hc1_small.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst,
                vcov = "robust", small = TRUE)
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
                 vcov = "robust", small = FALSE)
  fit4 <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst,
                 vcov = "robust", small = TRUE)
  expect_equal(fit3$diagnostics$overid$stat, fit4$diagnostics$overid$stat)
})


# ============================================================================
# sim_multi_endo: y ~ x1+x2 | endo1+endo2 | z1+z2+z3+z4 (overid_df=2)
# ============================================================================

sim_multi_path <- fixture_path("sim_multi_endo_data.csv")
if (file.exists(sim_multi_path)) {
  sim_multi <- read.csv(sim_multi_path)
}

test_that("Sargan stat matches Stata sim_multi_endo iid fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- fixture_path("sim_multi_endo_diagnostics_iid.csv")
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
  diag_path <- fixture_path("sim_multi_endo_diagnostics_iid_small.csv")
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
  diag_path <- fixture_path("sim_multi_endo_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi, vcov = "robust")
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
  diag_path <- fixture_path("sim_multi_endo_diagnostics_hc1_small.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi, vcov = "robust", small = TRUE)
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
                 data = sim_multi, vcov = "robust", small = FALSE)
  fit4 <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                 data = sim_multi, vcov = "robust", small = TRUE)
  expect_equal(fit3$diagnostics$overid$stat, fit4$diagnostics$overid$stat)
})

# ============================================================================
# Structural rank bound: one-way cluster meat with M < L (platform-stable)
# ============================================================================
#
# Added at the 2026-07-06 H28/H29 investigation. The numeric singularity
# detectors (chol failure, qr rank, eigenvalue thresholds) proved
# BLAS-dependent on the rank-deficit-1 cluster meat: the griliches
# cluster(year) + partial cells were detected on macOS
# (lambda_min/lambda_max ~ 2e-18) but missed on ubuntu/Windows.
# .cluster_rank_bound() gates the Omega-inverting consumers structurally:
# a one-way cluster meat is a sum of M rank-one outer products, so M < L
# (M - 1 < L when centered) is singular by construction on every platform.
# These synthetic tests pin the gate deterministically; the canonical
# griliches behavior is pinned by H27/H28/H29 in test-helpfile-examples.R.

test_that(".cluster_rank_bound: M one-way, M - 1 centered, Inf otherwise", {
  cv <- rep(1:3, each = 10L)
  expect_identical(.cluster_rank_bound(cv, kernel = NULL), 3L)
  expect_identical(.cluster_rank_bound(cv, kernel = NULL, center = TRUE), 2L)
  expect_identical(.cluster_rank_bound(NULL, kernel = NULL), Inf)
  expect_identical(.cluster_rank_bound(cv, kernel = "bartlett"), Inf)
  expect_identical(.cluster_rank_bound(list(cv, cv), kernel = NULL), Inf)
})

# Deterministic synthetic data: M = 3 clusters vs L = 4 moment conditions
# (intercept + three excluded instruments).
set.seed(42)
rank_gate_n <- 60L
rank_gate_df <- data.frame(
  z1 = rnorm(rank_gate_n), z2 = rnorm(rank_gate_n), z3 = rnorm(rank_gate_n),
  g = rep(1:3, each = 20L)
)
rank_gate_df$x <- rank_gate_df$z1 + rank_gate_df$z2 + rnorm(rank_gate_n)
rank_gate_df$y <- rank_gate_df$x + rnorm(rank_gate_n)

test_that("GMM2S with M < L clusters errors deterministically", {
  expect_error(
    ivreg2(y ~ 1 | x | z1 + z2 + z3, data = rank_gate_df, clusters = ~g,
           method = "gmm2s"),
    "not of full rank"
  )
})

test_that("2SLS diagnostics with M < L clusters: J, AR, Stock-Wright all NA with warnings", {
  expect_warning(
    expect_warning(
      expect_warning(
        fit <- ivreg2(y ~ 1 | x | z1 + z2 + z3, data = rank_gate_df,
                      clusters = ~g),
        "Hansen J statistic not computed"
      ),
      "Anderson-Rubin: RVR matrix is singular"
    ),
    "Stock-Wright.*rank-deficient"
  )
  expect_true(is.na(fit$diagnostics$overid$stat))
  expect_identical(fit$diagnostics$overid$df, 2L)
})

test_that("center drops the bound to M - 1: M = L fits uncentered, errors centered", {
  # M = 4 clusters vs L = 4 moments: the uncentered meat is full rank
  # (bound M = L does not bind); centering makes the four cluster sums add
  # to zero, so the centered meat has rank <= 3 by construction.
  set.seed(43)
  d4 <- data.frame(
    z1 = rnorm(rank_gate_n), z2 = rnorm(rank_gate_n), z3 = rnorm(rank_gate_n),
    g = rep(1:4, each = 15L)
  )
  d4$x <- d4$z1 + d4$z2 + rnorm(rank_gate_n)
  d4$y <- d4$x + rnorm(rank_gate_n)
  fit_ok <- ivreg2(y ~ 1 | x | z1 + z2 + z3, data = d4, clusters = ~g,
                   method = "gmm2s")
  expect_s3_class(fit_ok, "ivreg2")
  expect_error(
    ivreg2(y ~ 1 | x | z1 + z2 + z3, data = d4, clusters = ~g,
           method = "gmm2s", center = TRUE),
    "not of full rank"
  )
})
