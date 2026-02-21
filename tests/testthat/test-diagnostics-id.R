# ============================================================================
# Tests: Identification Diagnostics (Ticket D2)
# ============================================================================
#
# Anderson LM / Cragg-Donald F (IID) and Kleibergen-Paap rk LM / rk Wald F
# (robust/cluster) underidentification and weak identification tests.
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

sim_multi_path <- file.path(fixture_dir, "sim_multi_endo_data.csv")
if (file.exists(sim_multi_path)) {
  sim_multi <- read.csv(sim_multi_path)
}

sim_cluster_path <- file.path(fixture_dir, "sim_cluster_data.csv")
if (file.exists(sim_cluster_path)) {
  sim_cluster <- read.csv(sim_cluster_path)
}


# ============================================================================
# Anderson LM — IID
# ============================================================================

test_that("Anderson LM matches Stata card_just_id iid fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_just_id_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$test_name,
               "Anderson canon. corr. LM statistic")
  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("Anderson LM matches Stata card_overid iid fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_overid_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("Anderson LM matches Stata sim_multi_endo iid fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- file.path(fixture_dir, "sim_multi_endo_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})


# ============================================================================
# Cragg-Donald F — IID
# ============================================================================

test_that("Cragg-Donald F matches Stata card_just_id iid fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_just_id_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$test_name,
               "Cragg-Donald Wald F statistic")
  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})

test_that("Cragg-Donald F matches Stata card_overid iid fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_overid_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})

test_that("Cragg-Donald F matches Stata sim_multi_endo iid fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- file.path(fixture_dir, "sim_multi_endo_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})


# ============================================================================
# KP rk LM — HC1
# ============================================================================

test_that("KP rk LM matches Stata card_just_id hc1 fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_just_id_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$test_name,
               "Kleibergen-Paap rk LM statistic")
  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("KP rk LM matches Stata card_overid hc1 fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_overid_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("KP rk LM matches Stata sim_multi_endo hc1 fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- file.path(fixture_dir, "sim_multi_endo_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})


# ============================================================================
# KP rk Wald F — HC1
# ============================================================================

test_that("KP rk Wald F matches Stata card_just_id hc1 fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_just_id_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id_robust$test_name,
               "Kleibergen-Paap rk Wald F statistic")
  expect_equal(fit$diagnostics$weak_id_robust$stat, fixture$widstat,
               tolerance = stata_tol$stat)
})

test_that("KP rk Wald F matches Stata card_overid hc1 fixture", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_overid_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id_robust$stat, fixture$widstat,
               tolerance = stata_tol$stat)
})

test_that("KP rk Wald F matches Stata sim_multi_endo hc1 fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- file.path(fixture_dir, "sim_multi_endo_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id_robust$stat, fixture$widstat,
               tolerance = stata_tol$stat)
})


# ============================================================================
# KP rk LM — Cluster
# ============================================================================

test_that("KP rk LM matches Stata sim_cluster cl fixture", {
  skip_if(!file.exists(sim_cluster_path), "sim_cluster data not found")
  diag_path <- file.path(fixture_dir, "sim_cluster_diagnostics_cl.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$test_name,
               "Kleibergen-Paap rk LM statistic")
  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})


# ============================================================================
# KP rk Wald F — Cluster
# ============================================================================

test_that("KP rk Wald F matches Stata sim_cluster cl fixture", {
  skip_if(!file.exists(sim_cluster_path), "sim_cluster data not found")
  diag_path <- file.path(fixture_dir, "sim_cluster_diagnostics_cl.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id_robust$stat, fixture$widstat,
               tolerance = stata_tol$stat)
})


# ============================================================================
# CD F always present even with robust/cluster
# ============================================================================

test_that("Cragg-Donald F present alongside KP stats (HC1)", {
  skip_if(!file.exists(card_path), "card data not found")
  diag_path <- file.path(fixture_dir, "card_just_id_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$weak_id$test_name,
               "Cragg-Donald Wald F statistic")
})

test_that("Cragg-Donald F present alongside KP stats (cluster)", {
  skip_if(!file.exists(sim_cluster_path), "sim_cluster data not found")
  diag_path <- file.path(fixture_dir, "sim_cluster_diagnostics_cl.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})


# ============================================================================
# Test name dispatch: IID vs robust
# ============================================================================

test_that("IID reports Anderson LM, robust reports KP rk LM", {
  skip_if(!file.exists(card_path), "card data not found")

  fit_iid <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card)
  fit_hc1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card, vcov = "HC1")

  expect_equal(fit_iid$diagnostics$underid$test_name,
               "Anderson canon. corr. LM statistic")
  expect_equal(fit_hc1$diagnostics$underid$test_name,
               "Kleibergen-Paap rk LM statistic")
})


# ============================================================================
# weak_id_robust is NULL for IID
# ============================================================================

test_that("weak_id_robust is NULL for IID models", {
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  expect_null(fit$diagnostics$weak_id_robust)
})


# ============================================================================
# small flag invariance: stats identical for small=TRUE vs FALSE
# ============================================================================

test_that("small does not change Anderson LM statistic", {
  skip_if(!file.exists(card_path), "card data not found")

  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = TRUE)
  expect_equal(fit1$diagnostics$underid$stat, fit2$diagnostics$underid$stat)
  expect_equal(fit1$diagnostics$weak_id$stat, fit2$diagnostics$weak_id$stat)
})

test_that("small does not change KP rk LM statistic (HC1)", {
  skip_if(!file.exists(card_path), "card data not found")

  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1", small = TRUE)
  expect_equal(fit1$diagnostics$underid$stat, fit2$diagnostics$underid$stat)
  expect_equal(fit1$diagnostics$weak_id_robust$stat,
               fit2$diagnostics$weak_id_robust$stat)
})

test_that("small does not change KP stats (cluster)", {
  skip_if(!file.exists(sim_cluster_path), "sim_cluster data not found")

  fit1 <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                 data = sim_cluster, clusters = ~cluster_id, small = FALSE)
  fit2 <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                 data = sim_cluster, clusters = ~cluster_id, small = TRUE)
  expect_equal(fit1$diagnostics$underid$stat, fit2$diagnostics$underid$stat)
  expect_equal(fit1$diagnostics$weak_id_robust$stat,
               fit2$diagnostics$weak_id_robust$stat)
})


# ============================================================================
# HC0 = HC1 equivalence: KP stats identical
# ============================================================================

test_that("HC0 and HC1 produce identical KP rk LM and Wald F", {
  skip_if(!file.exists(card_path), "card data not found")

  fit_hc0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, vcov = "HC0")
  fit_hc1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, vcov = "HC1")

  expect_equal(fit_hc0$diagnostics$underid$stat,
               fit_hc1$diagnostics$underid$stat)
  expect_equal(fit_hc0$diagnostics$weak_id_robust$stat,
               fit_hc1$diagnostics$weak_id_robust$stat)
})


# ============================================================================
# OLS: identification diagnostics are NULL
# ============================================================================

test_that("OLS model has NULL identification diagnostics", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$diagnostics$underid)
  expect_null(fit$diagnostics$weak_id)
  expect_null(fit$diagnostics$weak_id_robust)
})


# ============================================================================
# CD F is VCE-invariant
# ============================================================================

test_that("Cragg-Donald F is identical across VCE types", {
  skip_if(!file.exists(card_path), "card data not found")

  fit_iid <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                    data = card)
  fit_hc1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                    data = card, vcov = "HC1")

  expect_equal(fit_iid$diagnostics$weak_id$stat,
               fit_hc1$diagnostics$weak_id$stat)
})


# ============================================================================
# sim_no_constant: y ~ 0 + x1 | endo1 | z1+z2 (noconstant)
# ============================================================================

sim_noconst_path <- file.path(fixture_dir, "sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

test_that("Anderson LM matches Stata sim_no_constant iid fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- file.path(fixture_dir, "sim_no_constant_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("Cragg-Donald F matches Stata sim_no_constant iid fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- file.path(fixture_dir, "sim_no_constant_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})

test_that("KP rk LM matches Stata sim_no_constant hc1 fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- file.path(fixture_dir, "sim_no_constant_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("KP rk Wald F matches Stata sim_no_constant hc1 fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- file.path(fixture_dir, "sim_no_constant_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id_robust$stat, fixture$widstat,
               tolerance = stata_tol$stat)
})

test_that("small does not change id stats (sim_no_constant)", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")

  fit1 <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, small = FALSE)
  fit2 <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, small = TRUE)
  expect_equal(fit1$diagnostics$underid$stat, fit2$diagnostics$underid$stat)
  expect_equal(fit1$diagnostics$weak_id$stat, fit2$diagnostics$weak_id$stat)
})

test_that("Cragg-Donald F present alongside KP stats sim_no_constant (HC1)", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- file.path(fixture_dir, "sim_no_constant_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})
