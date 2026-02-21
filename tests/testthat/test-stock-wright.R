# ============================================================================
# Tests: Stock-Wright LM S Statistic (Ticket J2)
# ============================================================================
#
# Weak-instrument-robust LM test of H0: all endogenous coefficients = 0
# and orthogonality conditions are valid. Verifies against Stata ivreg2
# fixtures (sstat, sstatp, sstatdf columns).

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

test_that("stock_wright is NULL for OLS models", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$diagnostics$stock_wright)
})

test_that("stock_wright has all expected fields", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  sw <- fit$diagnostics$stock_wright
  expect_type(sw, "list")
  expect_named(sw, c("stat", "p", "df"), ignore.order = TRUE)
})

test_that("stock_wright df = L1 (number of excluded instruments)", {
  skip_if(!file.exists(card_path), "card data not found")
  # Just identified: L1 = 1
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card)
  expect_identical(fit1$diagnostics$stock_wright$df, 1L)

  # Overidentified: L1 = 2
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                 data = card)
  expect_identical(fit2$diagnostics$stock_wright$df, 2L)
})

test_that("small=TRUE and small=FALSE produce identical S stats (IID)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = TRUE)
  expect_equal(fit1$diagnostics$stock_wright$stat,
               fit2$diagnostics$stock_wright$stat)
  expect_equal(fit1$diagnostics$stock_wright$p,
               fit2$diagnostics$stock_wright$p)
})

# S stat differs between IID and HC because the omega formula changes
test_that("IID and HC1 produce different S stats", {
  skip_if(!file.exists(card_path), "card data not found")
  fit_iid <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card, vcov = "iid")
  fit_hc1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card, vcov = "HC1")
  # These should differ (different omega formulas)
  expect_false(isTRUE(all.equal(
    fit_iid$diagnostics$stock_wright$stat,
    fit_hc1$diagnostics$stock_wright$stat
  )))
})

# HC0 and HC1 should produce identical S stats (no small-sample correction
# in the omega — both use Z' diag(e^2) Z / (N - dofminus))
test_that("HC0 and HC1 produce identical S stats", {
  skip_if(!file.exists(card_path), "card data not found")
  fit0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC0")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1")
  expect_equal(fit0$diagnostics$stock_wright$stat,
               fit1$diagnostics$stock_wright$stat)
  expect_equal(fit0$diagnostics$stock_wright$p,
               fit1$diagnostics$stock_wright$p)
})

test_that("stock_wright appears in glance()", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  g <- glance(fit)
  expect_true("stock_wright_stat" %in% names(g))
  expect_true("stock_wright_p" %in% names(g))
  expect_true("stock_wright_df" %in% names(g))
  expect_equal(g$stock_wright_stat, fit$diagnostics$stock_wright$stat)
})

test_that("glance() stock_wright columns are NA for OLS", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  g <- glance(fit)
  expect_true(is.na(g$stock_wright_stat))
  expect_true(is.na(g$stock_wright_p))
  expect_true(is.na(g$stock_wright_df))
})


# ============================================================================
# Helper: run fixture comparison for one spec/VCE combination
# ============================================================================

check_stock_wright <- function(fit, fixture_path, label) {
  fixture <- read_diagnostics(fixture_path)
  sw <- fit$diagnostics$stock_wright

  # Stata may leave sstat empty when omega is rank-deficient (e.g., M=2)
  fixture_missing <- is.na(fixture$sstat) || identical(fixture$sstat, "")

  if (fixture_missing) {
    test_that(paste("S stat is NA when Stata omits it", label), {
      expect_true(is.na(sw$stat))
    })
  } else {
    test_that(paste("S stat matches Stata", label), {
      expect_equal(sw$stat, fixture$sstat,
                   tolerance = stata_tol$stat)
    })

    test_that(paste("S p-value matches Stata", label), {
      expect_equal(sw$p, fixture$sstatp,
                   tolerance = stata_tol$pval)
    })

    test_that(paste("S df matches Stata", label), {
      expect_identical(sw$df, as.integer(fixture$sstatdf))
    })
  }
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
    check_stock_wright(fit, fixture_file, label)
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
    check_stock_wright(fit, fixture_file, label)
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
    check_stock_wright(fit, fixture_file, label)
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
    check_stock_wright(fit, fixture_file, label)
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
    check_stock_wright(fit, fixture_file, label)
  }
}


# ============================================================================
# dofminus fixtures: dofminus=1, sdofminus=1
# ============================================================================

# card_just_id with dofminus
for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid",      clusters = NULL),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small", clusters = NULL),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1",      clusters = NULL),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small", clusters = NULL),
  list(vcov = "iid",  small = FALSE, suffix = "cl",       clusters = ~smsa66),
  list(vcov = "iid",  small = TRUE,  suffix = "cl_small", clusters = ~smsa66)
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("card_just_id_dofminus_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_just_id_dofminus", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card, vcov = vce_combo$vcov, small = vce_combo$small,
                  clusters = vce_combo$clusters,
                  dofminus = 1L, sdofminus = 1L)
    check_stock_wright(fit, fixture_file, label)
  }
}

# card_overid with dofminus
for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid",      clusters = NULL),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small", clusters = NULL),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1",      clusters = NULL),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small", clusters = NULL),
  list(vcov = "iid",  small = FALSE, suffix = "cl",       clusters = ~smsa66),
  list(vcov = "iid",  small = TRUE,  suffix = "cl_small", clusters = ~smsa66)
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("card_overid_dofminus_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_overid_dofminus", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                  data = card, vcov = vce_combo$vcov, small = vce_combo$small,
                  clusters = vce_combo$clusters,
                  dofminus = 1L, sdofminus = 1L)
    check_stock_wright(fit, fixture_file, label)
  }
}

# sim_cluster with dofminus
for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid",      clusters = NULL),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small", clusters = NULL),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1",      clusters = NULL),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small", clusters = NULL),
  list(vcov = "iid",  small = FALSE, suffix = "cl",       clusters = ~cluster_id),
  list(vcov = "iid",  small = TRUE,  suffix = "cl_small", clusters = ~cluster_id)
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_cluster_dofminus_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_cluster_dofminus", vce_combo$suffix)

  if (file.exists(sim_cluster_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                  data = sim_cluster, vcov = vce_combo$vcov,
                  small = vce_combo$small, clusters = vce_combo$clusters,
                  dofminus = 1L, sdofminus = 1L)
    check_stock_wright(fit, fixture_file, label)
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
    check_stock_wright(fit, fixture_file, label)
  }
}
