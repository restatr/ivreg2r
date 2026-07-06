# ============================================================================
# Tests: Stock-Wright LM S Statistic (Ticket J2)
# ============================================================================
#
# Weak-instrument-robust LM test of H0: all endogenous coefficients = 0
# and orthogonality conditions are valid. Verifies against Stata ivreg2
# fixtures (sstat, sstatp, sstatdf columns).

# --- Load datasets ---
# card_data.csv is retained here ONLY for the dofminus block below (family
# M-29, not yet re-based); the card_just_id/card_overid diagnostics blocks
# that used to read it were retired to mroz at family M-10 (2026-07-06).
card_path <- fixture_path("card_data.csv")
if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

sim_multi_path <- fixture_path("sim_multi_endo_data.csv")
if (file.exists(sim_multi_path)) {
  sim_multi <- read.csv(sim_multi_path)
}

sim_cluster_path <- fixture_path("sim_cluster_data.csv")
if (file.exists(sim_cluster_path)) {
  sim_cluster <- read.csv(sim_cluster_path)
}

data(abdata, package = "ivreg2r")
data(mroz, package = "ivreg2r")


# ============================================================================
# Structural tests (no fixtures)
# ============================================================================

test_that("stock_wright is NULL for OLS models", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$diagnostics$stock_wright)
})

test_that("stock_wright has all expected fields", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz)
  sw <- fit$diagnostics$stock_wright
  expect_type(sw, "list")
  expect_named(sw, c("stat", "p", "df"), ignore.order = TRUE)
})

test_that("stock_wright df = L1 (number of excluded instruments)", {
  # Just identified: L1 = 1
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz)
  expect_identical(fit1$diagnostics$stock_wright$df, 1L)

  # Overidentified: L1 = 3
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                 data = mroz)
  expect_identical(fit2$diagnostics$stock_wright$df, 3L)
})

test_that("small=TRUE and small=FALSE produce identical S stats (IID)", {
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = TRUE)
  expect_equal(fit1$diagnostics$stock_wright$stat,
               fit2$diagnostics$stock_wright$stat)
  expect_equal(fit1$diagnostics$stock_wright$p,
               fit2$diagnostics$stock_wright$p)
})

# S stat differs between IID and HC because the omega formula changes
test_that("IID and HC1 produce different S stats", {
  fit_iid <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, vcov = "iid")
  fit_hc1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, vcov = "robust")
  # These should differ (different omega formulas)
  expect_false(isTRUE(all.equal(
    fit_iid$diagnostics$stock_wright$stat,
    fit_hc1$diagnostics$stock_wright$stat
  )))
})

# HC0 and HC1 should produce identical S stats (no small-sample correction
# in the omega — both use Z' diag(e^2) Z / (N - dofminus))
test_that("HC0 and HC1 produce identical S stats", {
  fit0 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, vcov = "robust")
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, vcov = "robust")
  expect_equal(fit0$diagnostics$stock_wright$stat,
               fit1$diagnostics$stock_wright$stat)
  expect_equal(fit0$diagnostics$stock_wright$p,
               fit1$diagnostics$stock_wright$p)
})

test_that("stock_wright appears in glance()", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz)
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
# mroz_justid: lwage ~ exper+expersq | educ | age (D5a disclosed
# single-instrument restriction; family M-10 re-base, replacing card_just_id)
# K1=1, L1=1, just-identified
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",    suffix = "iid"),
  list(vcov = "robust", suffix = "robust")
)) {
  fixture_file <- fixture_path(
    paste0("mroz_justid_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("mroz_justid", vce_combo$suffix)

  if (file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq | educ | age,
                  data = mroz, vcov = vce_combo$vcov)
    fit_small <- ivreg2(lwage ~ exper + expersq | educ | age,
                         data = mroz, vcov = vce_combo$vcov, small = TRUE)
    check_stock_wright(fit, fixture_file, label)
    check_stock_wright(fit_small, fixture_file, paste(label, "small"))

    test_that(paste("Stock-Wright direct equality: small vs non-small", label), {
      expect_equal(fit_small$diagnostics$stock_wright$stat,
                   fit$diagnostics$stock_wright$stat)
      expect_equal(fit_small$diagnostics$stock_wright$p,
                   fit$diagnostics$stock_wright$p)
      expect_identical(fit_small$diagnostics$stock_wright$df,
                        fit$diagnostics$stock_wright$df)
    })
  }
}


# ============================================================================
# mroz overid: lwage ~ exper+expersq | educ | age+kidslt6+kidsge6.
# iid anchor = H36 (help.txt:1305, same fit as H31 but with sstat populated,
# unlike H31); robust anchor = H41 (help.txt:1325). Family M-10 re-base,
# replacing card_overid.
# K1=1, L1=3, overidentified
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",    suffix = "H36"),
  list(vcov = "robust", suffix = "H41")
)) {
  fixture_file <- fixture_path(
    paste0("hf_mroz_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("mroz_overid", vce_combo$suffix)

  if (file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                  data = mroz, vcov = vce_combo$vcov)
    fit_small <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                         data = mroz, vcov = vce_combo$vcov, small = TRUE)
    check_stock_wright(fit, fixture_file, label)
    check_stock_wright(fit_small, fixture_file, paste(label, "small"))

    test_that(paste("Stock-Wright direct equality: small vs non-small", label), {
      expect_equal(fit_small$diagnostics$stock_wright$stat,
                   fit$diagnostics$stock_wright$stat)
      expect_equal(fit_small$diagnostics$stock_wright$p,
                   fit$diagnostics$stock_wright$p)
      expect_identical(fit_small$diagnostics$stock_wright$df,
                        fit$diagnostics$stock_wright$df)
    })
  }
}


# ============================================================================
# ab_cl: abdata H88-minus-gmm2s 2SLS, cluster(id) (family M-13, re-based;
# the iid/hc1 rows were retired -- unweighted iid/hc1 parity for these
# statistics is owned by the card blocks above, family M-10)
# ============================================================================

for (vce_combo in ab_cl_vce_cells) {
  fixture_file <- fixture_path(
    paste0("ab_cl_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("ab_cl", vce_combo$suffix)

  if (file.exists(fixture_file)) {
    fit <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                  clusters = ~id, small = vce_combo$small, endog = "w")
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
  list(vcov = "robust",  small = FALSE, suffix = "hc1"),
  list(vcov = "robust",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- fixture_path(
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
# card_fw: lwage ~ exper+expersq+black+south | educ | nearc4+nearc2,
# fweight fwt = mod(age,5)+1, cluster(region) M=9 (family M-11, re-based;
# card_wt from helper-fixtures.R)
# ============================================================================

for (vce_combo in card_fw_vce_cells) {
  fixture_file <- fixture_path(paste0("card_fw_diagnostics_", vce_combo$suffix, ".csv"))
  label <- paste("card_fw", vce_combo$suffix)

  if (file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                  data = card_wt, weights = fwt, weight_type = "fweight",
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
  list(vcov = "robust",  small = FALSE, suffix = "hc1",      clusters = NULL),
  list(vcov = "robust",  small = TRUE,  suffix = "hc1_small", clusters = NULL),
  list(vcov = "iid",  small = FALSE, suffix = "cl",       clusters = ~smsa66),
  list(vcov = "iid",  small = TRUE,  suffix = "cl_small", clusters = ~smsa66)
)) {
  fixture_file <- fixture_path(
    paste0("card_just_id_dofminus_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_just_id_dofminus", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    # M=2 clusters → expected rank-deficient diagnostics
    fit <- muffle_rank_warnings(
      ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
             data = card, vcov = vce_combo$vcov, small = vce_combo$small,
             clusters = vce_combo$clusters,
             dofminus = 1L, sdofminus = 1L)
    )
    check_stock_wright(fit, fixture_file, label)
  }
}

# card_overid with dofminus
for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid",      clusters = NULL),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small", clusters = NULL),
  list(vcov = "robust",  small = FALSE, suffix = "hc1",      clusters = NULL),
  list(vcov = "robust",  small = TRUE,  suffix = "hc1_small", clusters = NULL),
  list(vcov = "iid",  small = FALSE, suffix = "cl",       clusters = ~smsa66),
  list(vcov = "iid",  small = TRUE,  suffix = "cl_small", clusters = ~smsa66)
)) {
  fixture_file <- fixture_path(
    paste0("card_overid_dofminus_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_overid_dofminus", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    # M=2 clusters → expected rank-deficient diagnostics
    fit <- muffle_rank_warnings(
      ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
             data = card, vcov = vce_combo$vcov, small = vce_combo$small,
             clusters = vce_combo$clusters,
             dofminus = 1L, sdofminus = 1L)
    )
    check_stock_wright(fit, fixture_file, label)
  }
}

# sim_cluster with dofminus
for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid",      clusters = NULL),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small", clusters = NULL),
  list(vcov = "robust",  small = FALSE, suffix = "hc1",      clusters = NULL),
  list(vcov = "robust",  small = TRUE,  suffix = "hc1_small", clusters = NULL),
  list(vcov = "iid",  small = FALSE, suffix = "cl",       clusters = ~cluster_id),
  list(vcov = "iid",  small = TRUE,  suffix = "cl_small", clusters = ~cluster_id)
)) {
  fixture_file <- fixture_path(
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

sim_noconst_path <- fixture_path("sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "robust",  small = FALSE, suffix = "hc1"),
  list(vcov = "robust",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- fixture_path(
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
