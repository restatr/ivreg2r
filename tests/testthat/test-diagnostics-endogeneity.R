# ============================================================================
# Tests: Endogeneity Test / C-Statistic (Ticket E4)
# ============================================================================
#
# Tests H0: specified endogenous regressors are actually exogenous.
# Difference-in-Sargan (IID) or difference-in-J (robust/cluster) statistic.
# Verifies against Stata ivreg2 fixtures (e(estat), e(estatp), e(estatdf)).

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

sim_nc_path <- file.path(fixture_dir, "sim_no_constant_data.csv")
if (file.exists(sim_nc_path)) {
  sim_nc <- read.csv(sim_nc_path)
}


# ============================================================================
# Structural tests (no fixtures)
# ============================================================================

test_that("endogeneity is NULL for OLS models", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$diagnostics$endogeneity)
})

test_that("endogeneity has all expected fields", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  endog <- fit$diagnostics$endogeneity
  expect_type(endog, "list")
  expect_named(endog, c("stat", "p", "df", "test_name", "tested_vars"),
               ignore.order = TRUE)
  expect_equal(endog$tested_vars, "educ")
  expect_identical(endog$df, 1L)
})

test_that("small=TRUE and small=FALSE produce identical endogeneity stats (IID)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, small = TRUE)
  expect_equal(fit1$diagnostics$endogeneity$stat,
               fit2$diagnostics$endogeneity$stat)
  expect_equal(fit1$diagnostics$endogeneity$p,
               fit2$diagnostics$endogeneity$p)
})

test_that("small=TRUE and small=FALSE produce identical endogeneity stats (HC1)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1", small = TRUE)
  expect_equal(fit1$diagnostics$endogeneity$stat,
               fit2$diagnostics$endogeneity$stat)
  expect_equal(fit1$diagnostics$endogeneity$p,
               fit2$diagnostics$endogeneity$p)
})

test_that("HC0 and HC1 produce identical endogeneity stat", {
  skip_if(!file.exists(card_path), "card data not found")
  fit0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC0")
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = "HC1")
  expect_equal(fit0$diagnostics$endogeneity$stat,
               fit1$diagnostics$endogeneity$stat)
})

test_that("invalid endog variable name produces error", {
  skip_if(!file.exists(card_path), "card data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
           data = card, endog = "not_a_var"),
    "not in the endogenous list"
  )
})

test_that("endog = endo_names equals default endog = NULL", {
  skip_if(!file.exists(card_path), "card data not found")
  fit_default <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                        data = card)
  fit_explicit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                         data = card, endog = "educ")
  expect_equal(fit_default$diagnostics$endogeneity$stat,
               fit_explicit$diagnostics$endogeneity$stat)
  expect_equal(fit_default$diagnostics$endogeneity$p,
               fit_explicit$diagnostics$endogeneity$p)
})

test_that("endog must be character or NULL", {
  skip_if(!file.exists(card_path), "card data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
           data = card, endog = 1),
    "character vector or NULL"
  )
})

test_that("duplicate endog entries are deduplicated with warning", {
  skip_if(!file.exists(card_path), "card data not found")
  fit_single <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                       data = card, endog = "educ")
  expect_warning(
    fit_dup <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                      data = card, endog = c("educ", "educ")),
    "duplicate entries"
  )
  expect_equal(fit_dup$diagnostics$endogeneity$stat,
               fit_single$diagnostics$endogeneity$stat)
  expect_equal(fit_dup$diagnostics$endogeneity$df,
               fit_single$diagnostics$endogeneity$df)
})


# ============================================================================
# Helper: run fixture comparison for one spec/VCE combination
# ============================================================================

check_endogeneity <- function(fit, fixture_path, label) {
  fixture <- read_diagnostics(fixture_path)
  endog <- fit$diagnostics$endogeneity

  # Skip if Stata didn't compute the stat (e.g., too few clusters)
  if (is.na(fixture$estat)) return(invisible())

  test_that(paste("Endogeneity stat matches Stata", label), {
    expect_equal(endog$stat, fixture$estat,
                 tolerance = stata_tol$stat)
  })

  test_that(paste("Endogeneity p-value matches Stata", label), {
    expect_equal(endog$p, fixture$estatp,
                 tolerance = stata_tol$pval)
  })

  test_that(paste("Endogeneity df matches Stata", label), {
    expect_identical(endog$df, as.integer(fixture$estatdf))
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
    check_endogeneity(fit, fixture_file, label)
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
    check_endogeneity(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_cluster: y ~ x1 | endo1 | z1+z2, clusters=~cluster_id
# M=50
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
    check_endogeneity(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_multi_endo: y ~ x1+x2 | endo1+endo2 | z1+z2+z3+z4
# K1=2, L1=4, tests both endogenous
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
    check_endogeneity(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_no_constant: y ~ 0 + x1 | endo1 | z1+z2 (noconstant)
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "HC1",  small = FALSE, suffix = "hc1")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_no_constant_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_no_constant", vce_combo$suffix)

  if (file.exists(sim_nc_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2,
                  data = sim_nc, vcov = vce_combo$vcov,
                  small = vce_combo$small)
    check_endogeneity(fit, fixture_file, label)
  }
}


# ============================================================================
# card_just_id_weighted: lwage ~ exper+expersq+black+south | educ | nearc4
# with weights=weight
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
    check_endogeneity(fit, fixture_file, label)
  }
}
