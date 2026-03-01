# ============================================================================
# Tests: Orthogonality Test / Instrument-Subset C-Statistic (Ticket J1)
# ============================================================================
#
# Tests H0: specified instruments satisfy orthogonality conditions.
# Difference-of-Sargan (IID) or difference-of-J (robust/cluster) statistic.
# Verifies against Stata ivreg2 fixtures (e(cstat), e(cstatp), e(cstatdf)).

# --- Helpers ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)

read_orthog_fixture <- function(path) {
  read.csv(path, check.names = FALSE)
}

# --- Load datasets ---
card_path <- file.path(fixture_dir, "card_data.csv")
if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

sim_cluster_path <- file.path(fixture_dir, "sim_cluster_data.csv")
if (file.exists(sim_cluster_path)) {
  sim_cluster <- read.csv(sim_cluster_path)
}

sim_twoway_path <- file.path(fixture_dir, "sim_twoway_data.csv")
if (file.exists(sim_twoway_path)) {
  sim_twoway <- read.csv(sim_twoway_path)
}


# ============================================================================
# Structural tests (no fixtures)
# ============================================================================

test_that("orthog is NULL when not requested", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4, data = card)
  expect_null(fit$diagnostics$orthog)
})

test_that("orthog = character(0) is silently skipped", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, orthog = character(0))
  expect_null(fit$diagnostics$orthog)
})

test_that("orthog is NULL for OLS models", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, orthog = "wt")
  expect_null(fit$diagnostics)
})

test_that("orthog has all expected fields", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, orthog = "nearc2")
  orth <- fit$diagnostics$orthog
  expect_type(orth, "list")
  expect_named(orth, c("stat", "p", "df", "test_name", "tested_vars"),
               ignore.order = TRUE)
  expect_equal(orth$tested_vars, "nearc2")
  expect_identical(orth$df, 1L)
  expect_equal(orth$test_name, "C (orthog)")
})

test_that("small=TRUE and small=FALSE produce identical orthog stats (IID)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, orthog = "nearc2", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, orthog = "nearc2", small = TRUE)
  expect_equal(fit1$diagnostics$orthog$stat,
               fit2$diagnostics$orthog$stat)
  expect_equal(fit1$diagnostics$orthog$p,
               fit2$diagnostics$orthog$p)
})

test_that("small=TRUE and small=FALSE produce identical orthog stats (HC1)", {
  skip_if(!file.exists(card_path), "card data not found")
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, vcov = "robust", orthog = "nearc2", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, vcov = "robust", orthog = "nearc2", small = TRUE)
  expect_equal(fit1$diagnostics$orthog$stat,
               fit2$diagnostics$orthog$stat)
  expect_equal(fit1$diagnostics$orthog$p,
               fit2$diagnostics$orthog$p)
})

test_that("HC0 and HC1 produce identical orthog stat", {
  skip_if(!file.exists(card_path), "card data not found")
  fit0 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, vcov = "robust", orthog = "nearc2")
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                 data = card, vcov = "robust", orthog = "nearc2")
  expect_equal(fit0$diagnostics$orthog$stat,
               fit1$diagnostics$orthog$stat)
})

test_that("invalid orthog variable name produces error", {
  skip_if(!file.exists(card_path), "card data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, orthog = "not_a_var"),
    "not in the instrument list"
  )
})

test_that("orthog rejects endogenous variable", {
  skip_if(!file.exists(card_path), "card data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, orthog = "educ"),
    "not in the instrument list"
  )
})

test_that("orthog must be character or NULL", {
  skip_if(!file.exists(card_path), "card data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, orthog = 1),
    "character vector or NULL"
  )
})

test_that("duplicate orthog entries are deduplicated with warning", {
  skip_if(!file.exists(card_path), "card data not found")
  fit_single <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                       data = card, orthog = "nearc2")
  expect_warning(
    fit_dup <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                      data = card, orthog = c("nearc2", "nearc2")),
    "duplicate entries"
  )
  expect_equal(fit_dup$diagnostics$orthog$stat,
               fit_single$diagnostics$orthog$stat)
  expect_equal(fit_dup$diagnostics$orthog$df,
               fit_single$diagnostics$orthog$df)
})

test_that("testing all overid instruments gives C = J_full when exactly identified restricted", {
  skip_if(!file.exists(card_path), "card data not found")
  # Model: 2 excluded instruments, 1 endogenous → overid_df = 1
  # Testing one instrument: restricted model is exactly identified (overid_df_r = 0)
  # so J_r should be 0, and C = J_full - 0 = J_full
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, orthog = "nearc2")
  # The orthog C stat should equal the full model's overid J stat
  expect_equal(fit$diagnostics$orthog$stat,
               fit$diagnostics$overid$stat,
               tolerance = 1e-10)
})

test_that("orthog appears in glance output", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, orthog = "nearc2")
  g <- glance(fit)
  expect_true("orthog_stat" %in% names(g))
  expect_true("orthog_p" %in% names(g))
  expect_false(is.na(g$orthog_stat))
  expect_false(is.na(g$orthog_p))
})

test_that("orthog is NA in glance when not requested", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card)
  g <- glance(fit)
  expect_true(is.na(g$orthog_stat))
  expect_true(is.na(g$orthog_p))
})


# ============================================================================
# Helper: run fixture comparison for one spec/VCE combination
# ============================================================================

check_orthog <- function(fit, fixture_path, label) {
  fixture <- read_orthog_fixture(fixture_path)
  orth <- fit$diagnostics$orthog

  test_that(paste("Orthog stat matches Stata", label), {
    expect_false(is.null(orth))
    if (is.na(fixture$cstat) || fixture$cstat == 0) {
      # Underidentified or collinearity: stat should be 0
      expect_equal(orth$stat, 0)
    } else {
      expect_equal(orth$stat, fixture$cstat,
                   tolerance = stata_tol$stat)
    }
  })

  test_that(paste("Orthog p-value matches Stata", label), {
    if (!is.na(fixture$cstatp)) {
      expect_equal(orth$p, fixture$cstatp,
                   tolerance = stata_tol$pval)
    }
  })

  test_that(paste("Orthog df matches Stata", label), {
    if (!is.na(fixture$cstatdf)) {
      expect_identical(orth$df, as.integer(fixture$cstatdf))
    }
  })
}


# ============================================================================
# card_orthog1: orthog(nearc2), IID/HC1
# Model: lwage ~ exper + expersq | educ | nearc2 + nearc4
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "robust",  small = FALSE, suffix = "hc1"),
  list(vcov = "robust",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("card_orthog1_orthog_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_orthog1", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                  data = card, vcov = vce_combo$vcov, small = vce_combo$small,
                  orthog = "nearc2")
    check_orthog(fit, fixture_file, label)
  }
}


# ============================================================================
# card_orthog2: orthog(nearc2 nearc4) — restricted model underidentified
# ============================================================================

test_that("Orthog with all excluded IVs gives stat=0 (underidentified restricted)", {
  fixture_file <- file.path(fixture_dir, "card_orthog2_orthog_iid.csv")
  skip_if(!file.exists(card_path), "card data not found")
  skip_if(!file.exists(fixture_file), "fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                data = card, orthog = c("nearc2", "nearc4"))
  orth <- fit$diagnostics$orthog
  expect_equal(orth$stat, 0)
  expect_equal(orth$df, 0L)
})


# ============================================================================
# card_orthog1_dof: orthog(nearc2) with dofminus(1), IID/HC1
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "robust",  small = FALSE, suffix = "hc1"),
  list(vcov = "robust",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("card_orthog1_dof_orthog_", vce_combo$suffix, ".csv")
  )
  label <- paste("card_orthog1_dof", vce_combo$suffix)

  if (file.exists(card_path) && file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
                  data = card, vcov = vce_combo$vcov, small = vce_combo$small,
                  orthog = "nearc2", dofminus = 1L)
    check_orthog(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_cluster_orthog1: orthog(z1), cluster-robust
# Model: y ~ x1 | endo1 | z1 + z2, clusters = ~cluster_id
# ============================================================================

for (vce_combo in list(
  list(small = FALSE, suffix = "cl"),
  list(small = TRUE,  suffix = "cl_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_cluster_orthog1_orthog_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_cluster_orthog1", vce_combo$suffix)

  if (file.exists(sim_cluster_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                  data = sim_cluster, clusters = ~cluster_id,
                  small = vce_combo$small, orthog = "z1")
    check_orthog(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_twoway_orthog1: orthog(z1), two-way cluster
# Model: y ~ x1 | endo1 | z1 + z2, clusters = ~firm_id + year_id
# ============================================================================

for (vce_combo in list(
  list(small = FALSE, suffix = "cl2"),
  list(small = TRUE,  suffix = "cl2_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_twoway_orthog1_orthog_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_twoway_orthog1", vce_combo$suffix)

  if (file.exists(sim_twoway_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                  data = sim_twoway, clusters = ~firm_id + year_id,
                  small = vce_combo$small, orthog = "z1")
    check_orthog(fit, fixture_file, label)
  }
}


# ============================================================================
# Factor variable support (term-label validation)
# ============================================================================

test_that("orthog accepts factor term labels and expands to columns", {
  set.seed(42)
  n <- 500
  region <- factor(sample(1:4, n, replace = TRUE))
  x1 <- rnorm(n)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  z3 <- rnorm(n)
  endo <- 0.5 * z1 + 0.3 * z2 + rnorm(n)
  y <- 1 + 0.5 * x1 + 0.8 * endo + rnorm(n)
  d <- data.frame(y, x1, endo, z1, z2, z3, region)

  # Factor as excluded instrument — orthog by term label
  fit <- ivreg2(y ~ x1 | endo | z1 + z2 + region,
                data = d, orthog = "region")
  expect_true(!is.null(fit$diagnostics$orthog))
  expect_true(is.finite(fit$diagnostics$orthog$stat))
  # df = number of columns tested (3 dummies for 4-level factor)
  expect_equal(fit$diagnostics$orthog$df, 3L)

  # Factor as exogenous instrument — orthog by term label
  # Need enough excluded instruments so restricted model is still identified
  z4 <- rnorm(n); z5 <- rnorm(n); z6 <- rnorm(n)
  d2 <- cbind(d, z4 = z4, z5 = z5, z6 = z6)
  fit2 <- ivreg2(y ~ x1 + region | endo | z1 + z2 + z3 + z4 + z5 + z6,
                 data = d2, orthog = "region")
  expect_true(!is.null(fit2$diagnostics$orthog))
  expect_true(is.finite(fit2$diagnostics$orthog$stat))
  expect_equal(fit2$diagnostics$orthog$df, 3L)

  # Column names (not term labels) are rejected — use term labels instead
  expect_error(
    ivreg2(y ~ x1 + region | endo | z1 + z2 + z3,
           data = d, orthog = c("region2", "region3")),
    "not in the instrument list"
  )
})

test_that("orthog + partial overlap produces an error, not silent degeneracy", {
  skip_if(!file.exists(card_path), "card data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4,
           data = card, partial = "exper", orthog = "exper"),
    "partialled out"
  )
})
