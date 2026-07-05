# ============================================================================
# Tests: One-way cluster-robust VCV (Ticket C2)
# ============================================================================
#
# Family M-13 (re-based 2026-07-05): IV 2SLS x one-way cluster Stata parity
# on the abdata ab_cl cells (H88 model minus gmm2s, help.txt:1541;
# cluster(id) anchor H91, help.txt:1558; D5a). OLS x one-way cluster Stata
# parity is owned by the hf suite (H91 abdata, H103 nlswork -- compare_hf
# asserts coef/vcov/all diagnostics); the griliches cluster(year) M=7
# known-dirty warning-lesson is owned by hf H27/H28.
#
# Fixture naming: "cl" was generated with ivreg2's `cluster()` option (no
# small-sample correction); "cl_small" was generated with `cluster() small`
# (with the (N-1)/(N-K) * M/(M-1) correction).

data(abdata, package = "ivreg2r")

# Shared fits (M-17 hoisting precedent): each configuration is fit once here
# and reused across the parity, metadata, and property sections below.
fit_cl       <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                        clusters = ~id, endog = "w")
fit_cl_small <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                        clusters = ~id, small = TRUE, endog = "w")
fit_iid      <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                        endog = "w")
fit_hc1      <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                        vcov = "robust", endog = "w")


# ============================================================================
# Stata parity: cluster VCV and coefficients/SEs (small = FALSE)
# ============================================================================

test_that("2SLS cluster VCV matches Stata ab_cl cl fixture", {
  skip_if(!file.exists(fixture_path("ab_cl_vcov_cl.csv")), "VCV fixture not found")
  expect_vcov_fixture(fit_cl, "ab_cl_vcov_cl.csv")
})

test_that("2SLS cluster coef/SE match Stata ab_cl cl fixture", {
  skip_if(!file.exists(fixture_path("ab_cl_coef_cl.csv")), "Coef fixture not found")
  expect_coef_fixture(fit_cl, "ab_cl_coef_cl.csv")
})


# ============================================================================
# Stata parity: cluster VCV and coefficients/SEs (small = TRUE)
# ============================================================================

test_that("2SLS cluster VCV matches Stata ab_cl cl_small fixture", {
  skip_if(!file.exists(fixture_path("ab_cl_vcov_cl_small.csv")), "VCV fixture not found")
  expect_vcov_fixture(fit_cl_small, "ab_cl_vcov_cl_small.csv")
})

test_that("2SLS cluster coef/SE match Stata ab_cl cl_small fixture", {
  skip_if(!file.exists(fixture_path("ab_cl_coef_cl_small.csv")), "Coef fixture not found")
  expect_coef_fixture(fit_cl_small, "ab_cl_coef_cl_small.csv")
})


# ============================================================================
# Metadata parity
# ============================================================================

test_that("n_clusters is 140 for ab_cl data", {
  expect_identical(fit_cl$n_clusters, 140L)
})

test_that("df.residual is M-1 when clustered", {
  expect_identical(fit_cl$df.residual, 139L)
})

test_that("vcov_type is 'CL' when clustered", {
  expect_identical(fit_cl$vcov_type, "CL")
})

test_that("cluster_var is 'id' when clustered", {
  expect_identical(fit_cl$cluster_var, "id")
})

test_that("N matches Stata ab_cl diagnostics fixture", {
  diag_path <- fixture_path("ab_cl_diagnostics_cl.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")
  dx <- read_diagnostics(diag_path)
  expect_identical(nobs(fit_cl), as.integer(dx$N))
})


# ============================================================================
# Property: cluster VCV differs from HC1 VCV
# ============================================================================

test_that("cluster VCV differs from HC1 VCV", {
  expect_false(isTRUE(all.equal(fit_cl$vcov, fit_hc1$vcov)))
})


# ============================================================================
# Property: small=TRUE cluster VCV = small=FALSE * correction factor
# ============================================================================

test_that("small=TRUE cluster VCV equals small=FALSE times correction factor", {
  N <- fit_cl$nobs
  K <- length(coef(fit_cl))
  M <- fit_cl$n_clusters
  correction <- ((N - 1) / (N - K)) * (M / (M - 1))
  expect_equal(fit_cl_small$vcov, fit_cl$vcov * correction,
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Property: coefficients identical across iid / HC1 / CL
# ============================================================================

test_that("coefficients are identical for iid, HC1, and CL", {
  expect_equal(coef(fit_hc1), coef(fit_iid),
               tolerance = .Machine$double.eps^0.5)
  expect_equal(coef(fit_cl), coef(fit_iid),
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Property: cluster VCV is symmetric
# ============================================================================

test_that("cluster VCV is symmetric", {
  expect_equal(fit_cl$vcov, t(fit_cl$vcov),
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Input validation: non-formula clusters → error
# ============================================================================

test_that("non-formula clusters arg gives error", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, clusters = "cyl"),
    "one-sided formula"
  )
})


# ============================================================================
# Input validation: multi-variable clusters → error
# ============================================================================

test_that("three-variable clusters formula gives error", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, clusters = ~ cyl + gear + carb),
    "one or two variables"
  )
})


# ============================================================================
# Input validation: missing cluster variable → error
# ============================================================================

test_that("missing cluster variable gives error", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, clusters = ~ nonexistent),
    "not found in data"
  )
})


# ============================================================================
# Input validation: NA in cluster variable → error
# ============================================================================

test_that("NA in cluster variable gives error", {
  df <- mtcars
  df$clust <- c(NA, rep(1:3, length.out = nrow(df) - 1))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = df, clusters = ~ clust),
    "contains NA"
  )
})


# ============================================================================
# Robustness: tibble input works
# ============================================================================

test_that("cluster VCV works with tibble data", {
  skip_if_not_installed("tibble")

  d <- tibble::as_tibble(mtcars)
  fit <- ivreg2(mpg ~ wt + hp, data = d, clusters = ~cyl)
  expect_identical(fit$vcov_type, "CL")
  expect_identical(fit$n_clusters, length(unique(mtcars$cyl)))
})


# ============================================================================
# Robustness: character rownames (mtcars) work
# ============================================================================

test_that("cluster VCV works with character rownames", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, clusters = ~cyl)
  expect_identical(fit$vcov_type, "CL")
  expect_identical(fit$n_clusters, length(unique(mtcars$cyl)))
})
