# ============================================================================
# Tests: One-way cluster-robust VCV (Ticket C2)
# ============================================================================
#
# Fixture naming: Stata fixtures named "cl" were generated with ivreg2's
# `cluster()` option (no small-sample correction). Fixtures named "cl_small"
# were generated with `cluster() small` (with (N-1)/(N-K) * M/(M-1) correction).

# --- Helpers ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)

read_vcov_fixture <- function(path) {
  fixture <- read.csv(path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)
  vcov_cols <- grep("^vcov_", names(fixture), value = TRUE)
  V <- as.matrix(fixture[, vcov_cols])
  rownames(V) <- r_names
  col_stata <- sub("^vcov_", "", vcov_cols)
  colnames(V) <- ifelse(col_stata == "_cons", "(Intercept)", col_stata)
  V
}

expect_vcov_equal <- function(V_r, V_stata, tol = stata_tol$vcov) {
  shared <- intersect(rownames(V_r), rownames(V_stata))
  for (rn in shared) {
    for (cn in shared) {
      expect_equal(
        V_r[rn, cn], V_stata[rn, cn],
        tolerance = tol,
        info = paste("VCV mismatch:", rn, cn)
      )
    }
  }
}

# --- Load sim_cluster data ---
sim_path <- file.path(fixture_dir, "sim_cluster_data.csv")
if (file.exists(sim_path)) {
  sim_cluster <- read.csv(sim_path)
}


# ============================================================================
# Stata parity: cluster VCV (small=FALSE)
# ============================================================================

test_that("2SLS cluster VCV matches Stata sim_cluster cl fixture", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")
  vcov_path <- file.path(fixture_dir, "sim_cluster_vcov_cl.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})


# ============================================================================
# Stata parity: cluster VCV (small=TRUE)
# ============================================================================

test_that("2SLS cluster VCV matches Stata sim_cluster cl_small fixture", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")
  vcov_path <- file.path(fixture_dir, "sim_cluster_vcov_cl_small.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id, small = TRUE)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})


# ============================================================================
# Stata parity: cluster SEs (small=FALSE)
# ============================================================================

test_that("2SLS cluster SEs match Stata sim_cluster cl fixture", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")
  coef_path <- file.path(fixture_dir, "sim_cluster_coef_cl.csv")
  skip_if(!file.exists(coef_path), "Coef fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id)
  fixture <- read.csv(coef_path)
  fixture$r_name <- ifelse(fixture$term == "_cons", "(Intercept)", fixture$term)

  for (i in seq_len(nrow(fixture))) {
    nm <- fixture$r_name[i]
    expect_equal(
      sqrt(fit$vcov[nm, nm]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", nm)
    )
  }
})


# ============================================================================
# Stata parity: cluster SEs (small=TRUE)
# ============================================================================

test_that("2SLS cluster SEs match Stata sim_cluster cl_small fixture", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")
  coef_path <- file.path(fixture_dir, "sim_cluster_coef_cl_small.csv")
  skip_if(!file.exists(coef_path), "Coef fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id, small = TRUE)
  fixture <- read.csv(coef_path)
  fixture$r_name <- ifelse(fixture$term == "_cons", "(Intercept)", fixture$term)

  for (i in seq_len(nrow(fixture))) {
    nm <- fixture$r_name[i]
    expect_equal(
      sqrt(fit$vcov[nm, nm]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", nm)
    )
  }
})


# ============================================================================
# Stata parity: n_clusters = 50
# ============================================================================

test_that("n_clusters is 50 for sim_cluster data", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id)
  expect_identical(fit$n_clusters, 50L)
})


# ============================================================================
# Stata parity: df.residual = M - 1
# ============================================================================

test_that("df.residual is M-1 when clustered", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id)
  expect_identical(fit$df.residual, 49L)
})


# ============================================================================
# Property: cluster VCV differs from HC1 VCV
# ============================================================================

test_that("cluster VCV differs from HC1 VCV", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")

  fit_cl <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                   data = sim_cluster, clusters = ~cluster_id)
  fit_hc <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                   data = sim_cluster, vcov = "HC1")
  expect_false(isTRUE(all.equal(fit_cl$vcov, fit_hc$vcov)))
})


# ============================================================================
# Property: small=TRUE cluster VCV = small=FALSE * correction factor
# ============================================================================

test_that("small=TRUE cluster VCV equals small=FALSE times correction factor", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")

  fit_no <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                   data = sim_cluster, clusters = ~cluster_id, small = FALSE)
  fit_sm <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                   data = sim_cluster, clusters = ~cluster_id, small = TRUE)
  N <- fit_no$nobs
  K <- length(coef(fit_no))
  M <- fit_no$n_clusters
  correction <- ((N - 1) / (N - K)) * (M / (M - 1))
  expect_equal(fit_sm$vcov, fit_no$vcov * correction,
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Property: coefficients identical across iid / HC1 / CL
# ============================================================================

test_that("coefficients are identical for iid, HC1, and CL", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")

  fit_iid <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                    data = sim_cluster, vcov = "iid")
  fit_hc  <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                    data = sim_cluster, vcov = "HC1")
  fit_cl  <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                    data = sim_cluster, clusters = ~cluster_id)
  expect_equal(coef(fit_hc), coef(fit_iid),
               tolerance = .Machine$double.eps^0.5)
  expect_equal(coef(fit_cl), coef(fit_iid),
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Property: cluster VCV is symmetric
# ============================================================================

test_that("cluster VCV is symmetric", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id)
  expect_equal(fit$vcov, t(fit$vcov),
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Property: vcov_type is "CL" when clustered
# ============================================================================

test_that("vcov_type is 'CL' when clustered", {
  skip_if(!file.exists(sim_path), "sim_cluster data not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_cluster, clusters = ~cluster_id)
  expect_identical(fit$vcov_type, "CL")
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
