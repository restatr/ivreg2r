# ============================================================================
# Tests: Model F-test (Ticket F1)
# ============================================================================
#
# Wald F-test of H0: all slope coefficients = 0.
# Structural tests (OLS cross-check, edge cases) + Stata fixture tests.

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

sim_noconst_path <- file.path(fixture_dir, "sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

sim_cluster_path <- file.path(fixture_dir, "sim_cluster_data.csv")
if (file.exists(sim_cluster_path)) {
  sim_cluster <- read.csv(sim_cluster_path)
}


# ============================================================================
# Structural tests (no fixtures)
# ============================================================================

test_that("OLS model F matches summary(lm(...))$fstatistic", {
  fit_iv <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  lm_f <- summary(lm_fit)$fstatistic

  expect_equal(fit_iv$model_f, unname(lm_f["value"]), tolerance = 1e-10)
  expect_equal(fit_iv$model_f_df1, as.integer(lm_f["numdf"]))
  expect_equal(fit_iv$model_f_df2, as.integer(lm_f["dendf"]))
  expect_equal(fit_iv$model_f_p,
               unname(pf(lm_f["value"], lm_f["numdf"], lm_f["dendf"],
                         lower.tail = FALSE)),
               tolerance = 1e-10)
})

test_that("intercept-only model returns model_f = NA", {
  fit <- ivreg2(mpg ~ 1, data = mtcars)
  expect_true(is.na(fit$model_f))
  expect_true(is.na(fit$model_f_p))
  expect_true(is.na(fit$model_f_df1))
  expect_true(is.na(fit$model_f_df2))
})

test_that("df1 = K - 1 with intercept, df1 = K without", {
  # With intercept: K = 3, df1 = 2
  fit1 <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_identical(fit1$model_f_df1, 2L)

  # Without intercept: K = 2, df1 = 2
  fit2 <- ivreg2(mpg ~ 0 + wt + hp, data = mtcars)
  expect_identical(fit2$model_f_df1, 2L)
})

test_that("df2 = N - K for non-cluster models", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_identical(fit$model_f_df2, as.integer(nrow(mtcars) - 3L))
})

test_that("model_f is populated for OLS models", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_true(is.numeric(fit$model_f))
  expect_false(is.na(fit$model_f))
})

test_that("model_f is populated for IV models", {
  skip_if(!file.exists(card_path), "card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  expect_true(is.numeric(fit$model_f))
  expect_false(is.na(fit$model_f))
})

test_that("model F-stat is invariant to small (OLS IID)", {
  fit1 <- ivreg2(mpg ~ wt + hp, data = mtcars, small = FALSE)
  fit2 <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  expect_equal(fit1$model_f, fit2$model_f)
})

test_that("model F-stat is invariant to small (HC)", {
  fit1 <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC0")
  fit2 <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC1")
  expect_equal(fit1$model_f, fit2$model_f)
})


# ============================================================================
# Helper: run model F fixture comparison for one spec/VCE combination
# ============================================================================

check_model_f <- function(fit, fixture_path, label) {
  fixture <- read_diagnostics(fixture_path)

  test_that(paste("model F-stat matches Stata", label), {
    expect_equal(fit$model_f, fixture$F_stat,
                 tolerance = stata_tol$stat)
  })

  test_that(paste("model F p-value matches Stata", label), {
    expect_equal(fit$model_f_p, fixture$F_p,
                 tolerance = stata_tol$pval)
  })

  test_that(paste("model F df1 matches Stata", label), {
    expect_identical(fit$model_f_df1, as.integer(fixture$F_df1))
  })

  test_that(paste("model F df2 matches Stata", label), {
    expect_identical(fit$model_f_df2, as.integer(fixture$F_df2))
  })
}


# ============================================================================
# card_just_id: lwage ~ exper+expersq+black+south | educ | nearc4
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1"),
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
    check_model_f(fit, fixture_file, label)
  }
}


# ============================================================================
# card_overid: lwage ~ exper+expersq+black+south | educ | nearc2+nearc4
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1"),
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
    check_model_f(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_multi_endo: y ~ x1+x2 | endo1+endo2 | z1+z2+z3+z4
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1"),
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
    check_model_f(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_no_constant: y ~ 0+x1 | endo1 | z1+z2
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1"),
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
    check_model_f(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_cluster: y ~ x1 | endo1 | z1+z2, clusters=~cluster_id
# Includes IID, HC, and CL variants
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, clusters = NULL,          suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  clusters = NULL,          suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, clusters = NULL,          suffix = "hc1"),
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
    check_model_f(fit, fixture_file, label)
  }
}


# ============================================================================
# card_just_id_weighted: lwage ~ exper+expersq+black+south | educ | nearc4
# with weights=weight; cluster on smsa66 (M=2)
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, clusters = NULL,       suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  clusters = NULL,       suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, clusters = NULL,       suffix = "hc1"),
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
    check_model_f(fit, fixture_file, label)
  }
}


# ============================================================================
# Direct tests for .syminv_sweep()
# ============================================================================
# The sweep-based generalized inverse is only used when the VCV submatrix is
# rank-deficient (M < K clusters). These tests verify the algorithm directly
# rather than through the model F end-to-end path.

test_that(".syminv_sweep on full-rank matrix matches solve()", {
  set.seed(42)
  A <- crossprod(matrix(rnorm(9), 3, 3))  # 3x3 positive definite
  G <- ivreg2r:::.syminv_sweep(A)
  expect_equal(G, solve(A), tolerance = 1e-10)
})

test_that(".syminv_sweep on rank-1 matrix picks largest diagonal as pivot", {
  # V = v v' where v = (1, 2, 3). Diag = (1, 4, 9).
  # Largest diagonal is element 3 (value 9).
  # Sweep should pivot on element 3 only, giving G[3,3] = 1/9, rest = 0.
  v <- c(1, 2, 3)
  V <- outer(v, v)
  G <- ivreg2r:::.syminv_sweep(V)

  expected <- matrix(0, 3, 3)
  expected[3, 3] <- 1 / 9
  expect_equal(G, expected, tolerance = 1e-12)
})

test_that(".syminv_sweep on rank-2 matrix zeros out the right dimension", {
  # Build a rank-2 matrix from two outer products
  v1 <- c(3, 0, 0, 0)
  v2 <- c(0, 0, 0, 2)
  V <- outer(v1, v1) + outer(v2, v2)
  # Diag = (9, 0, 0, 4). Pivots should be elements 1 (9) and 4 (4).
  G <- ivreg2r:::.syminv_sweep(V)

  expected <- matrix(0, 4, 4)
  expected[1, 1] <- 1 / 9
  expected[4, 4] <- 1 / 4
  expect_equal(G, expected, tolerance = 1e-12)
})

test_that(".syminv_sweep returns NULL for zero matrix", {
  V <- matrix(0, 3, 3)
  expect_null(ivreg2r:::.syminv_sweep(V))
})

test_that(".syminv_sweep produces symmetric output", {
  # Rank-2 case with off-diagonal structure
  set.seed(123)
  U <- matrix(rnorm(12), 4, 3)  # 4x3, so U %*% t(U) has rank 3
  V <- tcrossprod(U)             # 4x4, rank 3
  G <- ivreg2r:::.syminv_sweep(V)
  expect_equal(G, t(G), tolerance = 1e-14)
})

test_that(".syminv_sweep satisfies A G A = A on pivoted subspace", {
  # For a valid g-inverse, A G A should equal A on the pivoted subspace.
  v <- c(1, 2, 3)
  V <- outer(v, v)  # rank 1
  G <- ivreg2r:::.syminv_sweep(V)
  AGA <- V %*% G %*% V

  # Only element (3,3) was pivoted, so AGA[3,3] should equal V[3,3]
  # and AGA should project V onto the pivoted subspace
  expect_equal(AGA[3, 3], V[3, 3], tolerance = 1e-12)
})

test_that(".syminv_sweep chi2 matches Stata for M=2 cluster case", {
  # This is the actual case that motivated the implementation:
  # Card data, weighted IV with M=2 clusters on smsa66.
  # Stata's `test` drops 4 of 5 constraints and reports chi2 = 1.546.
  # The surviving constraint is "black" (largest diagonal in VCV submatrix).
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight,
                clusters = ~smsa66, small = TRUE)

  # Verify the model F was computed (not NA) despite rank-deficient VCV
  expect_false(is.na(fit$model_f))
  expect_equal(fit$model_f_df2, 1L)  # M - 1 = 2 - 1

  # Cross-check: the chi2 should be beta_black^2 / V[black,black]
  # divided by df1=5 for the F-stat
  slope_idx <- setdiff(seq_along(fit$coefficients),
                       match("(Intercept)", names(fit$coefficients)))
  V_s <- fit$vcov[slope_idx, slope_idx]
  black_pos <- match("black", names(fit$coefficients[slope_idx]))
  chi2_manual <- unname(fit$coefficients[slope_idx][black_pos]^2 / V_s[black_pos, black_pos])
  expect_equal(fit$model_f, chi2_manual / length(slope_idx), tolerance = 1e-10)
})

test_that(".compute_model_f returns NA for indefinite VCV", {
  result <- ivreg2r:::.compute_model_f(
    coefficients = c("(Intercept)" = 1, x = 2),
    vcov = matrix(c(1, 0, 0, -1), 2, 2),  # indefinite
    N = 100L, K = 2L,
    has_intercept = TRUE,
    vcov_type = "iid",
    small = TRUE
  )
  expect_true(is.na(result$model_f))
})
