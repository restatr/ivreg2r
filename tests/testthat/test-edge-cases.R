# ============================================================================
# Tests: Edge Cases (Ticket L2)
# ============================================================================
#
# 1. Intercept-only exogenous model (y ~ 1 | endo1 | z1 + z2)
# 2. Near-singular KP matrix defensive guards
#
# Verifies against Stata ivreg2 fixtures and pure R unit tests.

# --- Helpers ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)

read_diagnostics <- function(path) {
  read.csv(path, check.names = FALSE)
}

read_firststage <- function(path) {
  read.csv(path, check.names = FALSE)
}

get_fs_value <- function(fixture, stat, endo_name) {
  as.numeric(fixture[fixture$statistic == stat, endo_name])
}

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


# --- Load dataset ---
intonly_path <- file.path(fixture_dir, "sim_intercept_only_data.csv")
if (file.exists(intonly_path)) {
  sim_intonly <- read.csv(intonly_path)
}


# ============================================================================
# Intercept-only: Coefficients
# ============================================================================

test_that("coefficients match Stata sim_intercept_only iid fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  coef_path <- file.path(fixture_dir, "sim_intercept_only_coef_iid.csv")
  skip_if(!file.exists(coef_path), "coefficient fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly)
  fixture <- read.csv(coef_path)

  for (i in seq_len(nrow(fixture))) {
    term <- fixture$term[i]
    r_name <- if (term == "_cons") "(Intercept)" else term
    expect_true(r_name %in% names(coef(fit)),
                info = paste("Missing coefficient:", r_name))
    expect_equal(
      unname(coef(fit)[r_name]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coefficient mismatch:", r_name)
    )
  }
})

test_that("coefficients match Stata sim_intercept_only iid_small fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  coef_path <- file.path(fixture_dir, "sim_intercept_only_coef_iid_small.csv")
  skip_if(!file.exists(coef_path), "coefficient fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly, small = TRUE)
  fixture <- read.csv(coef_path)

  for (i in seq_len(nrow(fixture))) {
    term <- fixture$term[i]
    r_name <- if (term == "_cons") "(Intercept)" else term
    expect_equal(
      unname(coef(fit)[r_name]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coefficient mismatch:", r_name)
    )
  }
})


# ============================================================================
# Intercept-only: VCV
# ============================================================================

test_that("IID VCV matches Stata sim_intercept_only fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  vcov_path <- file.path(fixture_dir, "sim_intercept_only_vcov_iid.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

test_that("HC0 VCV matches Stata sim_intercept_only robust fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  vcov_path <- file.path(fixture_dir, "sim_intercept_only_vcov_hc1.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly, vcov = "HC0")
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

test_that("HC1 VCV matches Stata sim_intercept_only robust small fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  vcov_path <- file.path(fixture_dir, "sim_intercept_only_vcov_hc1_small.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly, vcov = "HC1")
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})


# ============================================================================
# Intercept-only: Overidentification tests (overid_df = 1)
# ============================================================================

test_that("Sargan stat matches Stata sim_intercept_only iid fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  diag_path <- file.path(fixture_dir, "sim_intercept_only_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Sargan")
  expect_equal(fit$diagnostics$overid$stat, fixture$sargan,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$sarganp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$sargandf))
})

test_that("Hansen J stat matches Stata sim_intercept_only hc1 fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  diag_path <- file.path(fixture_dir, "sim_intercept_only_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$overid$p, fixture$jp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$overid$df, as.integer(fixture$jdf))
})


# ============================================================================
# Intercept-only: Identification tests
# ============================================================================

test_that("Anderson LM matches Stata sim_intercept_only iid fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  diag_path <- file.path(fixture_dir, "sim_intercept_only_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("Cragg-Donald F matches Stata sim_intercept_only iid fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  diag_path <- file.path(fixture_dir, "sim_intercept_only_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})

test_that("KP rk LM matches Stata sim_intercept_only hc1 fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  diag_path <- file.path(fixture_dir, "sim_intercept_only_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("KP rk Wald F matches Stata sim_intercept_only hc1 fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  diag_path <- file.path(fixture_dir, "sim_intercept_only_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id_robust$stat, fixture$widstat,
               tolerance = stata_tol$stat)
})


# ============================================================================
# Intercept-only: Anderson-Rubin test
# ============================================================================

check_anderson_rubin <- function(fit, fixture_path, label) {
  fixture <- read_diagnostics(fixture_path)
  ar <- fit$diagnostics$anderson_rubin

  test_that(paste("AR F-stat matches Stata", label), {
    expect_equal(ar$f_stat, fixture$arf,
                 tolerance = stata_tol$stat)
  })
  test_that(paste("AR F p-value matches Stata", label), {
    expect_equal(ar$f_p, fixture$arfp,
                 tolerance = stata_tol$pval)
  })
  test_that(paste("AR chi2 matches Stata", label), {
    expect_equal(ar$chi2_stat, fixture$archi2,
                 tolerance = stata_tol$stat)
  })
  test_that(paste("AR chi2 p-value matches Stata", label), {
    expect_equal(ar$chi2_p, fixture$archi2p,
                 tolerance = stata_tol$pval)
  })
  test_that(paste("AR F df1 matches Stata", label), {
    expect_identical(ar$f_df1, as.integer(fixture$ardf))
  })
  test_that(paste("AR F df2 matches Stata", label), {
    expect_identical(ar$f_df2, as.integer(fixture$ardf_r))
  })
  test_that(paste("AR chi2 df matches Stata", label), {
    expect_identical(ar$chi2_df, as.integer(fixture$ardf))
  })
}

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC1",  small = FALSE, suffix = "hc1"),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_intercept_only_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_intercept_only", vce_combo$suffix)

  if (file.exists(intonly_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ 1 | endo1 | z1 + z2,
                  data = sim_intonly, vcov = vce_combo$vcov,
                  small = vce_combo$small)
    check_anderson_rubin(fit, fixture_file, label)
  }
}


# ============================================================================
# Intercept-only: Stock-Wright S statistic
# ============================================================================

check_stock_wright <- function(fit, fixture_path, label) {
  fixture <- read_diagnostics(fixture_path)
  sw <- fit$diagnostics$stock_wright
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

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC1",  small = FALSE, suffix = "hc1"),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- file.path(
    fixture_dir,
    paste0("sim_intercept_only_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_intercept_only", vce_combo$suffix)

  if (file.exists(intonly_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ 1 | endo1 | z1 + z2,
                  data = sim_intonly, vcov = vce_combo$vcov,
                  small = vce_combo$small)
    check_stock_wright(fit, fixture_file, label)
  }
}


# ============================================================================
# Intercept-only: First-stage diagnostics
# ============================================================================

test_that("first-stage F matches Stata sim_intercept_only iid fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  fs_path <- file.path(fixture_dir, "sim_intercept_only_firststage_iid.csv")
  skip_if(!file.exists(fs_path), "firststage fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly)
  fixture <- read_firststage(fs_path)
  fs <- fit$first_stage$endo1

  expect_equal(fs$f_stat, get_fs_value(fixture, "F", "endo1"),
               tolerance = stata_tol$stat)
  expect_equal(fs$partial_r2, get_fs_value(fixture, "pr2", "endo1"),
               tolerance = stata_tol$coef)
  expect_equal(fs$shea_partial_r2, get_fs_value(fixture, "sheapr2", "endo1"),
               tolerance = stata_tol$coef)
  expect_equal(fs$rmse, get_fs_value(fixture, "rmse", "endo1"),
               tolerance = stata_tol$coef)
})

test_that("first-stage F matches Stata sim_intercept_only hc1 fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  fs_path <- file.path(fixture_dir, "sim_intercept_only_firststage_hc1.csv")
  skip_if(!file.exists(fs_path), "firststage fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly, vcov = "HC1")
  fixture <- read_firststage(fs_path)
  fs <- fit$first_stage$endo1

  expect_equal(fs$f_stat, get_fs_value(fixture, "F", "endo1"),
               tolerance = stata_tol$stat)
})


# ============================================================================
# Intercept-only: Model F
# ============================================================================

test_that("model F matches Stata sim_intercept_only iid fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  diag_path <- file.path(fixture_dir, "sim_intercept_only_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$model_f, fixture$F_stat, tolerance = stata_tol$stat)
  expect_equal(fit$model_f_p, fixture$F_p, tolerance = stata_tol$pval)
  expect_identical(fit$model_f_df1, as.integer(fixture$F_df1))
  expect_identical(fit$model_f_df2, as.integer(fixture$F_df2))
})

test_that("model F matches Stata sim_intercept_only hc1 fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  diag_path <- file.path(fixture_dir, "sim_intercept_only_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$model_f, fixture$F_stat, tolerance = stata_tol$stat)
  expect_equal(fit$model_f_p, fixture$F_p, tolerance = stata_tol$pval)
})


# ============================================================================
# Intercept-only: Summary statistics (R², adj R², RMSE)
# ============================================================================

test_that("R² matches Stata sim_intercept_only iid fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  diag_path <- file.path(fixture_dir, "sim_intercept_only_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$r.squared, fixture$r2, tolerance = stata_tol$coef)
  expect_equal(fit$adj.r.squared, fixture$r2_a, tolerance = stata_tol$coef)
  expect_equal(fit$sigma, fixture$rmse, tolerance = stata_tol$coef)
  expect_equal(fit$rss, fixture$rss, tolerance = stata_tol$coef)
})


# ============================================================================
# Intercept-only: Endogeneity test
# ============================================================================

test_that("endogeneity test matches Stata sim_intercept_only iid fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  diag_path <- file.path(fixture_dir, "sim_intercept_only_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$endogeneity$stat, fixture$estat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$endogeneity$p, fixture$estatp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$endogeneity$df, as.integer(fixture$estatdf))
})

test_that("endogeneity test matches Stata sim_intercept_only hc1 fixture", {
  skip_if(!file.exists(intonly_path), "sim_intercept_only data not found")
  diag_path <- file.path(fixture_dir, "sim_intercept_only_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 1 | endo1 | z1 + z2, data = sim_intonly, vcov = "HC1")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$endogeneity$stat, fixture$estat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$endogeneity$p, fixture$estatp,
               tolerance = stata_tol$pval)
})


# ============================================================================
# Near-singular KP matrix unit tests (no Stata fixtures)
# ============================================================================

test_that(".sym_sqrt clamps negative eigenvalues with warning", {
  # Construct a 2x2 matrix with a tiny negative eigenvalue
  # Start from eigendecomposition: eigenvalues (1, -1e-12)
  Q <- matrix(c(1, 1, 1, -1), nrow = 2) / sqrt(2)  # orthogonal
  A <- Q %*% diag(c(1, -1e-12)) %*% t(Q)

  expect_warning(
    result <- ivreg2r:::.sym_sqrt(A),
    "Negative eigenvalues clamped to 0"
  )
  # Result should be a valid matrix (no NaN)
  expect_false(any(is.nan(result)))
  expect_true(is.matrix(result))
})

test_that(".canonical_correlations returns NULL for zero instruments", {
  # Z1_perp = zero matrix → Qzz = 0 → Cholesky fails
  N <- 50
  Z1_perp <- matrix(0, nrow = N, ncol = 2)
  X1_perp <- matrix(rnorm(N), ncol = 1)

  expect_warning(
    result <- ivreg2r:::.canonical_correlations(X1_perp, Z1_perp, N,
                                                 K1 = 1L, L1 = 2L,
                                                 weights = NULL),
    "Cholesky factorization failed"
  )
  expect_null(result)
})

test_that(".canonical_correlations handles near-identical X1_perp and Z1_perp", {
  # When X1_perp == Z1_perp, mathematical cc = 1 but float precision
  # may produce cc slightly < 1. The function should either return NULL
  # with a warning (if cc >= 1) or a valid result with cc close to 1.
  set.seed(43)
  N <- 100
  z <- rnorm(N)
  Z1_perp <- matrix(z, ncol = 1)
  X1_perp <- matrix(z, ncol = 1)

  result <- tryCatch(
    withCallingHandlers(
      ivreg2r:::.canonical_correlations(X1_perp, Z1_perp, N,
                                         K1 = 1L, L1 = 1L,
                                         weights = NULL),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) "error"
  )

  if (is.null(result)) {
    # cc >= 1 triggered → returned NULL (expected edge case)
    succeed()
  } else {
    # cc < 1 due to float precision → valid result with cc ≈ 1
    expect_true(result$cc[1] > 0.999)
  }
})

test_that(".kp_rk_stat uses pseudoinverse for rank-deficient vlab", {
  # Construct inputs where vlab will be rank-deficient:
  # Use a synthetic setup with near-zero score covariance
  set.seed(44)
  N <- 100
  K1 <- 1L
  L1 <- 2L

  # Create valid cc_result components manually
  Z1_perp <- matrix(rnorm(N * L1), nrow = N, ncol = L1)
  X1_perp <- matrix(rnorm(N * K1), nrow = N, ncol = K1)
  Qzz <- crossprod(Z1_perp) / N
  Qzy <- crossprod(Z1_perp, X1_perp) / N
  Qyy <- crossprod(X1_perp) / N
  R_zz <- chol(Qzz)
  R_yy <- chol(Qyy)
  irQzz <- backsolve(R_zz, diag(L1))
  irQyy <- backsolve(R_yy, diag(K1))
  pihat <- chol2inv(R_zz) %*% Qzy
  theta <- R_zz %*% pihat %*% irQyy
  sv <- svd(theta, nu = L1, nv = K1)

  cc_result <- list(
    theta = theta, U = sv$u, cc = sv$d[1], V = sv$v,
    eval = sv$d[1]^2, pihat = pihat, irQyy = irQyy, irQzz = irQzz
  )

  # Create a rank-deficient shat0 (K1*L1 x K1*L1 = 2x2)
  v <- c(1, 0)
  shat0 <- outer(v, v)  # rank 1 (rank-deficient for 2x2)

  # Should produce a warning about rank deficiency and still compute
  expect_warning(
    result <- ivreg2r:::.kp_rk_stat(cc_result, shat0, N, K1, L1),
    "rank-deficient"
  )
  expect_true(is.finite(result$chi2))
  expect_true(is.integer(result$df))
})

test_that("ivreg2() with near-collinear instrument returns finite diagnostics", {
  set.seed(45)
  n <- 200

  # Exogenous regressor
  x <- rnorm(n)

  # Two excluded instruments, second nearly collinear with x
  z1 <- rnorm(n)
  z2 <- x + rnorm(n, sd = 0.01)  # nearly collinear with x

  # Structural errors and endogenous variable
  u <- rnorm(n)
  endo <- 0.5 * z1 + 0.3 * z2 + 0.6 * u

  # Outcome
  y <- 1 + x + 1.5 * endo + u
  d <- data.frame(y = y, x = x, endo = endo, z1 = z1, z2 = z2)

  # Should run without error; diagnostics should be finite
  fit <- ivreg2(y ~ x | endo | z1 + z2, data = d, vcov = "HC1")

  expect_true(is.finite(fit$diagnostics$underid$stat))
  expect_true(is.finite(fit$diagnostics$weak_id$stat))
  expect_true(is.finite(fit$diagnostics$weak_id_robust$stat))
  expect_true(is.finite(fit$diagnostics$overid$stat))
})
