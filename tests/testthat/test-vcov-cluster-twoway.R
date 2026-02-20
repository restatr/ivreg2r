# ============================================================================
# Tests: Two-way cluster-robust VCV (Ticket I1)
# ============================================================================
#
# Cameron-Gelbach-Miller (2006) two-way clustering:
# V_twoway = V_c1 + V_c2 - V_intersection
# Effective cluster count M = min(M1, M2) per Stata convention.

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

# --- Load sim_twoway data ---
sim_path <- file.path(fixture_dir, "sim_twoway_data.csv")
if (file.exists(sim_path)) {
  sim_twoway <- read.csv(sim_path)
}


# ============================================================================
# Section 1: Stata parity — coefficients, VCV, SEs (cl2)
# ============================================================================

test_that("2SLS two-way cluster VCV matches Stata cl2 fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  vcov_path <- file.path(fixture_dir, "sim_twoway_vcov_cl2.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

test_that("2SLS two-way cluster SEs match Stata cl2 fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  coef_path <- file.path(fixture_dir, "sim_twoway_coef_cl2.csv")
  skip_if(!file.exists(coef_path), "Coef fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id)
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

test_that("2SLS two-way cluster coefficients match Stata cl2 fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  coef_path <- file.path(fixture_dir, "sim_twoway_coef_cl2.csv")
  skip_if(!file.exists(coef_path), "Coef fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id)
  fixture <- read.csv(coef_path)
  fixture$r_name <- ifelse(fixture$term == "_cons", "(Intercept)", fixture$term)

  for (i in seq_len(nrow(fixture))) {
    nm <- fixture$r_name[i]
    expect_equal(
      unname(coef(fit)[nm]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", nm)
    )
  }
})


# ============================================================================
# Section 2: Stata parity — VCV and SEs (cl2_small)
# ============================================================================

test_that("2SLS two-way cluster VCV matches Stata cl2_small fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  vcov_path <- file.path(fixture_dir, "sim_twoway_vcov_cl2_small.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id,
                small = TRUE)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

test_that("2SLS two-way cluster SEs match Stata cl2_small fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  coef_path <- file.path(fixture_dir, "sim_twoway_coef_cl2_small.csv")
  skip_if(!file.exists(coef_path), "Coef fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id,
                small = TRUE)
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
# Section 3: Stata parity — diagnostics (cl2 and cl2_small)
# ============================================================================

test_that("two-way cluster diagnostics match Stata cl2 fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  diag_path <- file.path(fixture_dir, "sim_twoway_diagnostics_cl2.csv")
  skip_if(!file.exists(diag_path), "Diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id)
  stata <- read.csv(diag_path)

  # Hansen J
  expect_equal(fit$diagnostics$overid$stat, stata$j,
               tolerance = stata_tol$stat, info = "Hansen J stat")
  expect_equal(fit$diagnostics$overid$p, stata$jp,
               tolerance = stata_tol$pval, info = "Hansen J p")

  # KP rk LM (underidentification)
  expect_equal(fit$diagnostics$underid$stat, stata$idstat,
               tolerance = stata_tol$stat, info = "KP rk LM stat")
  expect_equal(fit$diagnostics$underid$p, stata$idp,
               tolerance = stata_tol$pval, info = "KP rk LM p")

  # KP rk Wald F (weak identification, robust)
  expect_equal(fit$diagnostics$weak_id_robust$stat, stata$widstat,
               tolerance = stata_tol$stat, info = "KP rk Wald F")

  # Anderson-Rubin
  expect_equal(fit$diagnostics$anderson_rubin$f_stat, stata$arf,
               tolerance = stata_tol$stat, info = "AR F stat")
  expect_equal(fit$diagnostics$anderson_rubin$f_p, stata$arfp,
               tolerance = stata_tol$pval, info = "AR F p")
  expect_equal(fit$diagnostics$anderson_rubin$chi2_stat, stata$archi2,
               tolerance = stata_tol$stat, info = "AR chi2 stat")
  expect_equal(fit$diagnostics$anderson_rubin$chi2_p, stata$archi2p,
               tolerance = stata_tol$pval, info = "AR chi2 p")

  # Stock-Wright S
  expect_equal(fit$diagnostics$stock_wright$stat, stata$sstat,
               tolerance = stata_tol$stat, info = "Stock-Wright S stat")
  expect_equal(fit$diagnostics$stock_wright$p, stata$sstatp,
               tolerance = stata_tol$pval, info = "Stock-Wright S p")

  # Endogeneity test
  expect_equal(fit$diagnostics$endogeneity$stat, stata$estat,
               tolerance = stata_tol$stat, info = "Endogeneity stat")
  expect_equal(fit$diagnostics$endogeneity$p, stata$estatp,
               tolerance = stata_tol$pval, info = "Endogeneity p")

  # Model F
  expect_equal(fit$model_f, stata$F_stat,
               tolerance = stata_tol$stat, info = "Model F stat")
  expect_equal(fit$model_f_p, stata$F_p,
               tolerance = stata_tol$pval, info = "Model F p")
})

test_that("two-way cluster diagnostics match Stata cl2_small fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  diag_path <- file.path(fixture_dir, "sim_twoway_diagnostics_cl2_small.csv")
  skip_if(!file.exists(diag_path), "Diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id,
                small = TRUE)
  stata <- read.csv(diag_path)

  # Model F (small)
  expect_equal(fit$model_f, stata$F_stat,
               tolerance = stata_tol$stat, info = "Model F stat (small)")
  expect_equal(fit$model_f_p, stata$F_p,
               tolerance = stata_tol$pval, info = "Model F p (small)")

  # N_clust
  expect_equal(as.integer(stata$N_clust), fit$n_clusters,
               info = "N_clust")
})


# ============================================================================
# Section 4: Stata parity — first-stage F
# ============================================================================

test_that("two-way cluster first-stage F matches Stata cl2 fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  fs_path <- file.path(fixture_dir, "sim_twoway_firststage_cl2.csv")
  skip_if(!file.exists(fs_path), "First-stage fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id)
  stata_fs <- read.csv(fs_path)

  # First-stage F for endo1
  stata_f <- stata_fs$endo1[stata_fs$statistic == "F"]
  expect_equal(fit$first_stage$endo1$f_stat, stata_f,
               tolerance = stata_tol$stat, info = "First-stage F")
})


# ============================================================================
# Section 5: Stata parity — OLS two-way cluster
# ============================================================================

test_that("OLS two-way cluster VCV matches Stata cl2_ols fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  vcov_path <- file.path(fixture_dir, "sim_twoway_vcov_cl2_ols.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ x1 + endo1,
                data = sim_twoway, clusters = ~ firm_id + year_id)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

test_that("OLS two-way cluster VCV matches Stata cl2_ols_small fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  vcov_path <- file.path(fixture_dir, "sim_twoway_vcov_cl2_ols_small.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ x1 + endo1,
                data = sim_twoway, clusters = ~ firm_id + year_id,
                small = TRUE)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})


# ============================================================================
# Section 6: Stata parity — weighted two-way cluster
# ============================================================================

test_that("weighted 2SLS two-way cluster VCV matches Stata cl2_wt fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  vcov_path <- file.path(fixture_dir, "sim_twoway_vcov_cl2_wt.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id,
                weights = wt)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

test_that("weighted 2SLS two-way cluster VCV matches Stata cl2_wt_small fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  vcov_path <- file.path(fixture_dir, "sim_twoway_vcov_cl2_wt_small.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id,
                weights = wt, small = TRUE)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

test_that("weighted two-way cluster diagnostics match Stata cl2_wt fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  diag_path <- file.path(fixture_dir, "sim_twoway_diagnostics_cl2_wt.csv")
  skip_if(!file.exists(diag_path), "Diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id,
                weights = wt)
  stata <- read.csv(diag_path)

  # Note: Stata warns "covariance matrix of moment conditions not of full
  # rank" for the weighted two-way cluster case. It reports J = NA and flags
  # all diagnostics as unreliable. We only compare diagnostics where Stata
  # reports a value AND does not issue a rank-deficiency warning.
  # VCV and model F are still valid since they use the main sandwich VCV.
  expect_equal(fit$model_f, stata$F_stat,
               tolerance = stata_tol$stat, info = "Model F (weighted)")
})


# ============================================================================
# Section 7: Stata parity — dofminus/sdofminus two-way cluster
# ============================================================================

test_that("two-way cluster with dofminus VCV matches Stata cl2_dof fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  vcov_path <- file.path(fixture_dir, "sim_twoway_vcov_cl2_dof.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id,
                dofminus = 1L, sdofminus = 1L)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

test_that("two-way cluster with dofminus VCV matches Stata cl2_dof_small fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  vcov_path <- file.path(fixture_dir, "sim_twoway_vcov_cl2_dof_small.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id,
                dofminus = 1L, sdofminus = 1L, small = TRUE)
  expect_vcov_equal(fit$vcov, read_vcov_fixture(vcov_path))
})

test_that("two-way cluster dofminus diagnostics match Stata cl2_dof fixture", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")
  diag_path <- file.path(fixture_dir, "sim_twoway_diagnostics_cl2_dof.csv")
  skip_if(!file.exists(diag_path), "Diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id,
                dofminus = 1L, sdofminus = 1L)
  stata <- read.csv(diag_path)

  expect_equal(fit$diagnostics$overid$stat, stata$j,
               tolerance = stata_tol$stat, info = "Hansen J (dofminus)")
  expect_equal(fit$diagnostics$underid$stat, stata$idstat,
               tolerance = stata_tol$stat, info = "KP rk LM (dofminus)")
  expect_equal(fit$model_f, stata$F_stat,
               tolerance = stata_tol$stat, info = "Model F (dofminus)")
})


# ============================================================================
# Section 8: Metadata
# ============================================================================

test_that("n_clusters = min(M1, M2) for two-way clustering", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id)

  M1 <- length(unique(sim_twoway$firm_id))
  M2 <- length(unique(sim_twoway$year_id))
  expect_identical(fit$n_clusters, as.integer(min(M1, M2)))
  expect_identical(fit$n_clusters1, as.integer(M1))
  expect_identical(fit$n_clusters2, as.integer(M2))
})

test_that("df.residual = M - 1 for two-way clustering", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id)

  M1 <- length(unique(sim_twoway$firm_id))
  M2 <- length(unique(sim_twoway$year_id))
  M <- min(M1, M2)
  expect_identical(fit$df.residual, as.integer(M - 1L))
})

test_that("vcov_type is 'CL' when two-way clustered", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id)
  expect_identical(fit$vcov_type, "CL")
})

test_that("cluster_var stores both variable names for two-way", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id)
  expect_identical(fit$cluster_var, c("firm_id", "year_id"))
})

test_that("n_clusters1 and n_clusters2 are NULL for one-way clustering", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id)
  expect_null(fit$n_clusters1)
  expect_null(fit$n_clusters2)
})

test_that("glance includes n_clusters1 and n_clusters2", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id)
  gl <- glance(fit)
  expect_true("n_clusters1" %in% names(gl))
  expect_true("n_clusters2" %in% names(gl))
  expect_equal(gl$n_clusters1, fit$n_clusters1)
  expect_equal(gl$n_clusters2, fit$n_clusters2)
})


# ============================================================================
# Section 9: Properties
# ============================================================================

test_that("two-way cluster VCV is symmetric", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")

  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                data = sim_twoway, clusters = ~ firm_id + year_id)
  expect_equal(fit$vcov, t(fit$vcov),
               tolerance = .Machine$double.eps^0.5)
})

test_that("coefficients are identical across iid / HC1 / one-way CL / two-way CL", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")

  fit_iid <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                    data = sim_twoway, vcov = "iid")
  fit_hc  <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                    data = sim_twoway, vcov = "HC1")
  fit_cl1 <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                    data = sim_twoway, clusters = ~ firm_id)
  fit_cl2 <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                    data = sim_twoway, clusters = ~ firm_id + year_id)
  expect_equal(coef(fit_hc), coef(fit_iid),
               tolerance = .Machine$double.eps^0.5)
  expect_equal(coef(fit_cl1), coef(fit_iid),
               tolerance = .Machine$double.eps^0.5)
  expect_equal(coef(fit_cl2), coef(fit_iid),
               tolerance = .Machine$double.eps^0.5)
})

test_that("two-way cluster VCV differs from one-way cluster VCV", {
  skip_if(!file.exists(sim_path), "sim_twoway data not found")

  fit_cl1 <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                    data = sim_twoway, clusters = ~ firm_id)
  fit_cl2 <- ivreg2(y ~ x1 | endo1 | z1 + z2,
                    data = sim_twoway, clusters = ~ firm_id + year_id)
  expect_false(isTRUE(all.equal(fit_cl1$vcov, fit_cl2$vcov)))
})


# ============================================================================
# Section 10: Input validation
# ============================================================================

test_that("3+ cluster variables gives error", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, clusters = ~ cyl + gear + carb),
    "one or two variables"
  )
})

test_that("missing second cluster variable gives error", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, clusters = ~ cyl + nonexistent),
    "not found in data"
  )
})

test_that("NA in second cluster variable gives error", {
  df <- mtcars
  df$clust2 <- c(NA, rep(1:3, length.out = nrow(df) - 1))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = df, clusters = ~ cyl + clust2),
    "contains NA"
  )
})

test_that("single-level second cluster variable gives error", {
  df <- mtcars
  df$const <- 1
  expect_error(
    ivreg2(mpg ~ wt + hp, data = df, clusters = ~ cyl + const),
    "At least 2 clusters"
  )
})

test_that("interaction operator in clusters formula gives error", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, clusters = ~ cyl:gear),
    "not `:` or `\\*`"
  )
})

test_that("star operator in clusters formula gives error", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, clusters = ~ cyl * gear),
    "not `:` or `\\*`"
  )
})
