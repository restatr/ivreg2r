# ============================================================================
# Tests: Two-step efficient GMM estimation (Ticket N1)
# ============================================================================

# --- Helpers ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)
card_path <- file.path(fixture_dir, "card_data.csv")

if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

read_vcov_fixture <- function(path) {
  fixture <- read.csv(path)
  # The VCV fixture has columns v1, v2, ..., vK — no term column
  V <- as.matrix(fixture)
  V
}

# ============================================================================
# 1. IID GMM2S — coefficients differ from 2SLS algebraically when IID omega
#    is used, because the IID omega (sigma^2 * Z'Z/N) gives the same
#    weighting as 2SLS. So GMM2S + IID should match 2SLS.
# ============================================================================

test_that("GMM2S IID matches 2SLS coefficients (algebraic equivalence)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card)
  fit_gmm <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                    data = card, method = "gmm2s")

  # Coefficients should be nearly identical (same algebra, different code path)
  expect_equal(coef(fit_gmm), coef(fit_2sls), tolerance = 1e-10)
})

test_that("GMM2S IID coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_coef_iid.csv")
  skip_if(!file.exists(fixture_path), "GMM2S IID fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s")
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("GMM2S IID diagnostics match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_diagnostics_iid.csv")
  skip_if(!file.exists(fixture_path), "GMM2S IID diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s")
  diag <- read.csv(fixture_path)

  # Overid (Sargan for IID)
  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = stata_tol$stat, info = "Sargan stat")
  expect_equal(fit$diagnostics$overid$p, diag$overid_p,
               tolerance = stata_tol$pval, info = "Sargan p")
  expect_equal(fit$diagnostics$overid$df, diag$overid_df, info = "Sargan df")
  expect_equal(fit$diagnostics$overid$test_name, "Sargan")

  # Underid
  expect_equal(fit$diagnostics$underid$stat, diag$underid_stat,
               tolerance = stata_tol$stat, info = "underid stat")

  # Weak id
  expect_equal(fit$diagnostics$weak_id$stat, diag$weak_id_cd_f,
               tolerance = stata_tol$stat, info = "CD F")

  # Summary stats
  expect_equal(fit$sigma, diag$sigma, tolerance = stata_tol$coef,
               info = "sigma")
  expect_equal(fit$rss, diag$rss, tolerance = stata_tol$coef,
               info = "rss")
  expect_equal(fit$r.squared, diag$r2, tolerance = stata_tol$coef,
               info = "r2")
  expect_equal(fit$nobs, diag$N, info = "N")
})


# ============================================================================
# 2. Robust GMM2S — coefficients differ from 2SLS
# ============================================================================

test_that("GMM2S robust coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_coef_robust.csv")
  skip_if(!file.exists(fixture_path), "GMM2S robust fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "HC0")
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("GMM2S robust coefficients differ from 2SLS", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "HC0")
  fit_gmm <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                    data = card, method = "gmm2s", vcov = "HC0")

  # Coefficients should differ for overidentified + heteroskedastic
  expect_false(isTRUE(all.equal(coef(fit_2sls)["educ"],
                                coef(fit_gmm)["educ"],
                                tolerance = 1e-4)))
})

test_that("GMM2S robust diagnostics match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_diagnostics_robust.csv")
  skip_if(!file.exists(fixture_path), "GMM2S robust diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "HC0")
  diag <- read.csv(fixture_path)

  # Overid (Hansen J for robust)
  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = stata_tol$stat, info = "Hansen J stat")
  expect_equal(fit$diagnostics$overid$p, diag$overid_p,
               tolerance = stata_tol$pval, info = "Hansen J p")
  expect_equal(fit$diagnostics$overid$df, diag$overid_df, info = "Hansen J df")
  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")

  # Underid (KP rk LM for robust)
  expect_equal(fit$diagnostics$underid$stat, diag$underid_stat,
               tolerance = stata_tol$stat, info = "KP rk LM")

  # Weak id
  expect_equal(fit$diagnostics$weak_id$stat, diag$weak_id_cd_f,
               tolerance = stata_tol$stat, info = "CD F")
  expect_equal(fit$diagnostics$weak_id_robust$stat, diag$weak_id_kp_f,
               tolerance = stata_tol$stat, info = "KP rk Wald F")

  # Summary stats
  expect_equal(fit$sigma, diag$sigma, tolerance = stata_tol$coef)
  expect_equal(fit$rss, diag$rss, tolerance = stata_tol$coef)
  expect_equal(fit$r.squared, diag$r2, tolerance = stata_tol$coef)

  # Model F
  expect_equal(fit$model_f, diag$model_f, tolerance = stata_tol$stat,
               info = "model F")
})


# ============================================================================
# 3. GMM2S robust small
# ============================================================================

test_that("GMM2S robust small coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_coef_robust_small.csv")
  skip_if(!file.exists(fixture_path), "GMM2S robust small fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "HC1", small = TRUE)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  # Coefficients same as non-small (small only affects SEs)
  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("GMM2S robust small diagnostics match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_diagnostics_robust_small.csv")
  skip_if(!file.exists(fixture_path), "GMM2S robust small diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "HC1", small = TRUE)
  diag <- read.csv(fixture_path)

  expect_equal(fit$sigma, diag$sigma, tolerance = stata_tol$coef, info = "sigma")
  expect_equal(fit$model_f, diag$model_f, tolerance = stata_tol$stat,
               info = "model F")
})


# ============================================================================
# 4. IID small
# ============================================================================

test_that("GMM2S IID small SEs match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_coef_iid_small.csv")
  skip_if(!file.exists(fixture_path), "GMM2S IID small fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", small = TRUE)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})


# ============================================================================
# 5. Cluster GMM2S
# ============================================================================

test_that("GMM2S cluster coefficients and SEs match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_coef_cluster.csv")
  skip_if(!file.exists(fixture_path), "GMM2S cluster fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", clusters = ~age)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("GMM2S cluster diagnostics match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_diagnostics_cluster.csv")
  skip_if(!file.exists(fixture_path), "GMM2S cluster diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", clusters = ~age)
  diag <- read.csv(fixture_path)

  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = stata_tol$stat, info = "J stat")
  expect_equal(fit$diagnostics$overid$p, diag$overid_p,
               tolerance = stata_tol$pval, info = "J p")
  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
  expect_equal(fit$sigma, diag$sigma, tolerance = stata_tol$coef)
  expect_equal(fit$rss, diag$rss, tolerance = stata_tol$coef)
})

test_that("GMM2S cluster small SEs match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_coef_cluster_small.csv")
  skip_if(!file.exists(fixture_path), "GMM2S cluster small fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", clusters = ~age, small = TRUE)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})


# ============================================================================
# 6. Just-identified GMM2S (should equal 2SLS)
# ============================================================================

test_that("Just-identified GMM2S equals 2SLS coefficients", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, vcov = "HC0")
  fit_gmm <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card, method = "gmm2s", vcov = "HC0")

  # Just-identified: GMM2S = 2SLS regardless of omega
  expect_equal(coef(fit_gmm), coef(fit_2sls), tolerance = 1e-10)
})

test_that("Just-identified GMM2S matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_justid_gmm2s_coef_robust.csv")
  skip_if(!file.exists(fixture_path), "GMM2S justid fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, method = "gmm2s", vcov = "HC0")
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
  }
})

test_that("Just-identified GMM2S J stat is zero", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, method = "gmm2s", vcov = "HC0")

  expect_equal(fit$diagnostics$overid$stat, 0)
  expect_equal(fit$diagnostics$overid$df, 0L)
  expect_true(is.na(fit$diagnostics$overid$p))
})


# ============================================================================
# 7. Weighted GMM2S
# ============================================================================

test_that("GMM2S aweight robust matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_weighted_coef_aw_robust.csv")
  skip_if(!file.exists(fixture_path), "GMM2S aweight fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "HC0", weights = age)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("GMM2S aweight cluster matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_weighted_coef_aw_cluster.csv")
  skip_if(!file.exists(fixture_path), "GMM2S aweight cluster fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", clusters = ~age, weights = weight)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("GMM2S pweight robust matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_weighted_coef_pw_robust.csv")
  skip_if(!file.exists(fixture_path), "GMM2S pweight fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", weight_type = "pweight", weights = age)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})


# ============================================================================
# 8. dofminus threading
# ============================================================================

test_that("GMM2S dofminus matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_dofminus_coef_robust.csv")
  skip_if(!file.exists(fixture_path), "GMM2S dofminus fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "HC0", dofminus = 2L)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("GMM2S dofminus diagnostics match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_dofminus_diagnostics_robust.csv")
  skip_if(!file.exists(fixture_path), "GMM2S dofminus diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "HC0", dofminus = 2L)
  diag <- read.csv(fixture_path)

  expect_equal(fit$sigma, diag$sigma, tolerance = stata_tol$coef,
               info = "sigma with dofminus")
  expect_equal(fit$rss, diag$rss, tolerance = stata_tol$coef,
               info = "rss with dofminus")
  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = stata_tol$stat, info = "J stat with dofminus")
})


# ============================================================================
# 9. HAC GMM2S
# ============================================================================

test_that("GMM2S HAC Bartlett bw=3 matches Stata fixture", {
  ts_path <- file.path(fixture_dir, "ts_gmm_data.csv")
  skip_if(!file.exists(ts_path), "ts_gmm_data not found")
  fixture_path <- file.path(fixture_dir, "ts_gmm2s_coef_hac_bartlett_bw3.csv")
  skip_if(!file.exists(fixture_path), "GMM2S HAC fixture not found")

  ts_data <- read.csv(ts_path)
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                method = "gmm2s", vcov = "HC0",
                kernel = "bartlett", bw = 3, tvar = "t")
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("GMM2S HAC diagnostics match Stata fixture", {
  ts_path <- file.path(fixture_dir, "ts_gmm_data.csv")
  skip_if(!file.exists(ts_path), "ts_gmm_data not found")
  fixture_path <- file.path(fixture_dir, "ts_gmm2s_diagnostics_hac_bartlett_bw3.csv")
  skip_if(!file.exists(fixture_path), "GMM2S HAC diagnostics fixture not found")

  ts_data <- read.csv(ts_path)
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                method = "gmm2s", vcov = "HC0",
                kernel = "bartlett", bw = 3, tvar = "t")
  diag <- read.csv(fixture_path)

  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = stata_tol$stat, info = "J stat HAC")
  expect_equal(fit$sigma, diag$sigma, tolerance = stata_tol$coef,
               info = "sigma HAC")
})


# ============================================================================
# 10. Endogeneity test
# ============================================================================

test_that("GMM2S endogeneity test matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- file.path(fixture_dir, "card_overid_gmm2s_endog_diagnostics_robust.csv")
  skip_if(!file.exists(fixture_path), "GMM2S endog fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "HC0", endog = "educ")
  diag <- read.csv(fixture_path)

  expect_equal(fit$diagnostics$endogeneity$stat, diag$endog_stat,
               tolerance = stata_tol$stat, info = "endog stat")
  expect_equal(fit$diagnostics$endogeneity$p, diag$endog_p,
               tolerance = stata_tol$pval, info = "endog p")
  expect_equal(fit$diagnostics$endogeneity$df, diag$endog_df,
               info = "endog df")
})


# ============================================================================
# 11. VCV matrix parity
# ============================================================================

test_that("GMM2S robust VCV matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  vcov_path <- file.path(fixture_dir, "card_overid_gmm2s_vcov_robust.csv")
  coef_path <- file.path(fixture_dir, "card_overid_gmm2s_coef_robust.csv")
  skip_if(!file.exists(vcov_path), "GMM2S VCV fixture not found")
  skip_if(!file.exists(coef_path), "GMM2S coef fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "HC0")
  V_stata <- unname(as.matrix(read.csv(vcov_path)))
  V_r <- vcov(fit)

  # The Stata VCV fixture uses v1..vK columns without term names.
 # Get Stata variable ordering from the coefficient fixture.
  coef_fixture <- read.csv(coef_path)
  stata_order <- ifelse(coef_fixture$term == "_cons", "(Intercept)", coef_fixture$term)

  # Reorder R's VCV to match Stata's variable ordering
  V_r_reordered <- V_r[stata_order, stata_order]

  # Match dimensions
  expect_equal(nrow(V_r_reordered), nrow(V_stata))
  expect_equal(ncol(V_r_reordered), ncol(V_stata))

  # Element-wise comparison
  for (i in seq_len(nrow(V_r_reordered))) {
    for (j in seq_len(ncol(V_r_reordered))) {
      expect_equal(unname(V_r_reordered[i, j]), V_stata[i, j],
                   tolerance = stata_tol$vcov,
                   info = paste("VCV mismatch:", i, j))
    }
  }
})


# ============================================================================
# 12. Input validation
# ============================================================================

test_that("GMM2S requires IV model", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, method = "gmm2s"),
    "requires an IV model"
  )
})

test_that("GMM2S incompatible with fuller", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "gmm2s", fuller = 1),
    "Cannot specify.*fuller.*gmm2s"
  )
})

test_that("GMM2S incompatible with kclass", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "gmm2s", kclass = 0.5),
    "Cannot specify.*kclass.*gmm2s"
  )
})

test_that("GMM2S method stored in return object", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "HC0")
  expect_equal(fit$method, "gmm2s")
})


# ============================================================================
# 13. Display
# ============================================================================

test_that("GMM2S print shows correct estimation label", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "HC0")
  output <- capture.output(print(summary(fit)))
  expect_true(any(grepl("2-Step GMM Estimation", output)))
  expect_true(any(grepl("heteroskedasticity", output)))
})


# ============================================================================
# 14. Stock-Yogo tables (should use IV tables for gmm2s)
# ============================================================================

test_that("GMM2S gets IV Stock-Yogo tables", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s")
  sy <- fit$diagnostics$weak_id_sy
  expect_false(is.null(sy))
  expect_true("IV size" %in% sy$type)
})
