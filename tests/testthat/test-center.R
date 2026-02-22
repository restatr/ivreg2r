# ============================================================================
# Tests: Center Option — Ticket N4
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
  as.matrix(read.csv(path))
}

# CUE tolerance (optimization noise; centering changes the GMM objective,
# so the optimizer may converge to a slightly different point than Stata's)
cue_tol <- list(
  coef = 2e-5,
  se   = 2e-5,
  vcov = 2e-5,
  stat = stata_tol$stat,
  pval = stata_tol$pval
)

# Helper: compare coefficients and SEs against fixture
check_coef_fixture <- function(fit, fixture_path, tol = stata_tol) {
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
}

# Helper: compare VCV against fixture
check_vcov_fixture <- function(fit, vcov_path, coef_path, tol = stata_tol) {
  V_stata <- read_vcov_fixture(vcov_path)
  stata_terms <- read.csv(coef_path)$term
  r_names <- ifelse(stata_terms == "_cons", "(Intercept)", stata_terms)
  r_order <- match(names(coef(fit)), r_names)
  V_stata <- V_stata[r_order, r_order]

  for (i in seq_len(nrow(V_stata))) {
    for (j in seq_len(ncol(V_stata))) {
      expect_equal(
        unname(vcov(fit)[i, j]), unname(V_stata[i, j]),
        tolerance = tol$vcov,
        info = paste0("VCV[", i, ",", j, "] mismatch")
      )
    }
  }
}

# Helper: compare diagnostics against fixture
check_diag_fixture <- function(fit, fixture_path, tol = stata_tol) {
  fixture <- read.csv(fixture_path)
  diag <- fit$diagnostics

  # Overid
  if (!is.null(diag$overid) && diag$overid$df > 0L && !is.na(fixture$overid_stat)) {
    expect_equal(diag$overid$stat, fixture$overid_stat,
                 tolerance = tol$stat, info = "overid stat")
    expect_equal(diag$overid$p, fixture$overid_p,
                 tolerance = tol$pval, info = "overid p")
  }

  # Underid
  if (!is.null(diag$underid) && !is.na(fixture$underid_stat)) {
    expect_equal(diag$underid$stat, fixture$underid_stat,
                 tolerance = tol$stat, info = "underid stat")
    expect_equal(diag$underid$p, fixture$underid_p,
                 tolerance = tol$pval, info = "underid p")
  }

  # Weak ID (Cragg-Donald)
  if (!is.null(diag$weak_id) && !is.na(fixture$weak_id_cd_f)) {
    expect_equal(diag$weak_id$stat, fixture$weak_id_cd_f,
                 tolerance = tol$stat, info = "CD F")
  }

  # Weak ID robust (KP)
  if (!is.null(diag$weak_id_robust) && !is.na(fixture$weak_id_kp_f)) {
    expect_equal(diag$weak_id_robust$stat, fixture$weak_id_kp_f,
                 tolerance = tol$stat, info = "KP F")
  }

  # Anderson-Rubin
  if (!is.null(diag$anderson_rubin) && !is.na(fixture$ar_f)) {
    expect_equal(diag$anderson_rubin$f_stat, fixture$ar_f,
                 tolerance = tol$stat, info = "AR F")
    expect_equal(diag$anderson_rubin$f_p, fixture$ar_f_p,
                 tolerance = tol$pval, info = "AR F p")
    expect_equal(diag$anderson_rubin$chi2_stat, fixture$ar_chi2,
                 tolerance = tol$stat, info = "AR chi2")
  }

  # Stock-Wright
  if (!is.null(diag$stock_wright) && !is.na(fixture$sw_stat)) {
    expect_equal(diag$stock_wright$stat, fixture$sw_stat,
                 tolerance = tol$stat, info = "SW stat")
    expect_equal(diag$stock_wright$p, fixture$sw_p,
                 tolerance = tol$pval, info = "SW p")
  }

  # Endogeneity
  if (!is.null(diag$endogeneity) && !is.na(fixture$endog_stat)) {
    expect_equal(diag$endogeneity$stat, fixture$endog_stat,
                 tolerance = tol$stat, info = "endogeneity stat")
    expect_equal(diag$endogeneity$p, fixture$endog_p,
                 tolerance = tol$pval, info = "endogeneity p")
  }

  # Model F
  if (!is.null(fit$model_f) && !is.na(fixture$model_f)) {
    expect_equal(fit$model_f, fixture$model_f,
                 tolerance = tol$stat, info = "model F")
  }

  # Sigma
  if (!is.na(fixture$sigma)) {
    expect_equal(fit$sigma, fixture$sigma,
                 tolerance = tol$coef, info = "sigma")
  }
}

# ============================================================================
# 1. Input validation
# ============================================================================

test_that("center must be TRUE or FALSE", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, center = 1),
    "`center` must be TRUE or FALSE"
  )
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, center = "yes"),
    "`center` must be TRUE or FALSE"
  )
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, center = NA),
    "`center` must be TRUE or FALSE"
  )
})


# ============================================================================
# 2. Warning for no-op configurations
# ============================================================================

test_that("center = TRUE + IID gives warning", {
  expect_warning(
    ivreg2(mpg ~ wt + hp, data = mtcars, center = TRUE),
    "center.*has no effect.*iid"
  )
})

test_that("center = TRUE + HC0 gives no warning", {
  expect_no_warning(
    ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC0", center = TRUE)
  )
})

test_that("center = FALSE gives no warning (iid)", {
  expect_no_warning(
    ivreg2(mpg ~ wt + hp, data = mtcars, center = FALSE)
  )
})


# ============================================================================
# 3. Regression: center = FALSE identical to omitting center
# ============================================================================

test_that("center = FALSE gives identical results to default", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_default <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC0"
  )
  fit_false <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC0", center = FALSE
  )

  expect_identical(coef(fit_default), coef(fit_false))
  expect_identical(vcov(fit_default), vcov(fit_false))
  expect_identical(fit_default$sigma, fit_false$sigma)
  expect_identical(fit_default$diagnostics$overid$stat,
                   fit_false$diagnostics$overid$stat)
})


# ============================================================================
# 4. Coefficient invariance: center affects VCE only, not point estimates
# ============================================================================

test_that("center does not change 2SLS coefficients", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_no <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC0", center = FALSE
  )
  fit_yes <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC0", center = TRUE
  )

  expect_identical(coef(fit_no), coef(fit_yes))
  expect_identical(fit_no$sigma, fit_yes$sigma)
  expect_identical(fit_no$rss, fit_yes$rss)
})


# ============================================================================
# 5. Storage and exposure
# ============================================================================

test_that("center is stored in fit object", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC0", center = TRUE)
  expect_true(fit$center)

  fit2 <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_false(fit2$center)
})

test_that("glance includes center column", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "HC0", center = TRUE)
  gl <- glance(fit)
  expect_true("center" %in% names(gl))
  expect_true(gl$center)
})


# ============================================================================
# 6. HC0 + center: VCE matches Stata
# ============================================================================

test_that("HC0 + center: coefficients and SEs match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fp <- file.path(fixture_dir, "card_overid_coef_center_hc0.csv")
  skip_if(!file.exists(fp), "HC0 center fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC0", center = TRUE
  )
  check_coef_fixture(fit, fp)
})

test_that("HC0 + center: VCV matches Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  vp <- file.path(fixture_dir, "card_overid_vcov_center_hc0.csv")
  cp <- file.path(fixture_dir, "card_overid_coef_center_hc0.csv")
  skip_if(!file.exists(vp), "HC0 center VCV fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC0", center = TRUE
  )
  check_vcov_fixture(fit, vp, cp)
})

test_that("HC0 + center: diagnostics match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  dp <- file.path(fixture_dir, "card_overid_diagnostics_center_hc0.csv")
  skip_if(!file.exists(dp), "HC0 center diag fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC0", center = TRUE
  )
  check_diag_fixture(fit, dp)
})


# ============================================================================
# 7. HC1 + small + center: VCE matches Stata
# ============================================================================

test_that("HC1 + small + center: coefs and SEs match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fp <- file.path(fixture_dir, "card_overid_coef_center_hc1_small.csv")
  skip_if(!file.exists(fp), "HC1 small center fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC1", small = TRUE, center = TRUE
  )
  check_coef_fixture(fit, fp)
})

test_that("HC1 + small + center: VCV matches Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  vp <- file.path(fixture_dir, "card_overid_vcov_center_hc1_small.csv")
  cp <- file.path(fixture_dir, "card_overid_coef_center_hc1_small.csv")
  skip_if(!file.exists(vp), "HC1 small center VCV fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC1", small = TRUE, center = TRUE
  )
  check_vcov_fixture(fit, vp, cp)
})


# ============================================================================
# 8. Cluster + center: VCE matches Stata
# ============================================================================

test_that("Cluster + center: coefs and SEs match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fp <- file.path(fixture_dir, "card_overid_coef_center_cl.csv")
  skip_if(!file.exists(fp), "Cluster center fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, clusters = ~smsa, center = TRUE
  )
  check_coef_fixture(fit, fp)
})

test_that("Cluster + center: VCV matches Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  vp <- file.path(fixture_dir, "card_overid_vcov_center_cl.csv")
  cp <- file.path(fixture_dir, "card_overid_coef_center_cl.csv")
  skip_if(!file.exists(vp), "Cluster center VCV fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, clusters = ~smsa, center = TRUE
  )
  check_vcov_fixture(fit, vp, cp)
})

test_that("Cluster + center: diagnostics match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  dp <- file.path(fixture_dir, "card_overid_diagnostics_center_cl.csv")
  skip_if(!file.exists(dp), "Cluster center diag fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, clusters = ~smsa, center = TRUE
  )
  check_diag_fixture(fit, dp)
})


# ============================================================================
# 9. Cluster + small + center
# ============================================================================

test_that("Cluster + small + center: coefs and SEs match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fp <- file.path(fixture_dir, "card_overid_coef_center_cl_small.csv")
  skip_if(!file.exists(fp), "Cluster small center fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, clusters = ~smsa, small = TRUE, center = TRUE
  )
  check_coef_fixture(fit, fp)
})


# ============================================================================
# 10. Just-identified HC0 + center
# ============================================================================

test_that("Just-identified HC0 + center: coefs and SEs match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fp <- file.path(fixture_dir, "card_justid_coef_center_hc0.csv")
  skip_if(!file.exists(fp), "Just-id center fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc4,
    data = card, vcov = "HC0", center = TRUE
  )
  check_coef_fixture(fit, fp)
})


# ============================================================================
# 11. GMM2S robust + center
# ============================================================================

test_that("GMM2S HC0 + center: coefs and SEs match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fp <- file.path(fixture_dir, "card_overid_coef_center_gmm2s_hc0.csv")
  skip_if(!file.exists(fp), "GMM2S HC0 center fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, method = "gmm2s", vcov = "HC0", center = TRUE
  )
  check_coef_fixture(fit, fp)
})

test_that("GMM2S HC0 + center: VCV matches Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  vp <- file.path(fixture_dir, "card_overid_vcov_center_gmm2s_hc0.csv")
  cp <- file.path(fixture_dir, "card_overid_coef_center_gmm2s_hc0.csv")
  skip_if(!file.exists(vp), "GMM2S HC0 center VCV fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, method = "gmm2s", vcov = "HC0", center = TRUE
  )
  check_vcov_fixture(fit, vp, cp)
})

test_that("GMM2S HC0 + center: diagnostics match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  dp <- file.path(fixture_dir, "card_overid_diagnostics_center_gmm2s_hc0.csv")
  skip_if(!file.exists(dp), "GMM2S HC0 center diag fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, method = "gmm2s", vcov = "HC0", center = TRUE
  )
  check_diag_fixture(fit, dp)
})


# ============================================================================
# 12. GMM2S cluster + center
# ============================================================================

test_that("GMM2S cluster + center: coefs and SEs match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fp <- file.path(fixture_dir, "card_overid_coef_center_gmm2s_cl.csv")
  skip_if(!file.exists(fp), "GMM2S cluster center fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, method = "gmm2s", clusters = ~age, center = TRUE
  )
  check_coef_fixture(fit, fp)
})

test_that("GMM2S cluster + center: diagnostics match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  dp <- file.path(fixture_dir, "card_overid_diagnostics_center_gmm2s_cl.csv")
  skip_if(!file.exists(dp), "GMM2S cluster center diag fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, method = "gmm2s", clusters = ~age, center = TRUE
  )
  check_diag_fixture(fit, dp)
})


# ============================================================================
# 13. CUE robust + center
# ============================================================================

test_that("CUE HC0 + center: coefs and SEs match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fp <- file.path(fixture_dir, "card_overid_coef_center_cue_hc0.csv")
  skip_if(!file.exists(fp), "CUE HC0 center fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, method = "cue", vcov = "HC0", center = TRUE
  )
  check_coef_fixture(fit, fp, tol = cue_tol)
})

test_that("CUE HC0 + center: diagnostics match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  dp <- file.path(fixture_dir, "card_overid_diagnostics_center_cue_hc0.csv")
  skip_if(!file.exists(dp), "CUE HC0 center diag fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, method = "cue", vcov = "HC0", center = TRUE
  )
  check_diag_fixture(fit, dp, tol = cue_tol)
})


# ============================================================================
# 14. CUE cluster + center
# ============================================================================

test_that("CUE cluster + center: coefs and SEs match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fp <- file.path(fixture_dir, "card_overid_coef_center_cue_cl.csv")
  skip_if(!file.exists(fp), "CUE cluster center fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, method = "cue", clusters = ~age, center = TRUE
  )
  check_coef_fixture(fit, fp, tol = cue_tol)
})


# ============================================================================
# 15. Endogeneity test + center
# ============================================================================

test_that("Endogeneity test + center: diagnostics match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  dp <- file.path(fixture_dir, "card_overid_diagnostics_center_endog.csv")
  skip_if(!file.exists(dp), "Endogeneity center fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC0", center = TRUE, endog = "educ"
  )
  check_diag_fixture(fit, dp)
})


# ============================================================================
# 16. Orthogonality test + center
# ============================================================================

test_that("Orthogonality test + center: diagnostics match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  dp <- file.path(fixture_dir, "card_overid_diagnostics_center_orthog.csv")
  skip_if(!file.exists(dp), "Orthogonality center fixture not found")

  fixture <- read.csv(dp)

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC0", center = TRUE, orthog = "nearc2"
  )

  # Orthogonality test is stored in diagnostics$orthog in the fixture
  # but extracted from e(cstat)/e(cstatp)/e(cstatdf) in Stata
  # The fixture extracts overid_stat which is the J stat, not the orthog C stat.
  # Check overall diagnostics first:
  check_diag_fixture(fit, dp)
})


# ============================================================================
# 17. dofminus + center
# ============================================================================

test_that("dofminus + center: coefs and SEs match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fp <- file.path(fixture_dir, "card_overid_coef_center_dofminus.csv")
  skip_if(!file.exists(fp), "dofminus center fixture not found")

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, vcov = "HC0", center = TRUE, dofminus = 1L
  )
  check_coef_fixture(fit, fp)
})


# ============================================================================
# 18. HAC Bartlett + center
# ============================================================================

test_that("HAC Bartlett + center: coefs and SEs match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fp <- file.path(fixture_dir, "card_overid_coef_center_hac_bartlett.csv")
  skip_if(!file.exists(fp), "HAC center fixture not found")

  card_ts <- card
  card_ts$t <- seq_len(nrow(card_ts))

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card_ts, kernel = "bartlett", bw = 3, tvar = "t",
    center = TRUE
  )
  check_coef_fixture(fit, fp)
})

test_that("HAC Bartlett + center: VCV matches Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  vp <- file.path(fixture_dir, "card_overid_vcov_center_hac_bartlett.csv")
  cp <- file.path(fixture_dir, "card_overid_coef_center_hac_bartlett.csv")
  skip_if(!file.exists(vp), "HAC center VCV fixture not found")

  card_ts <- card
  card_ts$t <- seq_len(nrow(card_ts))

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card_ts, kernel = "bartlett", bw = 3, tvar = "t",
    center = TRUE
  )
  check_vcov_fixture(fit, vp, cp)
})

test_that("HAC Bartlett + center: diagnostics match Stata", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  dp <- file.path(fixture_dir, "card_overid_diagnostics_center_hac_bartlett.csv")
  skip_if(!file.exists(dp), "HAC center diag fixture not found")

  card_ts <- card
  card_ts$t <- seq_len(nrow(card_ts))

  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card_ts, kernel = "bartlett", bw = 3, tvar = "t",
    center = TRUE
  )
  check_diag_fixture(fit, dp)
})
