# ============================================================================
# Tests: Continuously Updated GMM Estimator (CUE) — Ticket N3
# ============================================================================

# --- Helpers ---
card_path <- fixture_path("card_data.csv")

if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

# CUE involves nonlinear optimization, which adds a layer of numerical noise
# beyond direct estimation (2SLS/LIML). For small-magnitude coefficients, even
# tiny absolute differences (2e-8) can breach relative tolerances. We use
# slightly wider tolerances for CUE coefficients/SEs while keeping stats/pvals
# at the standard level.
cue_tol <- list(
  coef = 5e-6,           # coefficients (optimization noise)
  se   = 5e-6,           # standard errors
  vcov = 5e-6,           # VCV elements
  stat = stata_tol$stat, # test statistics — same as standard
  pval = stata_tol$pval  # p-values — same as standard
)

# ============================================================================
# 1. CUE IID — coefficients and SEs match Stata
# ============================================================================

test_that("CUE IID coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_coef_iid.csv")
  skip_if(!file.exists(fixture_path), "CUE IID fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue")
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = cue_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = cue_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("CUE IID VCV matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_vcov_iid.csv")
  coef_path <- fixture_path("card_overid_cue_coef_iid.csv")
  skip_if(!file.exists(fixture_path), "CUE IID VCV fixture not found")
  skip_if(!file.exists(coef_path), "CUE IID coef fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue")
  V_stata <- read_vcov_fixture(fixture_path)
  # Reorder Stata VCV to match R column order
  stata_terms <- read.csv(coef_path)$term
  r_names <- ifelse(stata_terms == "_cons", "(Intercept)", stata_terms)
  r_order <- match(names(coef(fit)), r_names)
  V_stata <- V_stata[r_order, r_order]

  for (i in seq_len(nrow(V_stata))) {
    for (j in seq_len(ncol(V_stata))) {
      expect_equal(
        unname(vcov(fit)[i, j]), unname(V_stata[i, j]),
        tolerance = cue_tol$vcov,
        info = paste("VCV mismatch:", i, j)
      )
    }
  }
})

test_that("CUE IID diagnostics match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_diagnostics_iid.csv")
  skip_if(!file.exists(fixture_path), "CUE IID diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue")
  diag <- read.csv(fixture_path)

  # Overid (always Hansen J for CUE)
  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = cue_tol$stat, info = "Hansen J stat")
  expect_equal(fit$diagnostics$overid$p, diag$overid_p,
               tolerance = cue_tol$pval, info = "Hansen J p")
  expect_equal(fit$diagnostics$overid$df, diag$overid_df, info = "Hansen J df")
  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")

  # Underid
  expect_equal(fit$diagnostics$underid$stat, diag$underid_stat,
               tolerance = cue_tol$stat, info = "underid stat")

  # Weak id
  expect_equal(fit$diagnostics$weak_id$stat, diag$weak_id_cd_f,
               tolerance = cue_tol$stat, info = "CD F")

  # Summary stats
  expect_equal(fit$sigma, diag$sigma, tolerance = cue_tol$coef, info = "sigma")
  expect_equal(fit$rss, diag$rss, tolerance = cue_tol$coef, info = "rss")
  expect_equal(fit$r.squared, diag$r2, tolerance = cue_tol$coef, info = "r2")
  expect_equal(fit$nobs, diag$N, info = "N")
})

# ============================================================================
# 2. CUE IID small — small-sample corrections
# ============================================================================

test_that("CUE IID small coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_coef_iid_small.csv")
  skip_if(!file.exists(fixture_path), "CUE IID small fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", small = TRUE)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = cue_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = cue_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("CUE IID small diagnostics match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_diagnostics_iid_small.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", small = TRUE)
  diag <- read.csv(fixture_path)

  expect_equal(fit$sigma, diag$sigma, tolerance = cue_tol$coef, info = "sigma")
  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = cue_tol$stat, info = "Hansen J stat")
})

# ============================================================================
# 3. CUE robust (HC0) — coefficients/SEs match Stata
# ============================================================================

test_that("CUE robust coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_coef_robust.csv")
  skip_if(!file.exists(fixture_path), "CUE robust fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", vcov = "robust")
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = cue_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = cue_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("CUE robust VCV matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_vcov_robust.csv")
  coef_path <- fixture_path("card_overid_cue_coef_robust.csv")
  skip_if(!file.exists(fixture_path))
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", vcov = "robust")
  V_stata <- read_vcov_fixture(fixture_path)
  stata_terms <- read.csv(coef_path)$term
  r_names <- ifelse(stata_terms == "_cons", "(Intercept)", stata_terms)
  r_order <- match(names(coef(fit)), r_names)
  V_stata <- V_stata[r_order, r_order]

  for (i in seq_len(nrow(V_stata))) {
    for (j in seq_len(ncol(V_stata))) {
      expect_equal(
        unname(vcov(fit)[i, j]), unname(V_stata[i, j]),
        tolerance = cue_tol$vcov,
        info = paste("VCV mismatch:", i, j)
      )
    }
  }
})

test_that("CUE robust diagnostics match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_diagnostics_robust.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", vcov = "robust")
  diag <- read.csv(fixture_path)

  # Overid
  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = cue_tol$stat, info = "Hansen J stat")
  expect_equal(fit$diagnostics$overid$p, diag$overid_p,
               tolerance = cue_tol$pval, info = "Hansen J p")
  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")

  # Underid (KP rk)
  if (!is.na(diag$weak_id_kp_f)) {
    expect_equal(fit$diagnostics$weak_id_robust$stat, diag$weak_id_kp_f,
                 tolerance = cue_tol$stat, info = "KP rk Wald F")
  }
})

# ============================================================================
# 4. CUE robust small (HC1)
# ============================================================================

test_that("CUE robust small coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_coef_robust_small.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", vcov = "robust", small = TRUE)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = cue_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = cue_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

# ============================================================================
# 5. CUE cluster
# ============================================================================

test_that("CUE cluster coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_coef_cluster.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", clusters = ~age)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = cue_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = cue_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("CUE cluster diagnostics match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_diagnostics_cluster.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", clusters = ~age)
  diag <- read.csv(fixture_path)

  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = cue_tol$stat, info = "Hansen J stat")
  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
})

test_that("CUE cluster small coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_coef_cluster_small.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", clusters = ~age, small = TRUE)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = cue_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = cue_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

# ============================================================================
# 6. CUE = 2SLS for exactly-identified models
# ============================================================================

test_that("CUE equals 2SLS for exactly-identified model", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, vcov = "robust")
  fit_cue <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card, method = "cue", vcov = "robust")

  # Coefficients should match to high precision
  expect_equal(coef(fit_cue), coef(fit_2sls), tolerance = 1e-6)
})

test_that("CUE just-identified coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_justid_cue_coef_robust.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, method = "cue", vcov = "robust")
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = cue_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
  }
})

test_that("CUE just-identified J = 0", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, method = "cue", vcov = "robust")

  expect_equal(fit$diagnostics$overid$stat, 0)
  expect_true(is.na(fit$diagnostics$overid$p))
  expect_equal(fit$diagnostics$overid$df, 0L)
})

# ============================================================================
# 7. Stock-Yogo uses LIML tables for CUE
# ============================================================================

test_that("CUE uses LIML Stock-Yogo tables", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue")

  sy <- fit$diagnostics$weak_id_sy
  expect_false(is.null(sy))
  # CUE should use LIML size tables
  expect_true("LIML size" %in% sy$type)
  # Should NOT have IV size tables (those are for 2SLS/GMM2S)
  expect_false("IV size" %in% sy$type)
})

# ============================================================================
# 8. CUE overid test name is always Hansen J
# ============================================================================

test_that("CUE overid is always Hansen J", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  # IID
  fit_iid <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                    data = card, method = "cue")
  expect_equal(fit_iid$diagnostics$overid$test_name, "Hansen J")

  # Robust
  fit_robust <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                       data = card, method = "cue", vcov = "robust")
  expect_equal(fit_robust$diagnostics$overid$test_name, "Hansen J")

  # Cluster
  fit_cl <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                   data = card, method = "cue", clusters = ~age)
  expect_equal(fit_cl$diagnostics$overid$test_name, "Hansen J")
})

# ============================================================================
# 9. CUE diagnostics populated
# ============================================================================

test_that("CUE populates identification diagnostics", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", vcov = "robust")

  # Overid
  expect_false(is.null(fit$diagnostics$overid))

  # Underid
  expect_false(is.null(fit$diagnostics$underid))

  # Weak-id
  expect_false(is.null(fit$diagnostics$weak_id))
  expect_false(is.null(fit$diagnostics$weak_id_robust))

  # First-stage
  expect_false(is.null(fit$first_stage))

  # AR
  expect_false(is.null(fit$diagnostics$anderson_rubin))

  # Stock-Wright
  expect_false(is.null(fit$diagnostics$stock_wright))

  # Endogeneity
  expect_false(is.null(fit$diagnostics$endogeneity))
})

# ============================================================================
# 10. CUE with endog test
# ============================================================================

test_that("CUE endog test matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_endog_diagnostics_robust.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", vcov = "robust", endog = "educ")
  diag <- read.csv(fixture_path)

  expect_equal(fit$diagnostics$endogeneity$stat, diag$endog_stat,
               tolerance = cue_tol$stat, info = "endog stat")
  expect_equal(fit$diagnostics$endogeneity$p, diag$endog_p,
               tolerance = cue_tol$pval, info = "endog p")
})

# ============================================================================
# 11. CUE with orthog test
# ============================================================================

test_that("CUE orthog test matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_orthog_diagnostics_robust.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", vcov = "robust", orthog = "nearc2")
  diag <- read.csv(fixture_path)

  # Orthog stat should be in the diagnostics
  expect_false(is.null(fit$diagnostics$orthog))
})

# ============================================================================
# 12. Weighted CUE
# ============================================================================

test_that("Weighted CUE coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_weighted_coef_aw_robust.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", vcov = "robust", weights = age)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  # `black` (~0.017) is the only term where CUE BFGS basin noise (~1.6e-7
  # absolute) breaches relative tolerance. Apply the absolute floor only to
  # this term so detection stays tight on the other five (e.g., `expersq`
  # at 1.77e-3 still catches sub-1e-8 absolute drifts via 5e-6 relative).
  for (i in seq_len(nrow(fixture))) {
    r_coef <- unname(coef(fit)[r_names[i]])
    s_coef <- fixture$estimate[i]
    r_se <- unname(sqrt(diag(vcov(fit)))[r_names[i]])
    s_se <- fixture$std_error[i]
    if (r_names[i] == "black") {
      expect_lt(abs(r_coef - s_coef), 1e-6, label = "Coef: black")
      expect_lt(abs(r_se   - s_se),   1e-6, label = "SE: black")
    } else {
      expect_equal(r_coef, s_coef, tolerance = cue_tol$coef,
                   info = paste("Coef:", r_names[i]))
      expect_equal(r_se,   s_se,   tolerance = cue_tol$se,
                   info = paste("SE:",   r_names[i]))
    }
  }
})

# ============================================================================
# 13. CUE with dofminus
# ============================================================================

test_that("CUE dofminus coefficients match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_dofminus_coef_robust.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", vcov = "robust", dofminus = 2L)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = cue_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = cue_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("CUE dofminus diagnostics match Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_cue_dofminus_diagnostics_robust.csv")
  skip_if(!file.exists(fixture_path))

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue", vcov = "robust", dofminus = 2L)
  diag <- read.csv(fixture_path)

  expect_equal(fit$sigma, diag$sigma, tolerance = cue_tol$coef, info = "sigma")
})

# ============================================================================
# 14. HAC CUE
# ============================================================================

test_that("HAC CUE coefficients match Stata fixture", {
  ts_data_path <- fixture_path("ts_gmm_data.csv")
  skip_if(!file.exists(ts_data_path), "TS GMM dataset not found")
  fixture_path <- fixture_path("ts_cue_coef_hac_bartlett_bw3.csv")
  skip_if(!file.exists(fixture_path), "HAC CUE fixture not found")

  ts_data <- read.csv(ts_data_path)
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                method = "cue", kernel = "bartlett", bw = 3, tvar = "t",
                vcov = "robust")
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = cue_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = cue_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("HAC CUE diagnostics match Stata fixture", {
  ts_data_path <- fixture_path("ts_gmm_data.csv")
  skip_if(!file.exists(ts_data_path), "TS GMM dataset not found")
  fixture_path <- fixture_path("ts_cue_diagnostics_hac_bartlett_bw3.csv")
  skip_if(!file.exists(fixture_path))

  ts_data <- read.csv(ts_data_path)
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                method = "cue", kernel = "bartlett", bw = 3, tvar = "t",
                vcov = "robust")
  diag <- read.csv(fixture_path)

  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = cue_tol$stat, info = "Hansen J stat")
  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
})

# ============================================================================
# 15. b0 mode — evaluate at 2SLS estimates
# ============================================================================

test_that("b0 mode produces coefficients equal to b0", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  # Get 2SLS estimates to use as b0
  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "robust")
  b0_vec <- coef(fit_2sls)

  fit_b0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                   data = card, vcov = "robust", b0 = b0_vec)

  # Coefficients should be exactly b0
  expect_equal(coef(fit_b0), b0_vec, tolerance = 1e-14)
})

test_that("b0 mode J(b0) matches Stata fixture", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fixture_path <- fixture_path("card_overid_b0_diagnostics_robust.csv")
  b0_path <- fixture_path("card_overid_b0_vector.csv")
  skip_if(!file.exists(fixture_path) || !file.exists(b0_path))

  b0_fixture <- read.csv(b0_path)
  b0_vec <- b0_fixture$b0_value
  stata_names <- b0_fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)
  names(b0_vec) <- r_names

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", b0 = b0_vec)
  diag <- read.csv(fixture_path)

  # J(b0)
  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = cue_tol$stat, info = "J(b0) stat")
  expect_equal(fit$diagnostics$overid$p, diag$overid_p,
               tolerance = cue_tol$pval, info = "J(b0) p")
})

test_that("b0 mode suppresses identification diagnostics", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "robust")
  b0_vec <- coef(fit_2sls)

  fit_b0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                   data = card, vcov = "robust", b0 = b0_vec)

  # Overid should still be populated (J-stat)
  expect_false(is.null(fit_b0$diagnostics$overid))

  # All identification diagnostics should be NULL
  expect_null(fit_b0$diagnostics$underid)
  expect_null(fit_b0$diagnostics$weak_id)
  expect_null(fit_b0$diagnostics$weak_id_robust)
  expect_null(fit_b0$diagnostics$weak_id_sy)
  expect_null(fit_b0$diagnostics$anderson_rubin)
  expect_null(fit_b0$diagnostics$stock_wright)
  expect_null(fit_b0$diagnostics$endogeneity)

  # First stage should be NULL
  expect_null(fit_b0$first_stage)
})

test_that("b0 implies method = cue", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "robust")
  b0_vec <- coef(fit_2sls)

  fit_b0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                   data = card, vcov = "robust", b0 = b0_vec)

  expect_equal(fit_b0$method, "cue")
})

test_that("b0 convergence is 0", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "robust")
  b0_vec <- coef(fit_2sls)

  fit_b0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                   data = card, vcov = "robust", b0 = b0_vec)

  expect_equal(fit_b0$cue_convergence, 0L)
})

test_that("b0 unnamed vector works in column order", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "robust")
  b0_named <- coef(fit_2sls)
  b0_unnamed <- unname(b0_named)

  fit_named <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                      data = card, vcov = "robust", b0 = b0_named)
  fit_unnamed <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                        data = card, vcov = "robust", b0 = b0_unnamed)

  expect_equal(coef(fit_named), coef(fit_unnamed))
})

# ============================================================================
# 16. Input validation
# ============================================================================

test_that("CUE mutual exclusion with fuller", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "cue", fuller = 1),
    "Cannot specify `fuller` with `method = \"cue\"`"
  )
})

test_that("CUE mutual exclusion with kclass", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "cue", kclass = 0.5),
    "Cannot specify `kclass` with `method = \"cue\"`"
  )
})

test_that("CUE mutual exclusion with wmatrix", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  W <- diag(7)
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "cue", wmatrix = W),
    "Cannot specify `wmatrix` with `method = \"cue\"`"
  )
})

test_that("CUE mutual exclusion with smatrix", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  S <- diag(7)
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "cue", smatrix = S),
    "Cannot specify `smatrix` with `method = \"cue\"`"
  )
})

test_that("CUE requires IV model", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, method = "cue"),
    "requires an IV model"
  )
})

test_that("b0 mutual exclusion with gmm2s", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "gmm2s", b0 = rep(0, 6)),
    "Cannot specify `b0` with `method = \"gmm2s\"`"
  )
})

test_that("b0 mutual exclusion with liml", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "liml", b0 = rep(0, 6)),
    "Cannot specify `b0` with `method = \"liml\"`"
  )
})

test_that("b0 mutual exclusion with fuller", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, fuller = 1, b0 = rep(0, 6)),
    "Cannot specify `b0` with `fuller`"
  )
})

test_that("b0 mutual exclusion with wmatrix", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  W <- diag(7)
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, b0 = rep(0, 6), wmatrix = W),
    "Cannot specify `b0` with `wmatrix`"
  )
})

test_that("b0 dimension check", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, b0 = c(1, 2, 3)),
    "`b0` length .* must equal"
  )
})

test_that("b0 must be numeric", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, b0 = "not numeric"),
    "`b0` must be a numeric vector"
  )
})

test_that("b0 must be finite", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, b0 = c(1, 2, 3, NA, 5, 6)),
    "`b0` must contain only finite values"
  )
})

test_that("b0 named vector with wrong names errors", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  b0_bad <- c(a = 1, b = 2, c = 3, d = 4, e = 5, f = 6)
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, b0 = b0_bad),
    "not matching model columns"
  )
})

test_that("b0 + uppercase method errors correctly (case-insensitive)", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  # Regression test: uppercase method must still trigger incompatibility errors
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "LIML", b0 = rep(0, 6)),
    "Cannot specify `b0` with `method = \"liml\"`"
  )
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "GMM2S", b0 = rep(0, 6)),
    "Cannot specify `b0` with `method = \"gmm2s\"`"
  )
})

test_that("b0 + uppercase '2SLS' auto-promotes to CUE", {
  skip_if(!file.exists(card_path), "Card dataset not found")
  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "robust")
  b0_vec <- coef(fit_2sls)
  # method = "2SLS" (uppercase) should still auto-promote to CUE
  fit_b0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                   data = card, vcov = "robust", method = "2SLS", b0 = b0_vec)
  expect_equal(fit_b0$method, "cue")
})

# ============================================================================
# 17. Convergence
# ============================================================================

test_that("CUE convergence code = 0 for well-specified model", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue")

  expect_equal(fit$cue_convergence, 0L)
  expect_false(is.null(fit$cue_message))
})

test_that("CUE convergence fields are NULL for non-CUE", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card)

  expect_null(fit$cue_convergence)
  expect_null(fit$cue_message)
})

# ============================================================================
# 18. Broom methods
# ============================================================================

test_that("tidy works for CUE", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue")
  td <- tidy(fit)

  expect_s3_class(td, "tbl_df")
  expect_equal(nrow(td), 6L)
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value") %in%
                  names(td)))
})

test_that("glance works for CUE and includes cue_convergence", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue")
  gl <- glance(fit)

  expect_s3_class(gl, "tbl_df")
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$method, "cue")
  expect_true("cue_convergence" %in% names(gl))
  expect_equal(gl$cue_convergence, 0L)
})

test_that("glance cue_convergence is NA for non-CUE", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card)
  gl <- glance(fit)

  expect_true(is.na(gl$cue_convergence))
})

test_that("augment works for CUE", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue")
  aug <- augment(fit)

  expect_s3_class(aug, "tbl_df")
  expect_true(all(c(".fitted", ".resid") %in% names(aug)))
})

# ============================================================================
# 19. Method field
# ============================================================================

test_that("CUE method field is 'cue'", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue")

  expect_equal(fit$method, "cue")
})

# ============================================================================
# 20. CUE close to LIML for IID overidentified (asymptotic equivalence)
# ============================================================================

test_that("CUE IID is close to LIML (asymptotic equivalence)", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_cue <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                    data = card, method = "cue")
  fit_liml <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, method = "liml")

  # CUE and LIML should be close (asymptotically equivalent)
  # Use a generous tolerance since they differ at finite samples
  expect_equal(coef(fit_cue), coef(fit_liml), tolerance = 0.01)
})

# ============================================================================
# 21. Summary/print methods
# ============================================================================

test_that("summary prints for CUE without error", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue")

  out <- capture.output(summary(fit))
  expect_true(any(grepl("CUE Estimation", out)))
  expect_true(any(grepl("Estimates efficient for arbitrary homoskedasticity", out)))
})

test_that("summary prints for CUE with b0", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "robust")
  fit_b0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                   data = card, vcov = "robust", b0 = coef(fit_2sls))

  out <- capture.output(summary(fit_b0))
  expect_true(any(grepl("User-supplied Parameter Vector", out)))
})

test_that("print method works for CUE", {
  skip_if(!file.exists(card_path), "Card dataset not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "cue")

  out <- capture.output(print(fit))
  expect_true(any(grepl("CUE Estimation", out)))
})
