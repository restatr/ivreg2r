# ============================================================================
# Tests: dofminus / sdofminus (Ticket K1)
# ============================================================================
#
# Verifies that dofminus(1) sdofminus(1) adjustments match Stata's ivreg2
# across all VCE types (iid, HC1, CL) x small for:
#   - Coefficients (unchanged — dofminus/sdofminus don't affect point estimates)
#   - Standard errors and VCV
#   - Diagnostic tests (overid, underid, weak id, AR, endogeneity)
#   - First-stage results
#   - df.residual

# --- Helpers ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)
card_path <- file.path(fixture_dir, "card_data.csv")
sim_path  <- file.path(fixture_dir, "sim_cluster_data.csv")

if (file.exists(card_path)) card <- read.csv(card_path)
if (file.exists(sim_path))  sim_cluster <- read.csv(sim_path)

# VCE configurations to loop over.
# Fixture naming: "hc1" = Stata's `robust` (= HC0 in sandwich terminology),
# "hc1_small" = Stata's `robust small` (= HC1 with correction).
vce_configs <- list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "HC0",  small = FALSE, suffix = "hc1"),
  list(vcov = "HC1",  small = TRUE,  suffix = "hc1_small"),
  list(vcov = "iid",  small = FALSE, suffix = "cl",       clusters = TRUE),
  list(vcov = "iid",  small = TRUE,  suffix = "cl_small", clusters = TRUE)
)


# ============================================================================
# Identity: dofminus=0, sdofminus=0 matches default
# ============================================================================

test_that("dofminus=0 sdofminus=0 is identical to default", {
  skip_if(!file.exists(card_path), "card data not found")

  fit_default <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                        data = card)
  fit_zero <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, dofminus = 0L, sdofminus = 0L)

  expect_identical(coef(fit_default), coef(fit_zero))
  expect_identical(fit_default$vcov, fit_zero$vcov)
  expect_identical(fit_default$sigma, fit_zero$sigma)
  expect_identical(fit_default$df.residual, fit_zero$df.residual)
  expect_identical(fit_default$diagnostics$overid$stat,
                   fit_zero$diagnostics$overid$stat)
})


# ============================================================================
# card_just_id_dofminus: coefficients & SEs (6 VCE combos)
# ============================================================================

for (cfg in vce_configs) {
  test_that(paste0("dofminus coefs match Stata: card_just_id ", cfg$suffix), {
    skip_if(!file.exists(card_path), "card data not found")
    coef_path <- file.path(fixture_dir,
                           paste0("card_just_id_dofminus_coef_", cfg$suffix, ".csv"))
    skip_if(!file.exists(coef_path), "fixture not found")

    args <- list(formula = lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = cfg$vcov, small = cfg$small,
                 dofminus = 1L, sdofminus = 1L)
    if (isTRUE(cfg$clusters)) args$clusters <- ~ smsa66
    fit <- do.call(ivreg2, args)
    fixture <- read.csv(coef_path)

    # Cluster VCV has ~5e-6 R-vs-Stata gap on Card data (pre-existing,
    # not caused by dofminus — confirmed VCV identical with/without dofminus).
    se_tol <- if (isTRUE(cfg$clusters)) 1e-5 else stata_tol$se
    for (i in seq_len(nrow(fixture))) {
      term <- fixture$term[i]
      r_name <- if (term == "_cons") "(Intercept)" else term
      expect_equal(unname(coef(fit)[r_name]), fixture$estimate[i],
                   tolerance = stata_tol$coef,
                   info = paste("coef", r_name, cfg$suffix))
      expect_equal(sqrt(fit$vcov[r_name, r_name]), fixture$std_error[i],
                   tolerance = se_tol,
                   info = paste("SE", r_name, cfg$suffix))
    }
  })
}


# ============================================================================
# card_just_id_dofminus: VCV (6 VCE combos)
# ============================================================================

for (cfg in vce_configs) {
  test_that(paste0("dofminus VCV matches Stata: card_just_id ", cfg$suffix), {
    skip_if(!file.exists(card_path), "card data not found")
    vcov_path <- file.path(fixture_dir,
                           paste0("card_just_id_dofminus_vcov_", cfg$suffix, ".csv"))
    skip_if(!file.exists(vcov_path), "fixture not found")

    args <- list(formula = lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = cfg$vcov, small = cfg$small,
                 dofminus = 1L, sdofminus = 1L)
    if (isTRUE(cfg$clusters)) args$clusters <- ~ smsa66
    fit <- do.call(ivreg2, args)
    fixture <- read.csv(vcov_path)

    stata_names <- fixture$term
    r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)
    vcov_cols <- grep("^vcov_", names(fixture), value = TRUE)
    V_stata <- as.matrix(fixture[, vcov_cols])
    rownames(V_stata) <- r_names
    col_stata <- sub("^vcov_", "", vcov_cols)
    colnames(V_stata) <- ifelse(col_stata == "_cons", "(Intercept)", col_stata)

    vcov_tol <- if (isTRUE(cfg$clusters)) 1e-5 else stata_tol$vcov
    shared <- intersect(r_names, rownames(fit$vcov))
    for (rn in shared) {
      for (cn in shared) {
        expect_equal(fit$vcov[rn, cn], V_stata[rn, cn],
                     tolerance = vcov_tol,
                     info = paste("VCV", rn, cn, cfg$suffix))
      }
    }
  })
}


# ============================================================================
# card_just_id_dofminus: diagnostics (6 VCE combos)
# ============================================================================

for (cfg in vce_configs) {
  test_that(paste0("dofminus diagnostics match Stata: card_just_id ", cfg$suffix), {
    skip_if(!file.exists(card_path), "card data not found")
    diag_path <- file.path(fixture_dir,
                           paste0("card_just_id_dofminus_diagnostics_", cfg$suffix, ".csv"))
    skip_if(!file.exists(diag_path), "fixture not found")

    args <- list(formula = lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = cfg$vcov, small = cfg$small,
                 dofminus = 1L, sdofminus = 1L)
    if (isTRUE(cfg$clusters)) args$clusters <- ~ smsa66
    fit <- do.call(ivreg2, args)
    d <- fit$diagnostics
    fx <- read.csv(diag_path)

    # Summary stats
    expect_equal(fit$sigma, sqrt(fx$sigmasq), tolerance = stata_tol$coef,
                 info = paste("sigma", cfg$suffix))
    expect_identical(fit$df.residual, as.integer(fx$F_df2),
                     info = paste("df.residual", cfg$suffix))

    # Underidentification (KP for robust/CL, Anderson for iid)
    if (!is.na(fx$idstat)) {
      expect_equal(d$underid$stat, fx$idstat, tolerance = stata_tol$stat,
                   info = paste("underid stat", cfg$suffix))
      expect_equal(d$underid$p, fx$idp, tolerance = stata_tol$pval,
                   info = paste("underid p", cfg$suffix))
    }

    # Weak identification: Cragg-Donald F
    if (!is.na(fx$cdf)) {
      expect_equal(d$weak_id$stat, fx$cdf, tolerance = stata_tol$stat,
                   info = paste("CD F", cfg$suffix))
    }

    # Weak identification: KP Wald F (robust/CL only; NULL for iid)
    if (!is.na(fx$widstat) && !is.null(d$weak_id_robust)) {
      expect_equal(d$weak_id_robust$stat, fx$widstat, tolerance = stata_tol$stat,
                   info = paste("KP Wald F", cfg$suffix))
    }

    # Anderson-Rubin F and chi2 (may be NA if RVR singular)
    if (!is.na(fx$arf) && !is.na(d$anderson_rubin$f_stat)) {
      expect_equal(d$anderson_rubin$f_stat, fx$arf, tolerance = stata_tol$stat,
                   info = paste("AR F", cfg$suffix))
      expect_equal(d$anderson_rubin$f_p, fx$arfp, tolerance = stata_tol$pval,
                   info = paste("AR F p", cfg$suffix))
    }
    if (!is.na(fx$archi2) && !is.na(d$anderson_rubin$chi2_stat)) {
      expect_equal(d$anderson_rubin$chi2_stat, fx$archi2, tolerance = stata_tol$stat,
                   info = paste("AR chi2", cfg$suffix))
      expect_equal(d$anderson_rubin$chi2_p, fx$archi2p, tolerance = stata_tol$pval,
                   info = paste("AR chi2 p", cfg$suffix))
    }

    # Endogeneity test (C-statistic)
    if (!is.na(fx$estat) && fx$estat > 0) {
      expect_equal(d$endogeneity$stat, fx$estat, tolerance = stata_tol$stat,
                   info = paste("endogeneity stat", cfg$suffix))
      expect_equal(d$endogeneity$p, fx$estatp, tolerance = stata_tol$pval,
                   info = paste("endogeneity p", cfg$suffix))
    }

    # Model F
    if (!is.na(fx$F_stat)) {
      expect_equal(fit$model_f, fx$F_stat, tolerance = stata_tol$stat,
                   info = paste("model F", cfg$suffix))
      expect_equal(fit$model_f_p, fx$F_p, tolerance = stata_tol$pval,
                   info = paste("model F p", cfg$suffix))
    }
  })
}


# ============================================================================
# card_just_id_dofminus: first-stage (6 VCE combos)
# ============================================================================

for (cfg in vce_configs) {
  test_that(paste0("dofminus first-stage matches Stata: card_just_id ", cfg$suffix), {
    skip_if(!file.exists(card_path), "card data not found")
    fs_path <- file.path(fixture_dir,
                         paste0("card_just_id_dofminus_firststage_", cfg$suffix, ".csv"))
    skip_if(!file.exists(fs_path), "fixture not found")

    args <- list(formula = lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, vcov = cfg$vcov, small = cfg$small,
                 dofminus = 1L, sdofminus = 1L)
    if (isTRUE(cfg$clusters)) args$clusters <- ~ smsa66
    fit <- do.call(ivreg2, args)
    fs <- fit$first_stage
    fx <- read.csv(fs_path)

    # Map fixture: rows are statistics, columns are endogenous variables
    fx_stats <- fx$statistic
    endo_cols <- setdiff(names(fx), "statistic")

    for (endo in endo_cols) {
      r_name <- endo

      # RMSE
      rmse_row <- which(fx_stats == "rmse")
      if (length(rmse_row) > 0 && !is.na(fx[rmse_row, endo])) {
        expect_equal(fs[[r_name]]$rmse, fx[rmse_row, endo],
                     tolerance = stata_tol$coef,
                     info = paste("FS RMSE", endo, cfg$suffix))
      }

      # Shea partial R²
      shea_row <- which(fx_stats == "sheapr2")
      if (length(shea_row) > 0 && !is.na(fx[shea_row, endo])) {
        expect_equal(fs[[r_name]]$shea_partial_r2, fx[shea_row, endo],
                     tolerance = stata_tol$stat,
                     info = paste("FS Shea PR2", endo, cfg$suffix))
      }

      # Partial R²
      pr2_row <- which(fx_stats == "pr2")
      if (length(pr2_row) > 0 && !is.na(fx[pr2_row, endo])) {
        expect_equal(fs[[r_name]]$partial_r2, fx[pr2_row, endo],
                     tolerance = stata_tol$stat,
                     info = paste("FS PR2", endo, cfg$suffix))
      }

      # F-stat
      f_row <- which(fx_stats == "F")
      if (length(f_row) > 0 && !is.na(fx[f_row, endo])) {
        expect_equal(fs[[r_name]]$f_stat, fx[f_row, endo],
                     tolerance = stata_tol$stat,
                     info = paste("FS F", endo, cfg$suffix))
      }

      # SW chi2
      sw_chi2_row <- which(fx_stats == "SWchi2")
      if (length(sw_chi2_row) > 0 && !is.na(fx[sw_chi2_row, endo])) {
        expect_equal(fs[[r_name]]$sw_chi2, fx[sw_chi2_row, endo],
                     tolerance = stata_tol$stat,
                     info = paste("FS SW chi2", endo, cfg$suffix))
      }

      # SW F
      sw_f_row <- which(fx_stats == "SWF")
      if (length(sw_f_row) > 0 && !is.na(fx[sw_f_row, endo])) {
        expect_equal(fs[[r_name]]$sw_f, fx[sw_f_row, endo],
                     tolerance = stata_tol$stat,
                     info = paste("FS SW F", endo, cfg$suffix))
      }
    }
  })
}


# ============================================================================
# card_overid_dofminus: diagnostics (overidentified, 6 VCE combos)
# ============================================================================

for (cfg in vce_configs) {
  test_that(paste0("dofminus diagnostics match Stata: card_overid ", cfg$suffix), {
    skip_if(!file.exists(card_path), "card data not found")
    diag_path <- file.path(fixture_dir,
                           paste0("card_overid_dofminus_diagnostics_", cfg$suffix, ".csv"))
    skip_if(!file.exists(diag_path), "fixture not found")

    args <- list(formula = lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                 data = card, vcov = cfg$vcov, small = cfg$small,
                 dofminus = 1L, sdofminus = 1L)
    if (isTRUE(cfg$clusters)) args$clusters <- ~ smsa66
    fit <- do.call(ivreg2, args)
    d <- fit$diagnostics
    fx <- read.csv(diag_path)

    # Overidentification (Sargan for iid, Hansen J for robust/CL)
    if (!is.na(fx$sargan) && fx$sargan > 0) {
      expect_equal(d$overid$stat, fx$sargan, tolerance = stata_tol$stat,
                   info = paste("Sargan", cfg$suffix))
      expect_equal(d$overid$p, fx$sarganp, tolerance = stata_tol$pval,
                   info = paste("Sargan p", cfg$suffix))
    }
    if (!is.na(fx$j) && fx$j > 0) {
      expect_equal(d$overid$stat, fx$j, tolerance = stata_tol$stat,
                   info = paste("Hansen J", cfg$suffix))
      expect_equal(d$overid$p, fx$jp, tolerance = stata_tol$pval,
                   info = paste("Hansen J p", cfg$suffix))
    }

    # Underidentification
    if (!is.na(fx$idstat)) {
      expect_equal(d$underid$stat, fx$idstat, tolerance = stata_tol$stat,
                   info = paste("underid", cfg$suffix))
    }

    # Weak id
    if (!is.na(fx$cdf)) {
      expect_equal(d$weak_id$stat, fx$cdf, tolerance = stata_tol$stat,
                   info = paste("CD F", cfg$suffix))
    }
    if (!is.na(fx$widstat) && !is.null(d$weak_id_robust)) {
      expect_equal(d$weak_id_robust$stat, fx$widstat, tolerance = stata_tol$stat,
                   info = paste("KP F", cfg$suffix))
    }

    # Anderson-Rubin
    if (!is.na(fx$arf) && !is.na(d$anderson_rubin$f_stat)) {
      expect_equal(d$anderson_rubin$f_stat, fx$arf, tolerance = stata_tol$stat,
                   info = paste("AR F", cfg$suffix))
    }

    # Endogeneity
    if (!is.na(fx$estat) && fx$estat > 0) {
      expect_equal(d$endogeneity$stat, fx$estat, tolerance = stata_tol$stat,
                   info = paste("endogeneity", cfg$suffix))
    }

    # Model F
    if (!is.na(fx$F_stat)) {
      expect_equal(fit$model_f, fx$F_stat, tolerance = stata_tol$stat,
                   info = paste("model F", cfg$suffix))
    }
  })
}


# ============================================================================
# sim_cluster_dofminus: cluster-robust tests
# ============================================================================

for (cfg in vce_configs) {
  test_that(paste0("dofminus coefs match Stata: sim_cluster ", cfg$suffix), {
    skip_if(!file.exists(sim_path), "sim_cluster data not found")
    coef_path <- file.path(fixture_dir,
                           paste0("sim_cluster_dofminus_coef_", cfg$suffix, ".csv"))
    skip_if(!file.exists(coef_path), "fixture not found")

    args <- list(formula = y ~ x1 | endo1 | z1 + z2,
                 data = sim_cluster, vcov = cfg$vcov, small = cfg$small,
                 dofminus = 1L, sdofminus = 1L)
    if (isTRUE(cfg$clusters)) args$clusters <- ~ cluster_id
    fit <- do.call(ivreg2, args)
    fixture <- read.csv(coef_path)

    se_tol <- if (isTRUE(cfg$clusters)) 1e-5 else stata_tol$se
    for (i in seq_len(nrow(fixture))) {
      term <- fixture$term[i]
      r_name <- if (term == "_cons") "(Intercept)" else term
      expect_equal(unname(coef(fit)[r_name]), fixture$estimate[i],
                   tolerance = stata_tol$coef,
                   info = paste("coef", r_name, cfg$suffix))
      expect_equal(sqrt(fit$vcov[r_name, r_name]), fixture$std_error[i],
                   tolerance = se_tol,
                   info = paste("SE", r_name, cfg$suffix))
    }
  })
}

for (cfg in vce_configs) {
  test_that(paste0("dofminus diagnostics match Stata: sim_cluster ", cfg$suffix), {
    skip_if(!file.exists(sim_path), "sim_cluster data not found")
    diag_path <- file.path(fixture_dir,
                           paste0("sim_cluster_dofminus_diagnostics_", cfg$suffix, ".csv"))
    skip_if(!file.exists(diag_path), "fixture not found")

    args <- list(formula = y ~ x1 | endo1 | z1 + z2,
                 data = sim_cluster, vcov = cfg$vcov, small = cfg$small,
                 dofminus = 1L, sdofminus = 1L)
    if (isTRUE(cfg$clusters)) args$clusters <- ~ cluster_id
    fit <- do.call(ivreg2, args)
    d <- fit$diagnostics
    fx <- read.csv(diag_path)

    # Summary stats
    expect_equal(fit$sigma, sqrt(fx$sigmasq), tolerance = stata_tol$coef,
                 info = paste("sigma", cfg$suffix))

    # Overid
    if (!is.na(fx$sargan) && fx$sargan > 0) {
      expect_equal(d$overid$stat, fx$sargan, tolerance = stata_tol$stat,
                   info = paste("Sargan", cfg$suffix))
    }
    if (!is.na(fx$j) && fx$j > 0) {
      expect_equal(d$overid$stat, fx$j, tolerance = stata_tol$stat,
                   info = paste("Hansen J", cfg$suffix))
    }

    # Underid
    if (!is.na(fx$idstat)) {
      expect_equal(d$underid$stat, fx$idstat, tolerance = stata_tol$stat,
                   info = paste("underid", cfg$suffix))
    }

    # Weak id
    if (!is.na(fx$cdf)) {
      expect_equal(d$weak_id$stat, fx$cdf, tolerance = stata_tol$stat,
                   info = paste("CD F", cfg$suffix))
    }
    if (!is.na(fx$widstat) && !is.null(d$weak_id_robust)) {
      expect_equal(d$weak_id_robust$stat, fx$widstat, tolerance = stata_tol$stat,
                   info = paste("KP F", cfg$suffix))
    }

    # AR
    if (!is.na(fx$arf) && !is.na(d$anderson_rubin$f_stat)) {
      expect_equal(d$anderson_rubin$f_stat, fx$arf, tolerance = stata_tol$stat,
                   info = paste("AR F", cfg$suffix))
    }

    # Endogeneity
    if (!is.na(fx$estat) && fx$estat > 0) {
      expect_equal(d$endogeneity$stat, fx$estat, tolerance = stata_tol$stat,
                   info = paste("endogeneity", cfg$suffix))
    }

    # Model F
    if (!is.na(fx$F_stat)) {
      expect_equal(fit$model_f, fx$F_stat, tolerance = stata_tol$stat,
                   info = paste("model F", cfg$suffix))
    }
  })
}


# ============================================================================
# df.residual correctness
# ============================================================================

test_that("df.residual = N - K - dofminus - sdofminus (non-cluster)", {
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, dofminus = 1L, sdofminus = 1L)

  # N=3010, K=6 (5 regressors + intercept), dofminus=1, sdofminus=1
  expect_identical(fit$df.residual, 3010L - 6L - 1L - 1L)
})

test_that("df.residual = M - 1 for cluster models with dofminus", {
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, clusters = ~ smsa66,
                dofminus = 1L, sdofminus = 1L)

  # For clustered models, df.residual = M - 1 (not affected by dofminus)
  expect_identical(fit$df.residual, as.integer(fit$n_clusters - 1L))
})


# ============================================================================
# Object stores dofminus/sdofminus
# ============================================================================

test_that("dofminus and sdofminus are stored in the object", {
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, dofminus = 3L, sdofminus = 2L)

  expect_identical(fit$dofminus, 3L)
  expect_identical(fit$sdofminus, 2L)
})

test_that("default dofminus/sdofminus are 0L", {
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)

  expect_identical(fit$dofminus, 0L)
  expect_identical(fit$sdofminus, 0L)
})


# ============================================================================
# Input validation
# ============================================================================

test_that("dofminus rejects invalid inputs", {
  skip_if(!file.exists(card_path), "card data not found")

  # Negative values
  expect_error(ivreg2(lwage ~ exper, data = card, dofminus = -1L),
               "non-negative")
  expect_error(ivreg2(lwage ~ exper, data = card, sdofminus = -1L),
               "non-negative")
  # NA
  expect_error(ivreg2(lwage ~ exper, data = card, dofminus = NA_integer_),
               "non-negative")
  # Fractional values (should be rejected, not silently truncated)
  expect_error(ivreg2(lwage ~ exper, data = card, dofminus = 1.5),
               "non-negative integer")
  expect_error(ivreg2(lwage ~ exper, data = card, sdofminus = 2.9),
               "non-negative integer")
  # Non-numeric
  expect_error(ivreg2(lwage ~ exper, data = card, dofminus = "1"),
               "non-negative integer")
  # Inf and overflow doubles
  expect_error(ivreg2(lwage ~ exper, data = card, dofminus = Inf),
               "non-negative integer")
  expect_error(ivreg2(lwage ~ exper, data = card, sdofminus = Inf),
               "non-negative integer")
  expect_error(ivreg2(lwage ~ exper, data = card, dofminus = 1e20),
               "non-negative integer")
})

test_that("dofminus rejects values too large for model dimensions", {
  skip_if(!file.exists(card_path), "card data not found")

  # dofminus >= N
  expect_error(ivreg2(lwage ~ exper, data = card, dofminus = 5000L),
               "must be less than N")
  # N - K - dofminus - sdofminus <= 0
  expect_error(ivreg2(lwage ~ exper, data = card, dofminus = 3000L, sdofminus = 10L),
               "too large")
  # IV model: N - L - dofminus - sdofminus <= 0
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
           data = card, dofminus = 3000L, sdofminus = 5L),
    "too large"
  )
})
