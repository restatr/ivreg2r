# ============================================================================
# Tests: Partial Option (FWL Projection) — Ticket O1
# ============================================================================

# --- Helpers ---
card_path <- fixture_path("card_data.csv")

if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

mroz_path <- fixture_path("mroz_data.csv")
if (file.exists(mroz_path)) {
  mroz_fix <- read.csv(mroz_path)
}

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
  if (!is.null(diag$overid) && diag$overid$df > 0L &&
      !is.na(fixture$overid_stat) && fixture$overid_stat != "") {
    overid_val <- as.numeric(fixture$overid_stat)
    if (!is.na(overid_val)) {
      expect_equal(diag$overid$stat, overid_val,
                   tolerance = tol$stat, info = "overid stat")
      expect_equal(diag$overid$p, as.numeric(fixture$overid_p),
                   tolerance = tol$pval, info = "overid p")
    }
  }

  # Underid
  if (!is.null(diag$underid) && !is.na(fixture$underid_stat) &&
      fixture$underid_stat != "") {
    underid_val <- as.numeric(fixture$underid_stat)
    if (!is.na(underid_val)) {
      expect_equal(diag$underid$stat, underid_val,
                   tolerance = tol$stat, info = "underid stat")
      expect_equal(diag$underid$p, as.numeric(fixture$underid_p),
                   tolerance = tol$pval, info = "underid p")
    }
  }

  # Weak ID (Cragg-Donald)
  if (!is.null(diag$weak_id) && !is.na(fixture$weak_id_cd_f) &&
      fixture$weak_id_cd_f != "") {
    cd_val <- as.numeric(fixture$weak_id_cd_f)
    if (!is.na(cd_val)) {
      expect_equal(diag$weak_id$stat, cd_val,
                   tolerance = tol$stat, info = "CD F")
    }
  }

  # Weak ID (Kleibergen-Paap rk Wald F)
  if (!is.null(diag$weak_id_robust) && !is.na(fixture$weak_id_kp_f) &&
      fixture$weak_id_kp_f != "") {
    kp_val <- as.numeric(fixture$weak_id_kp_f)
    if (!is.na(kp_val)) {
      expect_equal(diag$weak_id_robust$stat, kp_val,
                   tolerance = tol$stat, info = "KP rk Wald F")
    }
  }

  # Model F
  if (!is.null(fit$model_f) && !is.na(fixture$model_f) &&
      fixture$model_f != "") {
    mf_val <- as.numeric(fixture$model_f)
    if (!is.na(mf_val)) {
      expect_equal(fit$model_f, mf_val,
                   tolerance = tol$stat, info = "model F")
    }
  }

  # Model F df
  if (!is.na(fixture$model_f_df1) && fixture$model_f_df1 != "") {
    expect_equal(fit$model_f_df1, as.integer(fixture$model_f_df1),
                 info = "model F df1")
    expect_equal(fit$model_f_df2, as.integer(fixture$model_f_df2),
                 info = "model F df2")
  }

  # Sigma
  if (!is.na(fixture$sigma) && fixture$sigma != "") {
    expect_equal(fit$sigma, as.numeric(fixture$sigma),
                 tolerance = tol$coef, info = "sigma")
  }

  # RSS
  if (!is.na(fixture$rss) && fixture$rss != "") {
    expect_equal(fit$rss, as.numeric(fixture$rss),
                 tolerance = tol$coef, info = "rss")
  }

  # N, K, sdofminus, partial_ct, partialcons
  expect_equal(fit$nobs, as.integer(fixture$N), info = "N")
  expect_equal(length(coef(fit)), as.integer(fixture$K), info = "K")
  expect_equal(fit$sdofminus, as.integer(fixture$sdofminus), info = "sdofminus")
  expect_equal(fit$partial_ct, as.integer(fixture$partial_ct), info = "partial_ct")
  expect_equal(as.integer(fit$partialcons), as.integer(fixture$partialcons),
               info = "partialcons")
}


# ===========================================================================
# 1. Basic partial: partial(black south smsa) — IV, IID
# ===========================================================================
test_that("partial(black south smsa) matches Stata — IID", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"))

  # Coefficients and SEs
  check_coef_fixture(fit,
    fixture_path("card_partial_basic_coef_iid.csv"))

  # VCV
  check_vcov_fixture(fit,
    fixture_path("card_partial_basic_vcov_iid.csv"),
    fixture_path("card_partial_basic_coef_iid.csv"))

  # Diagnostics
  check_diag_fixture(fit,
    fixture_path("card_partial_basic_diagnostics_iid.csv"))

  # Metadata
  expect_equal(fit$partial_ct, 4L)
  expect_true(fit$partialcons)
  expect_equal(sort(fit$partial_names), sort(c("black", "south", "smsa")))
  expect_equal(fit$sdofminus, 4L)

  # Only non-partialled coefficients reported
  expect_equal(length(coef(fit)), 3L)
  expect_true(all(c("educ", "exper", "expersq") %in% names(coef(fit))))
  expect_false("(Intercept)" %in% names(coef(fit)))
  expect_false("black" %in% names(coef(fit)))
})


# ===========================================================================
# 1b. Basic partial — IID small
# ===========================================================================
test_that("partial(black south smsa) matches Stata — IID small", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"),
                small = TRUE)

  check_coef_fixture(fit,
    fixture_path("card_partial_basic_coef_iid_small.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_basic_diagnostics_iid_small.csv"))
})


# ===========================================================================
# 2. partial(_cons) — demean only
# ===========================================================================
test_that("partial(_cons) matches Stata — IID", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(
    lwage ~ exper + expersq + black + south + smsa | educ | nearc2 + nearc4,
    data = card, partial = "_cons")

  check_coef_fixture(fit,
    fixture_path("card_partial_cons_coef_iid.csv"))
  check_vcov_fixture(fit,
    fixture_path("card_partial_cons_vcov_iid.csv"),
    fixture_path("card_partial_cons_coef_iid.csv"))

  diag_fix <- read.csv(
    fixture_path("card_partial_cons_diagnostics_iid.csv"))
  expect_equal(fit$partial_ct, as.integer(diag_fix$partial_ct))
  expect_true(fit$partialcons)

  # 6 coefficients reported (all exog + endo, no intercept)
  expect_equal(length(coef(fit)), 6L)
  expect_false("(Intercept)" %in% names(coef(fit)))
})

# ===========================================================================
# 2b. (Intercept) accepted as synonym for _cons
# ===========================================================================
test_that("partial accepts both '_cons' and '(Intercept)'", {
  skip_if_not(file.exists(card_path))
  fit1 <- ivreg2(
    lwage ~ exper + expersq + black + south + smsa | educ | nearc2 + nearc4,
    data = card, partial = "_cons")
  fit2 <- ivreg2(
    lwage ~ exper + expersq + black + south + smsa | educ | nearc2 + nearc4,
    data = card, partial = "(Intercept)")
  expect_equal(coef(fit1), coef(fit2))
  expect_equal(vcov(fit1), vcov(fit2))
})


# ===========================================================================
# 3. partial(_all) — partial all exogenous
# ===========================================================================
test_that("partial(_all) matches Stata — IID", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(
    lwage ~ exper + expersq + black + south + smsa |
      educ | nearc2 + nearc4,
    data = card, partial = "_all")

  check_coef_fixture(fit,
    fixture_path("card_partial_all_coef_iid.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_all_diagnostics_iid.csv"))

  # Only endogenous variable(s) remain
  expect_equal(length(coef(fit)), 1L)
  expect_equal(names(coef(fit)), "educ")
})

test_that("partial(_all) id tests equal partial(_cons) id tests (FWL invariance)", {
  skip_if_not(file.exists(card_path))
  fit_all <- ivreg2(
    lwage ~ exper + expersq + black + south + smsa |
      educ | nearc2 + nearc4,
    data = card, partial = "_all")
  fit_cons <- ivreg2(
    lwage ~ exper + expersq + black + south + smsa |
      educ | nearc2 + nearc4,
    data = card, partial = "_cons")

  expect_equal(fit_all$diagnostics$underid$stat,
               fit_cons$diagnostics$underid$stat)
  expect_equal(fit_all$diagnostics$weak_id$stat,
               fit_cons$diagnostics$weak_id$stat)
})


# ===========================================================================
# 3b. F5 regression: id tests with an empty exogenous block
#     (partial(_all) with K1 > 0, and noconstant with K2 = 0)
# ===========================================================================
test_that("partial(_all) id tests match Stata — mroz, IID", {
  skip_if_not(file.exists(mroz_path))
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_fix, partial = "_all")

  check_coef_fixture(fit,
    fixture_path("mroz_partial_all_coef_iid.csv"))
  check_diag_fixture(fit,
    fixture_path("mroz_partial_all_diagnostics_iid.csv"))
})

test_that("partial(_all) id tests match Stata — mroz, robust (KP rk)", {
  skip_if_not(file.exists(mroz_path))
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_fix, partial = "_all", vcov = "robust")

  check_coef_fixture(fit,
    fixture_path("mroz_partial_all_coef_robust.csv"))
  check_diag_fixture(fit,
    fixture_path("mroz_partial_all_diagnostics_robust.csv"))
})

test_that("noconstant with no exogenous regressors id tests match Stata — mroz, IID", {
  skip_if_not(file.exists(mroz_path))
  fit <- ivreg2(lwage ~ 0 | educ | age + kidslt6 + kidsge6, data = mroz_fix)

  check_coef_fixture(fit,
    fixture_path("mroz_nocons_coef_iid.csv"))
  check_diag_fixture(fit,
    fixture_path("mroz_nocons_diagnostics_iid.csv"))
})

test_that("weighted partial(_all) id tests match Stata — mroz, IID", {
  skip_if_not(file.exists(mroz_path))
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_fix, weights = hours + 1, partial = "_all")

  check_coef_fixture(fit,
    fixture_path("mroz_partial_all_w_coef_iid.csv"))
  check_diag_fixture(fit,
    fixture_path("mroz_partial_all_w_diagnostics_iid.csv"))
})

test_that("weighted partial(_all) id tests match Stata — mroz, robust (KP rk)", {
  skip_if_not(file.exists(mroz_path))
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_fix, weights = hours + 1, partial = "_all",
                vcov = "robust")

  check_coef_fixture(fit,
    fixture_path("mroz_partial_all_w_coef_robust.csv"))
  check_diag_fixture(fit,
    fixture_path("mroz_partial_all_w_diagnostics_robust.csv"))
})


# ===========================================================================
# 4. Weighted partial
# ===========================================================================
test_that("weighted partial matches Stata — IID", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, weights = weight,
                partial = c("black", "south", "smsa"))

  check_coef_fixture(fit,
    fixture_path("card_partial_weighted_coef_iid.csv"))
  check_vcov_fixture(fit,
    fixture_path("card_partial_weighted_vcov_iid.csv"),
    fixture_path("card_partial_weighted_coef_iid.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_weighted_diagnostics_iid.csv"))
})


# ===========================================================================
# 5. nopartialsmall
# ===========================================================================
test_that("nopartialsmall suppresses sdofminus increment — IID", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"),
                nopartialsmall = TRUE)

  expect_equal(fit$sdofminus, 0L)
  expect_equal(fit$partial_ct, 4L)  # count unchanged; only sdofminus suppressed

  check_coef_fixture(fit,
    fixture_path("card_partial_nosmall_coef_iid.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_nosmall_diagnostics_iid.csv"))
})

test_that("nopartialsmall matches Stata — IID small", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"),
                nopartialsmall = TRUE, small = TRUE)

  check_coef_fixture(fit,
    fixture_path("card_partial_nosmall_coef_iid_small.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_nosmall_diagnostics_iid_small.csv"))
})


# ===========================================================================
# 6. OLS with partial
# ===========================================================================
test_that("OLS partial matches Stata — IID", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa,
                data = card, partial = c("black", "south"))

  check_coef_fixture(fit,
    fixture_path("card_partial_ols_coef_iid.csv"))
  check_vcov_fixture(fit,
    fixture_path("card_partial_ols_vcov_iid.csv"),
    fixture_path("card_partial_ols_coef_iid.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_ols_diagnostics_iid.csv"))
})


# ===========================================================================
# 7. LIML + partial
# ===========================================================================
test_that("LIML + partial matches Stata — IID", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"),
                method = "liml")

  check_coef_fixture(fit,
    fixture_path("card_partial_liml_coef_iid.csv"))
  check_vcov_fixture(fit,
    fixture_path("card_partial_liml_vcov_iid.csv"),
    fixture_path("card_partial_liml_coef_iid.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_liml_diagnostics_iid.csv"))
})


# ===========================================================================
# 8. GMM2S + partial
# ===========================================================================
test_that("GMM2S + partial matches Stata — robust", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"),
                method = "gmm2s", vcov = "robust")

  check_coef_fixture(fit,
    fixture_path("card_partial_gmm2s_coef_robust.csv"))
  check_vcov_fixture(fit,
    fixture_path("card_partial_gmm2s_vcov_robust.csv"),
    fixture_path("card_partial_gmm2s_coef_robust.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_gmm2s_diagnostics_robust.csv"))
})


# ===========================================================================
# 9. Robust + partial
# ===========================================================================
test_that("robust partial matches Stata — HC0", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"),
                vcov = "robust")

  check_coef_fixture(fit,
    fixture_path("card_partial_robust_coef_hc0.csv"))
  check_vcov_fixture(fit,
    fixture_path("card_partial_robust_vcov_hc0.csv"),
    fixture_path("card_partial_robust_coef_hc0.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_robust_diagnostics_hc0.csv"))
})

test_that("robust small partial matches Stata — HC1 small", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"),
                vcov = "robust", small = TRUE)

  check_coef_fixture(fit,
    fixture_path("card_partial_robust_coef_hc1_small.csv"))
  check_vcov_fixture(fit,
    fixture_path("card_partial_robust_vcov_hc1_small.csv"),
    fixture_path("card_partial_robust_coef_hc1_small.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_robust_diagnostics_hc1_small.csv"))
})


# ===========================================================================
# 10. Cluster + partial
# ===========================================================================
test_that("cluster partial matches Stata — CL", {
  skip_if_not(file.exists(card_path))
  fit <- muffle_rank_warnings(ivreg2(
    lwage ~ exper + expersq + black + south + smsa | educ | nearc2 + nearc4,
    data = card, partial = c("black", "south", "smsa"),
    clusters = ~smsa66))

  check_coef_fixture(fit,
    fixture_path("card_partial_cluster_coef_cl.csv"))
  check_vcov_fixture(fit,
    fixture_path("card_partial_cluster_vcov_cl.csv"),
    fixture_path("card_partial_cluster_coef_cl.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_cluster_diagnostics_cl.csv"))
})

test_that("cluster small partial matches Stata — CL small", {
  skip_if_not(file.exists(card_path))
  fit <- muffle_rank_warnings(ivreg2(
    lwage ~ exper + expersq + black + south + smsa | educ | nearc2 + nearc4,
    data = card, partial = c("black", "south", "smsa"),
    clusters = ~smsa66, small = TRUE))

  check_coef_fixture(fit,
    fixture_path("card_partial_cluster_coef_cl_small.csv"))
  check_diag_fixture(fit,
    fixture_path("card_partial_cluster_diagnostics_cl_small.csv"))
})


# ===========================================================================
# Predict restrictions
# ===========================================================================
test_that("predict() returns fitted values with message for partial models", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"))
  expect_message(predict(fit), "non-partialled regressors")
})

test_that("predict(newdata) errors for partial models", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"))
  expect_error(predict(fit, newdata = card[1:5, ]),
               "Cannot predict on new data after partialling")
})


# ===========================================================================
# Metadata tests
# ===========================================================================
test_that("partial_ct and partial_names are correct", {
  skip_if_not(file.exists(card_path))

  # partial(black south smsa) → partial_ct = 4 (3 vars + cons)
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                   educ | nearc2 + nearc4,
                 data = card, partial = c("black", "south", "smsa"))
  expect_equal(fit1$partial_ct, 4L)
  expect_true(fit1$partialcons)
  expect_equal(sort(fit1$partial_names), sort(c("black", "south", "smsa")))

  # partial(_cons) → partial_ct = 1
  fit2 <- ivreg2(
    lwage ~ exper + expersq + black + south + smsa | educ | nearc2 + nearc4,
    data = card, partial = "_cons")
  expect_equal(fit2$partial_ct, 1L)
  expect_true(fit2$partialcons)
  expect_equal(fit2$partial_names, character(0L))

  # partial(_all) → partial_ct = 6 (5 exog vars + cons)
  fit3 <- ivreg2(
    lwage ~ exper + expersq + black + south + smsa |
      educ | nearc2 + nearc4,
    data = card, partial = "_all")
  expect_equal(fit3$partial_ct, 6L)
  expect_true(fit3$partialcons)
})

test_that("no partial → partial_ct = 0", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4, data = card)
  expect_equal(fit$partial_ct, 0L)
  expect_false(fit$partialcons)
  expect_equal(fit$partial_names, character(0L))
})


# ===========================================================================
# Error tests
# ===========================================================================
test_that("partial with invalid variable names errors", {
  skip_if_not(file.exists(card_path))
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south + smsa |
             educ | nearc2 + nearc4,
           data = card, partial = c("black", "nonexistent")),
    "not in the exogenous regressor list"
  )
})

test_that("partial _cons with noconstant model errors", {
  skip_if_not(file.exists(card_path))
  expect_error(
    ivreg2(lwage ~ 0 + exper + expersq + black + south + smsa |
             educ | nearc2 + nearc4,
           data = card, partial = "_cons"),
    "without an intercept"
  )
})

test_that("partial must be character or NULL", {
  skip_if_not(file.exists(card_path))
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south + smsa |
             educ | nearc2 + nearc4,
           data = card, partial = 42),
    "character vector or NULL"
  )
})

test_that("nopartialsmall must be logical", {
  skip_if_not(file.exists(card_path))
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south + smsa |
             educ | nearc2 + nearc4,
           data = card, partial = c("black"), nopartialsmall = "yes"),
    "TRUE or FALSE"
  )
})

test_that("partial cannot include endogenous variables", {
  skip_if_not(file.exists(card_path))
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south + smsa |
             educ | nearc2 + nearc4,
           data = card, partial = c("educ")),
    "not in the exogenous regressor list"
  )
})


# ===========================================================================
# Summary footer
# ===========================================================================
test_that("summary prints partial footer", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"))
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Partialled out", out)))
  expect_true(any(grepl("_cons", out)))
  expect_true(any(grepl("partial-model", out)))
})


# ===========================================================================
# Broom methods with partial
# ===========================================================================
test_that("glance includes partial_ct", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"))
  gl <- glance(fit)
  expect_equal(gl$partial_ct, 4L)
})

test_that("tidy works with partial model", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"))
  td <- tidy(fit)
  expect_equal(nrow(td), 3L)
  expect_true(all(c("educ", "exper", "expersq") %in% td$term))
})


# ===========================================================================
# FWL coefficient invariance (pure R test)
# ===========================================================================
test_that("FWL theorem: partial coefficients match full model coefficients", {
  skip_if_not(file.exists(card_path))
  # Run full model
  fit_full <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                       educ | nearc2 + nearc4, data = card)
  # Run partial model
  fit_partial <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                          educ | nearc2 + nearc4, data = card,
                        partial = c("black", "south", "smsa"))
  # Coefficients on non-partialled variables should be identical
  shared <- intersect(names(coef(fit_full)), names(coef(fit_partial)))
  for (nm in shared) {
    expect_equal(unname(coef(fit_partial)[nm]), unname(coef(fit_full)[nm]),
                 tolerance = 1e-10,
                 info = paste("FWL invariance:", nm))
  }
})

test_that("FWL theorem: partial(_cons) coefficients match full model", {
  skip_if_not(file.exists(card_path))
  fit_full <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                       educ | nearc2 + nearc4, data = card)
  fit_cons <- ivreg2(
    lwage ~ exper + expersq + black + south + smsa | educ | nearc2 + nearc4,
    data = card, partial = "_cons")
  shared <- intersect(names(coef(fit_full)), names(coef(fit_cons)))
  for (nm in shared) {
    expect_equal(unname(coef(fit_cons)[nm]), unname(coef(fit_full)[nm]),
                 tolerance = 1e-10,
                 info = paste("FWL cons invariance:", nm))
  }
})


# ===========================================================================
# Interaction with sdofminus parameter
# ===========================================================================
# ===========================================================================
# CUE + partial: warning and non-invariance
# ===========================================================================
test_that("CUE + partial emits FWL non-invariance warning", {
  skip_if_not(file.exists(card_path))
  expect_warning(
    ivreg2(lwage ~ exper + expersq + black + south + smsa |
             educ | nearc2 + nearc4,
           data = card, method = "cue", partial = "smsa"),
    "FWL invariance does not hold for CUE"
  )
})

test_that("CUE + partial is NOT FWL-invariant with robust VCE", {
  skip_if_not(file.exists(card_path))
  # Under IID, CUE = LIML which IS FWL-invariant.
  # Non-invariance requires a heteroskedasticity-robust S matrix.
  fit_full <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                       educ | nearc2 + nearc4,
                     data = card, method = "cue", vcov = "robust")
  expect_warning(
    fit_partial <- ivreg2(
      lwage ~ exper + expersq + black + south + smsa |
        educ | nearc2 + nearc4,
      data = card, method = "cue", vcov = "robust", partial = "smsa"),
    "FWL invariance does not hold for CUE"
  )
  shared <- intersect(names(coef(fit_full)), names(coef(fit_partial)))
  # At least one coefficient should differ (non-invariance)
  diffs <- vapply(shared, function(nm) {
    abs(coef(fit_full)[nm] - coef(fit_partial)[nm])
  }, numeric(1))
  expect_true(max(diffs) > 1e-6,
              info = "CUE + partial should NOT be FWL-invariant with robust VCE")
})


test_that("CUE b0 + partial works (b0 validated after partialling)", {
  skip_if_not(file.exists(card_path))
  # Regression test: b0 must be validated AFTER partialling removes columns.
  # Before the fix, b0 was validated against pre-partial K, causing a
  # dimension mismatch when partial removes exogenous regressors.
  # After partialling "smsa", remaining regressors are:
  #   educ, exper, expersq, black, south (K=5, no intercept)
  # CUE + IID = LIML, so we can compare against LIML + partial.
  fit_liml <- ivreg2(
    lwage ~ exper + expersq + black + south + smsa | educ | nearc2 + nearc4,
    data = card, method = "liml", partial = "smsa")
  b0_vec <- coef(fit_liml)
  fit_cue <- ivreg2(
    lwage ~ exper + expersq + black + south + smsa | educ | nearc2 + nearc4,
    data = card, method = "cue", partial = "smsa", b0 = b0_vec)
  # CUE + IID should match LIML
  for (nm in names(b0_vec)) {
    expect_equal(unname(coef(fit_cue)[nm]), unname(coef(fit_liml)[nm]),
                 tolerance = 1e-4,
                 info = paste("b0+partial CUE vs LIML:", nm))
  }
})

test_that("partial + sdofminus stack correctly", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa |
                  educ | nearc2 + nearc4,
                data = card, partial = c("black", "south", "smsa"),
                sdofminus = 2L)
  # sdofminus = 2 (user) + 4 (partial_ct) = 6
  expect_equal(fit$sdofminus, 6L)
})
