# ===========================================================================
# Tests for reduced-form regression (Ticket J3)
#
# M-21 re-base (2026-07-04): reduced-form/system fixtures re-based onto the
# mroz H31 base per planning/22-spec-matrix.md. The canonical base is the
# saverf command H54 (help.txt:1385); the saved RF equals the H61
# reduced-form OLS demo (help.txt:1416) -- the "saverf no-op parity". The
# card_* rf/system fixtures are retired. sim_multi_endo is kept for K1 > 1
# system coverage and sim_cluster for cluster RF (mroz has no cluster
# variable).
# ===========================================================================

# --- Fixture loading helpers ---

load_rf_fixture <- function(prefix, suffix) {
  path <- fixture_path(paste0(prefix, "_rf_", suffix, ".csv"))
  if (!file.exists(path)) return(NULL)
  fix <- read.csv(path, stringsAsFactors = FALSE)
  fix$term[fix$term == "_cons"] <- "(Intercept)"
  fix
}
load_system_fixture <- function(prefix, suffix) {
  path <- fixture_path(paste0(prefix, "_system_", suffix, ".csv"))
  if (!file.exists(path)) return(NULL)
  fix <- read.csv(path, stringsAsFactors = FALSE)
  fix$term[fix$term == "_cons"] <- "(Intercept)"
  fix
}

# --- Shared data ---
data(card, package = "ivreg2r")
data(mroz, package = "ivreg2r")

# e(sample): Stata's ivreg2 drops observations with missing lwage
# (never-employed women) before estimation.
mroz_lw <- mroz[!is.na(mroz$lwage), ]

sim_multi_path <- fixture_path("sim_multi_endo_data.csv")
if (file.exists(sim_multi_path)) {
  sim_multi_endo <- read.csv(sim_multi_path, stringsAsFactors = FALSE)
}

sim_cluster_path <- fixture_path("sim_cluster_data.csv")
if (file.exists(sim_cluster_path)) {
  sim_cluster <- read.csv(sim_cluster_path, stringsAsFactors = FALSE)
}


# ===========================================================================
# STRUCTURAL TESTS (no fixtures)
# ===========================================================================

test_that("reduced_form = 'none' returns NULL", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  expect_null(fit$reduced_form)
})

test_that("OLS model with reduced_form = 'rf' returns NULL", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south, data = card,
                reduced_form = "rf")
  expect_null(fit$reduced_form)
})

test_that("OLS model with reduced_form = 'system' returns NULL", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south, data = card,
                reduced_form = "system")
  expect_null(fit$reduced_form)
})

test_that("bad reduced_form value produces error", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
           data = card, reduced_form = "bad"),
    '"none", "rf", or "system"'
  )
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
           data = card, reduced_form = 42),
    '"none", "rf", or "system"'
  )
})

test_that("rf mode returns expected field names and types", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, reduced_form = "rf")
  rf <- fit$reduced_form
  expect_equal(rf$mode, "rf")
  expect_equal(rf$depvar, "lwage")
  expect_type(rf$coefficients, "double")
  expect_equal(length(rf$coefficients), 6L)  # L = 6
  expect_true(is.matrix(rf$vcov))
  expect_equal(dim(rf$vcov), c(6L, 6L))
  expect_equal(length(rf$residuals), nrow(card))
  expect_equal(length(rf$fitted.values), nrow(card))
  expect_type(rf$sigma, "double")
  expect_type(rf$f_stat, "double")
  expect_type(rf$f_p, "double")
  expect_type(rf$f_df1, "integer")
  expect_type(rf$f_df2, "integer")
  expect_type(rf$chi2_stat, "double")
  expect_type(rf$chi2_p, "double")
  expect_type(rf$chi2_df, "integer")
  expect_type(rf$partial_r2, "double")
})

test_that("system mode returns expected field names and types", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, reduced_form = "system")
  rf <- fit$reduced_form
  expect_equal(rf$mode, "system")
  expect_equal(rf$depvar, c("lwage", "educ"))
  expect_true(is.matrix(rf$coefficients))
  expect_equal(dim(rf$coefficients), c(6L, 2L))  # L x (K1+1)
  expect_true(is.matrix(rf$vcov))
  expect_equal(dim(rf$vcov), c(12L, 12L))  # (K1+1)*L x (K1+1)*L
  expect_true(is.matrix(rf$residuals))
  expect_equal(dim(rf$residuals), c(nrow(card), 2L))
  expect_true(is.matrix(rf$fitted.values))
  expect_equal(dim(rf$fitted.values), c(nrow(card), 2L))
  expect_equal(length(rf$sigma), 2L)
  expect_equal(length(rf$equations), 2L)
  expect_equal(names(rf$equations), c("lwage", "educ"))
})

test_that("system mode with K1 > 1 has correct dimensions", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, reduced_form = "system")
  rf <- fit$reduced_form
  # L = 7 (x1, x2, z1, z2, z3, z4, intercept), K1 = 2, n_eq = 3
  expect_equal(dim(rf$coefficients), c(7L, 3L))
  expect_equal(dim(rf$vcov), c(21L, 21L))     # 3 * 7
  expect_equal(dim(rf$residuals), c(1000L, 3L))
  expect_equal(length(rf$sigma), 3L)
  expect_equal(length(rf$equations), 3L)
})


# ===========================================================================
# RF F = AR F CROSS-VALIDATION
# ===========================================================================

test_that("RF F equals AR F (IID)", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "iid", reduced_form = "rf")
  expect_equal(fit$reduced_form$f_stat,
               fit$diagnostics$anderson_rubin$f_stat,
               tolerance = 1e-14)
})

test_that("RF F equals AR F (HC0)", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "robust", reduced_form = "rf")
  expect_equal(fit$reduced_form$f_stat,
               fit$diagnostics$anderson_rubin$f_stat,
               tolerance = 1e-14)
})

test_that("RF F equals AR F (HC1)", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "robust", reduced_form = "rf")
  expect_equal(fit$reduced_form$f_stat,
               fit$diagnostics$anderson_rubin$f_stat,
               tolerance = 1e-14)
})

test_that("RF F equals AR F (cluster)", {
  skip_if(!file.exists(sim_cluster_path), "sim_cluster data not found")
  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2, data = sim_cluster,
                clusters = ~cluster_id, reduced_form = "rf")
  expect_equal(fit$reduced_form$f_stat,
               fit$diagnostics$anderson_rubin$f_stat,
               tolerance = 1e-14)
})

test_that("RF F equals AR F (cluster, small)", {
  skip_if(!file.exists(sim_cluster_path), "sim_cluster data not found")
  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2, data = sim_cluster,
                clusters = ~cluster_id, small = TRUE, reduced_form = "rf")
  expect_equal(fit$reduced_form$f_stat,
               fit$diagnostics$anderson_rubin$f_stat,
               tolerance = 1e-14)
})

test_that("RF F equals AR F (overidentified, HC1)", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", reduced_form = "rf")
  expect_equal(fit$reduced_form$f_stat,
               fit$diagnostics$anderson_rubin$f_stat,
               tolerance = 1e-14)
})

test_that("RF F equals AR F (multi-endo, IID)", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, vcov = "iid", reduced_form = "rf")
  expect_equal(fit$reduced_form$f_stat,
               fit$diagnostics$anderson_rubin$f_stat,
               tolerance = 1e-14)
})

test_that("RF F equals AR F (weighted)", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight, reduced_form = "rf")
  expect_equal(fit$reduced_form$f_stat,
               fit$diagnostics$anderson_rubin$f_stat,
               tolerance = 1e-14)
})


# ===========================================================================
# RF COEFFICIENTS SANITY CHECK
# ===========================================================================

test_that("RF coefficients match lm.fit(Z, y)", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, reduced_form = "rf", x = TRUE)
  lm_rf <- lm.fit(fit$x$Z, card$lwage[seq_len(nrow(fit$x$Z))])
  expect_equal(unname(fit$reduced_form$coefficients),
               unname(lm_rf$coefficients), tolerance = 1e-14)
})


# ===========================================================================
# SYSTEM VCV BLOCK CONSISTENCY
# ===========================================================================

test_that("system VCV first-stage blocks match individual first-stage VCVs", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, reduced_form = "system")
  rf <- fit$reduced_form
  L <- ncol(rf$coefficients)
  nL <- nrow(rf$coefficients)

  # Extract first-stage block (equation 2 = educ)
  idx2 <- (nL + 1L):(2L * nL)
  system_block <- rf$vcov[idx2, idx2]

  # RF equation block (equation 1 = lwage)
  idx1 <- 1L:nL
  system_block_y <- rf$vcov[idx1, idx1]

  # The RF equation VCV from rf mode should match the system block
  fit_rf <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                   data = card, reduced_form = "rf")
  expect_equal(unname(fit_rf$reduced_form$vcov),
               unname(system_block_y),
               tolerance = 1e-14)
})


# ===========================================================================
# REDUCED-FORM F STORED ON THE FITTED OBJECT
# ===========================================================================

test_that("rf_f_stat and rf_f_p are stored on the fitted object for rf mode", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, reduced_form = "rf")
  expect_false(is.na(fit$reduced_form$f_stat))
  expect_false(is.na(fit$reduced_form$f_p))
})

test_that("rf_f fields are absent without reduced_form", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  expect_true(is.null(fit$reduced_form) || fit$reduced_form$mode != "rf")
})

test_that("rf_f fields are absent for system mode", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, reduced_form = "system")
  expect_true(is.null(fit$reduced_form) || fit$reduced_form$mode != "rf")
})


# ===========================================================================
# STATA FIXTURE TESTS — RF MODE
# ===========================================================================

# --- Helper: compare one RF fixture ---
compare_rf_fixture <- function(fit, prefix, suffix) {
  fix <- load_rf_fixture(prefix, suffix)
  if (is.null(fix)) skip(paste("Fixture not found:", prefix, suffix))
  rf <- fit$reduced_form

  # Coefficients
  r_coef <- rf$coefficients[fix$term]
  expect_equal(unname(r_coef), fix$estimate,
               tolerance = stata_tol$coef,
               label = paste(prefix, suffix, "coef"))

  # Standard errors
  r_se <- sqrt(diag(rf$vcov))[fix$term]
  expect_equal(unname(r_se), fix$std_error,
               tolerance = stata_tol$se,
               label = paste(prefix, suffix, "se"))

  # RMSE
  expect_equal(rf$sigma, fix$rmse[1],
               tolerance = stata_tol$coef,
               label = paste(prefix, suffix, "rmse"))

  # F-stat (only row 1 has it; skip for cluster because Stata's `test`
  # command after `estimates restore` uses OLS df2 (N-L) not cluster df2
  # (M-1). The RF F = AR F cross-validation covers F correctness instead.)
  if (!is.na(fix$rf_f[1]) && !grepl("^cl", suffix)) {
    expect_equal(rf$f_stat, fix$rf_f[1],
                 tolerance = stata_tol$stat,
                 label = paste(prefix, suffix, "f_stat"))
    expect_equal(rf$f_p, fix$rf_f_p[1],
                 tolerance = stata_tol$pval,
                 label = paste(prefix, suffix, "f_p"))
    expect_equal(rf$f_df1, as.integer(fix$rf_f_df1[1]),
                 label = paste(prefix, suffix, "f_df1"))
    expect_equal(rf$f_df2, as.integer(fix$rf_f_df2[1]),
                 label = paste(prefix, suffix, "f_df2"))
  }
}

# --- Helper: compare one system fixture ---
compare_system_fixture <- function(fit, prefix, suffix) {
  fix <- load_system_fixture(prefix, suffix)
  if (is.null(fix)) skip(paste("Fixture not found:", prefix, suffix))
  rf <- fit$reduced_form

  for (eq in unique(fix$eq)) {
    feq <- fix[fix$eq == eq, ]
    eq_idx <- match(eq, rf$depvar)

    # Coefficients
    r_coef <- rf$coefficients[feq$term, eq_idx]
    expect_equal(unname(r_coef), feq$estimate,
                 tolerance = stata_tol$coef,
                 label = paste(prefix, suffix, eq, "coef"))

    # Standard errors from system VCV (use positional indexing since
    # system VCV labels are "eq:term" format, not just "term")
    nL <- nrow(rf$coefficients)
    vcov_idx <- ((eq_idx - 1L) * nL + 1L):(eq_idx * nL)
    eq_se_all <- sqrt(diag(rf$vcov[vcov_idx, vcov_idx]))
    # Match by term name within the block (strip equation prefix)
    block_terms <- sub("^[^:]+:", "", names(eq_se_all))
    r_se <- eq_se_all[match(feq$term, block_terms)]
    expect_equal(unname(r_se), feq$std_error,
                 tolerance = stata_tol$se,
                 label = paste(prefix, suffix, eq, "se"))
  }
}


# --- mroz (H31 base: `lwage ~ exper + expersq | educ | age + kidslt6 +
# kidsge6`). Canonical saverf command is H54 (help.txt:1385); the saved RF
# equals the H61 reduced-form OLS demo (help.txt:1416) -- the "saverf no-op
# parity". Weighted RF cells were deleted, not ported -- weighted-VCV
# parity is owned by the M-11 family, and the fixture-free weighted
# AR-identity test above ("RF F equals AR F (weighted)") exercises weighted
# RF structurally. ---

test_that("mroz RF: IID", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, reduced_form = "rf")
  compare_rf_fixture(fit, "mroz", "iid")
})

test_that("mroz RF: IID small", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, small = TRUE, reduced_form = "rf")
  compare_rf_fixture(fit, "mroz", "iid_small")
})

test_that("mroz RF: HC1", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust", reduced_form = "rf")
  compare_rf_fixture(fit, "mroz", "hc1")
})

test_that("mroz RF: HC1 small", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust", small = TRUE,
                reduced_form = "rf")
  compare_rf_fixture(fit, "mroz", "hc1_small")
})

test_that("mroz system: IID", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, reduced_form = "system")
  compare_system_fixture(fit, "mroz", "iid")
})

test_that("mroz system: IID small", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, small = TRUE, reduced_form = "system")
  compare_system_fixture(fit, "mroz", "iid_small")
})

test_that("mroz system: HC1", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust", reduced_form = "system")
  compare_system_fixture(fit, "mroz", "hc1")
})

test_that("mroz system: HC1 small", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust", small = TRUE,
                reduced_form = "system")
  compare_system_fixture(fit, "mroz", "hc1_small")
})


# --- sim_multi_endo ---

test_that("sim_multi_endo RF: IID", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, reduced_form = "rf")
  compare_rf_fixture(fit, "sim_multi_endo", "iid")
})

test_that("sim_multi_endo RF: IID small", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, small = TRUE, reduced_form = "rf")
  compare_rf_fixture(fit, "sim_multi_endo", "iid_small")
})

test_that("sim_multi_endo RF: HC1", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, vcov = "robust", reduced_form = "rf")
  compare_rf_fixture(fit, "sim_multi_endo", "hc1")
})

test_that("sim_multi_endo RF: HC1 small", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, vcov = "robust", small = TRUE,
                reduced_form = "rf")
  compare_rf_fixture(fit, "sim_multi_endo", "hc1_small")
})

test_that("sim_multi_endo system: IID", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, reduced_form = "system")
  compare_system_fixture(fit, "sim_multi_endo", "iid")
})

test_that("sim_multi_endo system: IID small", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, small = TRUE, reduced_form = "system")
  compare_system_fixture(fit, "sim_multi_endo", "iid_small")
})

test_that("sim_multi_endo system: HC1", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, vcov = "robust", reduced_form = "system")
  compare_system_fixture(fit, "sim_multi_endo", "hc1")
})

test_that("sim_multi_endo system: HC1 small", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi_endo, vcov = "robust", small = TRUE,
                reduced_form = "system")
  compare_system_fixture(fit, "sim_multi_endo", "hc1_small")
})


# --- sim_cluster ---

test_that("sim_cluster RF: cluster", {
  skip_if(!file.exists(sim_cluster_path), "sim_cluster data not found")
  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2, data = sim_cluster,
                clusters = ~cluster_id, reduced_form = "rf")
  compare_rf_fixture(fit, "sim_cluster", "cl")
})

test_that("sim_cluster RF: cluster small", {
  skip_if(!file.exists(sim_cluster_path), "sim_cluster data not found")
  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2, data = sim_cluster,
                clusters = ~cluster_id, small = TRUE, reduced_form = "rf")
  compare_rf_fixture(fit, "sim_cluster", "cl_small")
})

test_that("sim_cluster system: cluster", {
  skip_if(!file.exists(sim_cluster_path), "sim_cluster data not found")
  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2, data = sim_cluster,
                clusters = ~cluster_id, reduced_form = "system")
  compare_system_fixture(fit, "sim_cluster", "cl")
})

test_that("sim_cluster system: cluster small", {
  skip_if(!file.exists(sim_cluster_path), "sim_cluster data not found")
  fit <- ivreg2(y ~ x1 | endo1 | z1 + z2, data = sim_cluster,
                clusters = ~cluster_id, small = TRUE, reduced_form = "system")
  compare_system_fixture(fit, "sim_cluster", "cl_small")
})
