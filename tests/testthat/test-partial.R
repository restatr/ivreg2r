# ============================================================================
# Tests: Partial Option (FWL Projection) — Ticket O1 (M-24 fixture re-base)
#
# The M-24 re-base moves the partial()/FWL Stata-parity cells off the Card
# fixtures onto the griliches76 H28-minus-cluster base (help.txt:1253):
# `ivreg2 lw s expr tenure rns smsa _I* (iq=med kww age), partial(_I*)`, with
# cluster(year) dropped. Stata documents partial() with worked examples only
# in the cluster context (H27-H29), so every Stata-parity cell below is a
# D5a option-variation on that base (planning/22-spec-matrix.md).
#
# NOT duplicated here: cluster(year) x partial(_I*) parity is hf-owned
# (hf griliches H28, with the M = 7 clusters vs L = 8 instruments
# rank-deficiency warning as the documented lesson); the H29
# gmm2s+cluster+partial command is a documented help-file bug, also hf-owned
# (Stata rejects it with r(506)). card_partial_all was deleted -- partial(_all)
# Stata parity is owned by M-04's mroz cells, which stay in this file. The
# card cluster(smsa66) cells were deleted outright (the binary M=2 cluster
# anti-pattern on the spec-matrix delete table).
#
# Invariance retirement (byte-identity verified in the retired card fixtures,
# 2026-07-06): nopartialsmall under plain IID leaves e(b) and e(V)
# byte-identical to the basic cell, so the nosmall cell exports diagnostics
# only and the coef/vcov identity is pinned fixture-free. The nosmall x small
# cell is retired compositionally: under IID small the basic and nosmall fits
# share coefficients, rss, and (X'X)^-1, so their VCVs differ by exactly the
# ratio df_r(nosmall)/df_r(basic); pinned as a fixture-free
# correction-factor identity (M-14 precedent).
# ============================================================================

data(griliches, package = "ivreg2r")
data(card, package = "ivreg2r")
data(mroz, package = "ivreg2r")

# Base model for the griliches partial() arc: help-file H28 (help.txt:1253)
# minus cluster(year). The shared gril_formula in helper-fixtures.R is the
# H06 model WITH mrt; this family's base model drops mrt, hence the local
# formula object here.
gril_partial_formula <- lw ~ s + expr + tenure + rns + smsa + factor(year) |
  iq | med + kww + age

# One-part (OLS) variant of the same base model, used by the OLS x partial
# cell (no endogenous regressor, no excluded instruments).
gril_partial_ols_formula <- lw ~ s + expr + tenure + rns + smsa + factor(year)

# Compare partial bookkeeping (sdofminus, partial_ct, partialcons, K, and,
# where posted, the model F degrees of freedom) against a diagnostics
# fixture -- the M-24 addition to the shared diagnostics schema.
expect_partial_bookkeeping <- function(fit, diag_file) {
  dx <- read_diagnostics(fixture_path(diag_file))
  expect_equal(fit$sdofminus, as.integer(dx$sdofminus), info = "sdofminus")
  expect_equal(fit$partial_ct, as.integer(dx$partial_ct), info = "partial_ct")
  expect_equal(as.integer(fit$partialcons), as.integer(dx$partialcons),
               info = "partialcons")
  expect_equal(length(coef(fit)), as.integer(dx$K), info = "K")
  # Stata posts e(df_r) only under small, so the fixture column is populated
  # only in the _small cells; R stores df.residual for every fit.
  if (!is.na(dx$df_r)) {
    expect_identical(fit$df.residual, as.integer(dx$df_r))
  }
  if (!is.na(dx$model_f_df1)) {
    expect_equal(fit$model_f_df1, as.integer(dx$model_f_df1),
                 info = "model F df1")
    expect_equal(fit$model_f_df2, as.integer(dx$model_f_df2),
                 info = "model F df2")
  }
}

# Shared fits (M-17 hoisting precedent): the plain iid fit, the
# partial(_cons) fit, and the unpartialled full-model fit are each reused by
# several sections below (Stata-parity, invariance, metadata, FWL invariance).
fit_basic       <- ivreg2(gril_partial_formula, data = griliches,
                           partial = "factor(year)")
fit_basic_small <- ivreg2(gril_partial_formula, data = griliches,
                           partial = "factor(year)", small = TRUE)
fit_cons        <- ivreg2(gril_partial_formula, data = griliches,
                           partial = "_cons")
fit_full        <- ivreg2(gril_partial_formula, data = griliches)


# ============================================================================
# Stata-parity cells (9 D5a option-variations on the H28-minus-cluster base)
# ============================================================================

test_that("partial(factor(year)) matches Stata -- griliches H28 minus cluster, IID", {
  expect_coef_fixture(fit_basic, "gril_partial_coef_iid.csv")
  expect_vcov_fixture(fit_basic, "gril_partial_vcov_iid.csv")
  expect_diagnostics_fixture(fit_basic, "gril_partial_diagnostics_iid.csv")
  expect_partial_bookkeeping(fit_basic, "gril_partial_diagnostics_iid.csv")

  expect_equal(fit_basic$partial_names, "factor(year)")
  expect_true(fit_basic$partialcons)
  expect_false("(Intercept)" %in% names(coef(fit_basic)))
  expect_setequal(names(coef(fit_basic)),
                  c("iq", "s", "expr", "tenure", "rns", "smsa"))
})

test_that("partial(factor(year)) matches Stata -- griliches H28 minus cluster, IID small", {
  expect_coef_fixture(fit_basic_small, "gril_partial_coef_iid_small.csv")
  expect_vcov_fixture(fit_basic_small, "gril_partial_vcov_iid_small.csv")
  expect_diagnostics_fixture(fit_basic_small,
                              "gril_partial_diagnostics_iid_small.csv")
  expect_partial_bookkeeping(fit_basic_small,
                              "gril_partial_diagnostics_iid_small.csv")
})

test_that("partial(factor(year)) matches Stata -- griliches H28 minus cluster, robust", {
  fit <- ivreg2(gril_partial_formula, data = griliches,
                partial = "factor(year)", vcov = "robust")
  expect_coef_fixture(fit, "gril_partial_coef_robust.csv")
  expect_vcov_fixture(fit, "gril_partial_vcov_robust.csv")
  expect_diagnostics_fixture(fit, "gril_partial_diagnostics_robust.csv")
  expect_partial_bookkeeping(fit, "gril_partial_diagnostics_robust.csv")
})

test_that("partial(factor(year)) matches Stata -- griliches H28 minus cluster, robust small", {
  fit <- ivreg2(gril_partial_formula, data = griliches,
                partial = "factor(year)", vcov = "robust", small = TRUE)
  expect_coef_fixture(fit, "gril_partial_coef_robust_small.csv")
  expect_vcov_fixture(fit, "gril_partial_vcov_robust_small.csv")
  expect_diagnostics_fixture(fit, "gril_partial_diagnostics_robust_small.csv")
  expect_partial_bookkeeping(fit, "gril_partial_diagnostics_robust_small.csv")
})

test_that("partial(factor(year)) matches Stata -- griliches H28 minus cluster, aweight", {
  fit <- ivreg2(gril_partial_formula, data = griliches_awt, weights = awt,
                partial = "factor(year)")
  expect_coef_fixture(fit, "gril_partial_aw_coef_iid.csv")
  expect_vcov_fixture(fit, "gril_partial_aw_vcov_iid.csv")
  expect_diagnostics_fixture(fit, "gril_partial_aw_diagnostics_iid.csv")
  expect_partial_bookkeeping(fit, "gril_partial_aw_diagnostics_iid.csv")
})

test_that("partial(factor(year)) matches Stata -- griliches H28 minus cluster, LIML", {
  fit <- ivreg2(gril_partial_formula, data = griliches,
                partial = "factor(year)", method = "liml")
  expect_coef_fixture(fit, "gril_partial_liml_coef_iid.csv")
  expect_vcov_fixture(fit, "gril_partial_liml_vcov_iid.csv")
  expect_diagnostics_fixture(fit, "gril_partial_liml_diagnostics_iid.csv")
  expect_partial_bookkeeping(fit, "gril_partial_liml_diagnostics_iid.csv")
})

test_that("partial(factor(year)) matches Stata -- griliches H28 minus cluster, GMM2S robust", {
  # Feasible with robust (contrast H29's cluster infeasibility, r(506),
  # owned by the hf suite: #clusters < #moments there, not here).
  fit <- ivreg2(gril_partial_formula, data = griliches,
                partial = "factor(year)", method = "gmm2s", vcov = "robust")
  expect_coef_fixture(fit, "gril_partial_gmm2s_coef_robust.csv")
  expect_vcov_fixture(fit, "gril_partial_gmm2s_vcov_robust.csv")
  expect_diagnostics_fixture(fit, "gril_partial_gmm2s_diagnostics_robust.csv")
  expect_partial_bookkeeping(fit, "gril_partial_gmm2s_diagnostics_robust.csv")
})

test_that("partial(_cons) matches Stata -- griliches H28 minus cluster, IID", {
  expect_coef_fixture(fit_cons, "gril_partial_cons_coef_iid.csv")
  expect_vcov_fixture(fit_cons, "gril_partial_cons_vcov_iid.csv")
  expect_diagnostics_fixture(fit_cons, "gril_partial_cons_diagnostics_iid.csv")
  expect_partial_bookkeeping(fit_cons, "gril_partial_cons_diagnostics_iid.csv")

  # partial(_cons) only demeans -- the year dummies stay as regressors.
  expect_equal(fit_cons$partial_ct, 1L)
  expect_true(fit_cons$partialcons)
  expect_equal(fit_cons$partial_names, character(0L))
  expect_equal(length(coef(fit_cons)), 12L)
  expect_false("(Intercept)" %in% names(coef(fit_cons)))
})

test_that("OLS + partial(factor(year)) matches Stata -- griliches H28 minus cluster", {
  fit <- ivreg2(gril_partial_ols_formula, data = griliches,
                partial = "factor(year)")
  expect_coef_fixture(fit, "gril_partial_ols_coef_iid.csv")
  expect_vcov_fixture(fit, "gril_partial_ols_vcov_iid.csv")

  # No endogenous part, so no id tests are posted and the fixture's id-test
  # fields are empty; the fixture's overid fields carry Stata's degenerate
  # 0/0 sargan sentinel, which the shared helper's df > 0 guard skips. Model
  # F, sigma, rss, r2, and N are real quantities and are asserted.
  expect_diagnostics_fixture(fit, "gril_partial_ols_diagnostics_iid.csv")
  expect_partial_bookkeeping(fit, "gril_partial_ols_diagnostics_iid.csv")
})


# ============================================================================
# Invariance retirement: nopartialsmall (fixture-free identities)
# ============================================================================

fit_nosmall       <- ivreg2(gril_partial_formula, data = griliches,
                             partial = "factor(year)", nopartialsmall = TRUE)
fit_nosmall_small <- ivreg2(gril_partial_formula, data = griliches,
                             partial = "factor(year)", nopartialsmall = TRUE,
                             small = TRUE)

test_that("nopartialsmall leaves coefficients and VCV byte-identical to the basic fit", {
  # Stata evidence: byte-identical e(b)/e(V) in the retired card fixtures
  # (2026-07-06); nopartialsmall only suppresses the sdofminus increment, so
  # it moves df-dependent diagnostics (CD F, model F, R-sq) but not the point
  # estimates or their IID VCV.
  expect_equal(coef(fit_nosmall), coef(fit_basic), tolerance = 0)
  expect_equal(vcov(fit_nosmall), vcov(fit_basic), tolerance = 0)

  expect_diagnostics_fixture(fit_nosmall, "gril_partial_nosmall_diagnostics_iid.csv")
  expect_partial_bookkeeping(fit_nosmall, "gril_partial_nosmall_diagnostics_iid.csv")
})

test_that("nosmall x small compositional identity: VCVs differ by exactly the df ratio", {
  # Under IID small, fit_basic_small and fit_nosmall_small share coefficients,
  # rss, and bread; only df_r differs (sdofminus is zeroed by nopartialsmall),
  # so the VCVs differ by exactly the df ratio (M-14 correction-factor-identity
  # precedent). This is the fixture-free retirement of the old
  # card_partial_nosmall small cell.
  expect_equal(coef(fit_nosmall_small), coef(fit_basic_small), tolerance = 0)

  K <- length(coef(fit_basic_small))
  df_basic   <- nobs(fit_basic_small) - K - fit_basic_small$sdofminus
  df_nosmall <- nobs(fit_nosmall_small) - K
  expect_equal(vcov(fit_nosmall_small) * df_nosmall,
               vcov(fit_basic_small) * df_basic,
               tolerance = 1e-12)
})


# ============================================================================
# M-04 mroz sections: partial(_all) id tests, F5 regression (kept, reworked)
# ============================================================================

test_that("partial(_all) id tests match Stata -- mroz, IID", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz, partial = "_all")
  expect_coef_fixture(fit, "mroz_partial_all_coef_iid.csv")
  expect_diagnostics_fixture(fit, "mroz_partial_all_diagnostics_iid.csv")
  expect_partial_bookkeeping(fit, "mroz_partial_all_diagnostics_iid.csv")
})

test_that("partial(_all) id tests match Stata -- mroz, robust (KP rk)", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz, partial = "_all", vcov = "robust")
  expect_coef_fixture(fit, "mroz_partial_all_coef_robust.csv")
  expect_diagnostics_fixture(fit, "mroz_partial_all_diagnostics_robust.csv")
  expect_partial_bookkeeping(fit, "mroz_partial_all_diagnostics_robust.csv")
})

test_that("noconstant with no exogenous regressors id tests match Stata -- mroz, IID", {
  fit <- ivreg2(lwage ~ 0 | educ | age + kidslt6 + kidsge6, data = mroz)
  expect_coef_fixture(fit, "mroz_nocons_coef_iid.csv")
  expect_diagnostics_fixture(fit, "mroz_nocons_diagnostics_iid.csv")
  expect_partial_bookkeeping(fit, "mroz_nocons_diagnostics_iid.csv")
})

test_that("weighted partial(_all) id tests match Stata -- mroz, IID", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz, weights = hours + 1, partial = "_all")
  expect_coef_fixture(fit, "mroz_partial_all_w_coef_iid.csv")
  expect_diagnostics_fixture(fit, "mroz_partial_all_w_diagnostics_iid.csv")
  expect_partial_bookkeeping(fit, "mroz_partial_all_w_diagnostics_iid.csv")
})

test_that("weighted partial(_all) id tests match Stata -- mroz, robust (KP rk)", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz, weights = hours + 1, partial = "_all",
                vcov = "robust")
  expect_coef_fixture(fit, "mroz_partial_all_w_coef_robust.csv")
  expect_diagnostics_fixture(fit, "mroz_partial_all_w_diagnostics_robust.csv")
  expect_partial_bookkeeping(fit, "mroz_partial_all_w_diagnostics_robust.csv")
})


# ============================================================================
# partial(_all) vs partial(_cons) id-test invariance (fixture-free, griliches)
# ============================================================================

test_that("partial(_all) id tests equal partial(_cons) id tests (FWL invariance)", {
  fit_all <- ivreg2(gril_partial_formula, data = griliches, partial = "_all")

  expect_equal(fit_all$diagnostics$underid$stat, fit_cons$diagnostics$underid$stat)
  expect_equal(fit_all$diagnostics$weak_id$stat, fit_cons$diagnostics$weak_id$stat)
})


# ============================================================================
# '_cons' / '(Intercept)' synonym test (griliches)
# ============================================================================

test_that("partial accepts both '_cons' and '(Intercept)'", {
  fit2 <- ivreg2(gril_partial_formula, data = griliches, partial = "(Intercept)")
  expect_equal(coef(fit_cons), coef(fit2))
  expect_equal(vcov(fit_cons), vcov(fit2))
})


# ============================================================================
# FWL coefficient invariance (pure R tests, griliches)
# ============================================================================

test_that("FWL theorem: partial coefficients match full model coefficients", {
  shared <- intersect(names(coef(fit_full)), names(coef(fit_basic)))
  for (nm in shared) {
    expect_equal(unname(coef(fit_basic)[nm]), unname(coef(fit_full)[nm]),
                 tolerance = 1e-10,
                 info = paste("FWL invariance:", nm))
  }
})

test_that("FWL theorem: partial(_cons) coefficients match full model", {
  shared <- intersect(names(coef(fit_full)), names(coef(fit_cons)))
  for (nm in shared) {
    expect_equal(unname(coef(fit_cons)[nm]), unname(coef(fit_full)[nm]),
                 tolerance = 1e-10,
                 info = paste("FWL cons invariance:", nm))
  }
})


# ============================================================================
# Metadata tests (griliches)
# ============================================================================

test_that("partial_ct and partial_names are correct for a variable subset", {
  # partial(rns smsa) -> partial_ct = 3 (2 vars + cons). The fit_basic and
  # fit_cons bookkeeping is already asserted in their Stata-parity cells.
  fit1 <- ivreg2(gril_partial_formula, data = griliches,
                 partial = c("rns", "smsa"))
  expect_equal(fit1$partial_ct, 3L)
  expect_true(fit1$partialcons)
  expect_equal(sort(fit1$partial_names), sort(c("rns", "smsa")))
})

test_that("no partial -> partial_ct = 0", {
  expect_equal(fit_full$partial_ct, 0L)
  expect_false(fit_full$partialcons)
  expect_equal(fit_full$partial_names, character(0L))
})


# ============================================================================
# Error tests (griliches; iq is the endogenous-variable case)
# ============================================================================

test_that("partial with invalid variable names errors", {
  expect_error(
    ivreg2(gril_partial_formula, data = griliches,
           partial = c("s", "nonexistent")),
    "not in the exogenous regressor list"
  )
})

test_that("partial _cons with noconstant model errors", {
  expect_error(
    ivreg2(lw ~ 0 + s + expr + tenure + rns + smsa + factor(year) |
             iq | med + kww + age,
           data = griliches, partial = "_cons"),
    "without an intercept"
  )
})

test_that("partial must be character or NULL", {
  expect_error(
    ivreg2(gril_partial_formula, data = griliches, partial = 42),
    "character vector or NULL"
  )
})

test_that("nopartialsmall must be logical", {
  expect_error(
    ivreg2(gril_partial_formula, data = griliches, partial = c("s"),
           nopartialsmall = "yes"),
    "TRUE or FALSE"
  )
})

test_that("partial cannot include endogenous variables", {
  expect_error(
    ivreg2(gril_partial_formula, data = griliches, partial = c("iq")),
    "not in the exogenous regressor list"
  )
})


# ============================================================================
# Predict restrictions (griliches)
# ============================================================================

test_that("predict() returns fitted values with message for partial models", {
  expect_message(predict(fit_basic), "non-partialled regressors")
})

test_that("predict(newdata) errors for partial models", {
  expect_error(predict(fit_basic, newdata = griliches[1:5, ]),
               "Cannot predict on new data after partialling")
})


# ============================================================================
# Summary footer (griliches)
# ============================================================================

test_that("summary prints partial footer", {
  out <- capture.output(print(summary(fit_basic)))
  expect_true(any(grepl("Partialled out", out)))
  expect_true(any(grepl("_cons", out)))
  expect_true(any(grepl("partial-model", out)))
})


# ============================================================================
# Broom methods with partial (griliches)
# ============================================================================

test_that("partial_ct is stored on the fitted object", {
  expect_equal(fit_basic$partial_ct, 7L)
})

test_that("tidy works with partial model", {
  td <- tidy(fit_basic)
  expect_equal(nrow(td), 6L)
  expect_true(all(c("iq", "s", "expr") %in% td$term))
})


# ============================================================================
# CUE + partial: warning and non-invariance (kept on card, bundled data)
# ============================================================================

test_that("CUE + partial emits FWL non-invariance warning", {
  expect_warning(
    ivreg2(lwage ~ exper + expersq + black + south + smsa |
             educ | nearc2 + nearc4,
           data = card, method = "cue", partial = "smsa"),
    "FWL invariance does not hold for CUE"
  )
})

test_that("CUE + partial is NOT FWL-invariant with robust VCE", {
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


# ============================================================================
# Interaction with sdofminus parameter (griliches)
# ============================================================================

test_that("partial + sdofminus stack correctly", {
  fit <- ivreg2(gril_partial_formula, data = griliches,
                partial = "factor(year)", sdofminus = 2L)
  # sdofminus = 2 (user) + 7 (partial_ct) = 9
  expect_equal(fit$sdofminus, 9L)
})
