# ============================================================================
# Tests: Two-step efficient GMM estimation (Ticket N1; M-16 fixture re-base)
#
# The M-16 fixture re-base moved the GMM2S Stata-parity cells off Card/synthetic fixtures onto griliches H06 (robust plus aweight/dofminus/endog option-variations), abdata H88 (cluster plus aweight), and phillips D5a (kernel-without-robust AC) option-variation cells. The canonical gmm2s parity lives elsewhere and is not re-asserted here: gmm2s robust in the hf suite (test-helpfile-examples.R H06), gmm2s + orthog in the hf suite (H25/H26), gmm2s cluster in tsop_ab (test-ts-operators.R H88), and gmm2s + HAC in tsop_phil (H84/H85) plus the M-19 ts_wmat_autobiw cells. Card is retained only for fixture-free identity and behavioral tests.
# ============================================================================

data(card, package = "ivreg2r")
data(griliches, package = "ivreg2r")
data(abdata, package = "ivreg2r")
data(phillips, package = "ivreg2r")

# H06 model formula (help.txt:1154); mirrors the inline formula in
# test-helpfile-examples.R's H06 test -- keep the two in sync.
gril_formula <- lw ~ s + expr + tenure + rns + smsa + factor(year) | iq |
  med + kww + age + mrt

# D5a phillips model, shared with the M-19/tsop_phil GMM2S+HAC cells.
phil_formula <- cinf ~ 1 | unem | l(unem, 1:3)

# Compare the full VCV against a Stata gmm2s vcov fixture. These fixtures
# carry no term column -- just v1..vK columns whose order matches the
# companion coef fixture's term order -- so this reads that coef fixture for
# term order (translating both ts-operator and xi names), reorders vcov(fit)
# to match, and compares the whole matrix in one call.
expect_gmm_vcov_fixture <- function(fit, vcov_file, coef_file, tol = stata_tol$vcov) {
  coef_fx <- read.csv(fixture_path(coef_file))
  term_order <- translate_stata_xi_names(translate_stata_ts_names(coef_fx$term))
  V_r <- vcov(fit)[term_order, term_order]
  V_stata <- unname(as.matrix(read.csv(fixture_path(vcov_file))))
  expect_equal(V_r, V_stata, tolerance = tol, ignore_attr = TRUE)
}


# ============================================================================
# 1. IID GMM2S — algebraic equivalence to 2SLS (fixture-free identity)
# ============================================================================

test_that("GMM2S IID matches 2SLS coefficients (compositional-retirement pin)", {
  # Pins the retired card_overid_gmm2s_coef/vcov_iid(_small) fixtures: the
  # Stata-side coefficient agreement in those CSVs was ~2.7e-11 (2SLS parity
  # itself is owned by M-10), so GMM2S+IID's Stata parity holds transitively
  # through this algebraic identity plus 2SLS's own parity coverage.
  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card)
  fit_gmm <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                    data = card, method = "gmm2s")

  expect_equal(coef(fit_gmm), coef(fit_2sls), tolerance = 1e-10)
  expect_equal(vcov(fit_gmm), vcov(fit_2sls), tolerance = 1e-10)
  expect_equal(fit_gmm$diagnostics$overid$stat, fit_2sls$diagnostics$overid$stat,
               tolerance = 1e-10)
  expect_equal(fit_gmm$diagnostics$overid$test_name, "Sargan")
  expect_equal(fit_2sls$diagnostics$overid$test_name, "Sargan")

  fit_2sls_small <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                           data = card, small = TRUE)
  fit_gmm_small <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                          data = card, method = "gmm2s", small = TRUE)
  expect_equal(coef(fit_gmm_small), coef(fit_2sls_small), tolerance = 1e-10)
  expect_equal(vcov(fit_gmm_small), vcov(fit_2sls_small), tolerance = 1e-10)
})


# ============================================================================
# 2. Robust GMM2S — coefficients differ from 2SLS (behavioral, fixture-free)
# ============================================================================

test_that("GMM2S robust coefficients differ from 2SLS", {
  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "robust")
  fit_gmm <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                    data = card, method = "gmm2s", vcov = "robust")

  # Coefficients should differ for overidentified + heteroskedastic
  expect_false(isTRUE(all.equal(coef(fit_2sls)["educ"],
                                coef(fit_gmm)["educ"],
                                tolerance = 1e-4)))
})


# ============================================================================
# 3. Just-identified GMM2S — equals 2SLS regardless of omega (fixture-free)
# ============================================================================

test_that("Just-identified GMM2S equals 2SLS coefficients", {
  # Pins the retired card_justid_gmm2s_coef/vcov_robust(_small) fixtures: the
  # Stata-side coefficient agreement in those CSVs was ~7e-12.
  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, vcov = "robust")
  fit_gmm <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card, method = "gmm2s", vcov = "robust")

  expect_equal(coef(fit_gmm), coef(fit_2sls), tolerance = 1e-10)
  expect_equal(vcov(fit_gmm), vcov(fit_2sls), tolerance = 1e-10)
})

test_that("Just-identified GMM2S J stat is zero", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, method = "gmm2s", vcov = "robust")

  expect_equal(fit$diagnostics$overid$stat, 0)
  expect_equal(fit$diagnostics$overid$df, 0L)
  expect_true(is.na(fit$diagnostics$overid$p))
})


# ============================================================================
# Retired-fixture disposition (M-16 re-base)
# ============================================================================
# gmm2s+orthog parity lives in the hf suite (test-helpfile-examples.R H25/H26); gmm2s+HAC lives in tsop_phil (test-ts-operators.R H84/H85) plus the M-19 ts_wmat_autobiw cells; gmm2s robust lives in the hf suite H06; gmm2s cluster lives in tsop_ab (test-ts-operators.R H88). No dedicated orthog test previously existed in this file to delete; this comment documents the owner for anyone searching here.


# ============================================================================
# 4. Griliches H06 GMM2S parity (M-16 re-base)
# ============================================================================

test_that("Griliches GMM2S robust small: coef, vcov, diagnostics match Stata", {
  skip_if(!file.exists(fixture_path("gril_gmm2s_coef_robust_small.csv")),
          "GMM2S robust-small fixture not found")
  fit <- ivreg2(gril_formula, data = griliches, method = "gmm2s",
                vcov = "robust", small = TRUE)
  expect_coef_fixture(fit, "gril_gmm2s_coef_robust_small.csv")
  expect_gmm_vcov_fixture(fit, "gril_gmm2s_vcov_robust_small.csv",
                          "gril_gmm2s_coef_robust_small.csv")

  dx <- read_diagnostics(fixture_path("gril_gmm2s_diagnostics_robust_small.csv"))
  expect_equal(fit$sigma, dx$sigma, tolerance = stata_tol$coef, info = "sigma")
  expect_equal(fit$rss, dx$rss, tolerance = stata_tol$coef, info = "rss")
  expect_equal(fit$model_f, dx$model_f, tolerance = stata_tol$stat, info = "model F")
  expect_equal(fit$diagnostics$overid$stat, dx$overid_stat,
               tolerance = stata_tol$stat, info = "Hansen J stat")
  expect_equal(fit$diagnostics$overid$p, dx$overid_p,
               tolerance = stata_tol$pval, info = "Hansen J p")
})

test_that("Griliches GMM2S aweight robust: coef, vcov, diagnostics match Stata; pweight equivalence", {
  skip_if(!file.exists(fixture_path("gril_gmm2s_aw_coef_robust.csv")),
          "GMM2S aweight fixture not found")
  fit_aw <- ivreg2(gril_formula, data = griliches_awt, method = "gmm2s",
                   vcov = "robust", weights = awt)
  expect_coef_fixture(fit_aw, "gril_gmm2s_aw_coef_robust.csv")
  expect_gmm_vcov_fixture(fit_aw, "gril_gmm2s_aw_vcov_robust.csv",
                          "gril_gmm2s_aw_coef_robust.csv")

  dx <- read_diagnostics(fixture_path("gril_gmm2s_aw_diagnostics_robust.csv"))
  expect_equal(fit_aw$sigma, dx$sigma, tolerance = stata_tol$coef, info = "sigma")
  expect_equal(fit_aw$rss, dx$rss, tolerance = stata_tol$coef, info = "rss")
  expect_equal(fit_aw$diagnostics$overid$stat, dx$overid_stat,
               tolerance = stata_tol$stat, info = "overid stat")

  # Stata `[aw] gmm2s robust` == `[pw] gmm2s` was verified byte-identical in
  # the retired card fixtures (M-11 invariance rule), so the pweight variant
  # is pinned by this direct equality plus the aweight fixture above.
  fit_pw <- ivreg2(gril_formula, data = griliches_awt, method = "gmm2s",
                   weights = awt, weight_type = "pweight")
  expect_equal(coef(fit_pw), coef(fit_aw))
  expect_equal(vcov(fit_pw), vcov(fit_aw))
  expect_equal(fit_pw$diagnostics$overid$stat, fit_aw$diagnostics$overid$stat)
})

test_that("Griliches GMM2S dofminus robust: coef, vcov, diagnostics match Stata", {
  skip_if(!file.exists(fixture_path("gril_gmm2s_dof_coef_robust.csv")),
          "GMM2S dofminus fixture not found")
  fit <- ivreg2(gril_formula, data = griliches, method = "gmm2s",
                vcov = "robust", dofminus = 2L)
  expect_coef_fixture(fit, "gril_gmm2s_dof_coef_robust.csv")
  expect_gmm_vcov_fixture(fit, "gril_gmm2s_dof_vcov_robust.csv",
                          "gril_gmm2s_dof_coef_robust.csv")

  dx <- read_diagnostics(fixture_path("gril_gmm2s_dof_diagnostics_robust.csv"))
  expect_equal(fit$sigma, dx$sigma, tolerance = stata_tol$coef,
               info = "sigma with dofminus")
  expect_equal(fit$rss, dx$rss, tolerance = stata_tol$coef,
               info = "rss with dofminus")
  expect_equal(fit$diagnostics$overid$stat, dx$overid_stat,
               tolerance = stata_tol$stat, info = "overid stat with dofminus")
})

test_that("Griliches GMM2S endog(iq) robust: endogeneity test, coef, vcov match Stata", {
  skip_if(!file.exists(fixture_path("gril_gmm2s_endog_coef_robust.csv")),
          "GMM2S endog fixture not found")
  fit <- ivreg2(gril_formula, data = griliches, method = "gmm2s",
                vcov = "robust", endog = "iq")
  expect_coef_fixture(fit, "gril_gmm2s_endog_coef_robust.csv")
  expect_gmm_vcov_fixture(fit, "gril_gmm2s_endog_vcov_robust.csv",
                          "gril_gmm2s_endog_coef_robust.csv")

  dx <- read_diagnostics(fixture_path("gril_gmm2s_endog_diagnostics_robust.csv"))
  expect_equal(fit$diagnostics$endogeneity$stat, dx$endog_stat,
               tolerance = stata_tol$stat, info = "endog stat")
  expect_equal(fit$diagnostics$endogeneity$p, dx$endog_p,
               tolerance = stata_tol$pval, info = "endog p")
  expect_equal(fit$diagnostics$endogeneity$df, dx$endog_df, info = "endog df")
})


# ============================================================================
# 5. abdata H88 GMM2S parity (M-16 re-base)
# ============================================================================

test_that("abdata GMM2S cluster small: coef, vcov, diagnostics match Stata", {
  skip_if(!file.exists(fixture_path("ab_gmm2s_coef_cl_small.csv")),
          "GMM2S cluster-small fixture not found")
  fit <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                method = "gmm2s", clusters = ~id, small = TRUE)
  expect_coef_fixture(fit, "ab_gmm2s_coef_cl_small.csv")
  expect_gmm_vcov_fixture(fit, "ab_gmm2s_vcov_cl_small.csv",
                          "ab_gmm2s_coef_cl_small.csv")

  dx <- read_diagnostics(fixture_path("ab_gmm2s_diagnostics_cl_small.csv"))
  expect_equal(nobs(fit), dx$N, info = "N")
  expect_equal(fit$model_f, dx$model_f, tolerance = stata_tol$stat, info = "model F")
  expect_equal(fit$diagnostics$overid$stat, dx$overid_stat,
               tolerance = stata_tol$stat, info = "overid stat")
  expect_equal(fit$diagnostics$overid$p, dx$overid_p,
               tolerance = stata_tol$pval, info = "overid p")
  expect_equal(fit$diagnostics$underid$stat, dx$underid_stat,
               tolerance = stata_tol$stat, info = "underid stat")
})

test_that("abdata GMM2S aweight cluster: coef, vcov, diagnostics match Stata", {
  skip_if(!file.exists(fixture_path("ab_gmm2s_aw_coef_cl.csv")),
          "GMM2S aweight-cluster fixture not found")
  fit <- ivreg2(ab_formula, data = abdata_awt, tvar = "year", ivar = "id",
                method = "gmm2s", clusters = ~id, weights = abwt)
  expect_coef_fixture(fit, "ab_gmm2s_aw_coef_cl.csv")
  expect_gmm_vcov_fixture(fit, "ab_gmm2s_aw_vcov_cl.csv",
                          "ab_gmm2s_aw_coef_cl.csv")

  dx <- read_diagnostics(fixture_path("ab_gmm2s_aw_diagnostics_cl.csv"))
  expect_equal(nobs(fit), dx$N, info = "N")
  expect_equal(fit$diagnostics$overid$stat, dx$overid_stat,
               tolerance = stata_tol$stat, info = "overid stat")
  expect_equal(fit$sigma, dx$sigma, tolerance = stata_tol$coef, info = "sigma")
})


# ============================================================================
# 6. Phillips GMM2S-AC parity (M-16 re-base)
# ============================================================================
# Kernel without robust must dispatch the AC (homoskedastic Kronecker) omega.
# This cell is the regression guard for the GMM2S-AC weighting-matrix bug
# (commit 0a405e8), which once produced ~15% coefficient errors.

test_that("Phillips GMM2S AC bw=3: coef, vcov, diagnostics match Stata", {
  skip_if(!file.exists(fixture_path("phil_gmm2s_coef_ac_bw3.csv")),
          "GMM2S AC fixture not found")
  fit <- ivreg2(phil_formula, data = phillips, tvar = "year",
                method = "gmm2s", kernel = "bartlett", bw = 3)
  expect_coef_fixture(fit, "phil_gmm2s_coef_ac_bw3.csv")
  expect_gmm_vcov_fixture(fit, "phil_gmm2s_vcov_ac_bw3.csv",
                          "phil_gmm2s_coef_ac_bw3.csv")

  dx <- read_diagnostics(fixture_path("phil_gmm2s_diagnostics_ac_bw3.csv"))
  expect_equal(fit$diagnostics$overid$stat, dx$overid_stat,
               tolerance = stata_tol$stat, info = "overid stat")
  expect_equal(fit$sigma, dx$sigma, tolerance = stata_tol$coef, info = "sigma")
})


# ============================================================================
# 7. Input validation (fixture-free)
# ============================================================================

test_that("GMM2S requires IV model", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, method = "gmm2s"),
    "requires an IV model"
  )
})

test_that("GMM2S incompatible with fuller", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "gmm2s", fuller = 1),
    "Cannot specify.*fuller.*gmm2s"
  )
})

test_that("GMM2S incompatible with kclass", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "gmm2s", kclass = 0.5),
    "Cannot specify.*kclass.*gmm2s"
  )
})

test_that("GMM2S method stored in return object", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "robust")
  expect_equal(fit$method, "gmm2s")
})


# ============================================================================
# 8. Display (fixture-free)
# ============================================================================

test_that("GMM2S print shows correct estimation label", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s", vcov = "robust")
  output <- capture.output(print(summary(fit)))
  expect_true(any(grepl("2-Step GMM Estimation", output)))
  expect_true(any(grepl("heteroskedasticity", output)))
})


# ============================================================================
# 9. Stock-Yogo tables (should use IV tables for gmm2s; fixture-free)
# ============================================================================

test_that("GMM2S gets IV Stock-Yogo tables", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "gmm2s")
  sy <- fit$diagnostics$weak_id_sy
  expect_false(is.null(sy))
  expect_true("IV size" %in% sy$type)
})
