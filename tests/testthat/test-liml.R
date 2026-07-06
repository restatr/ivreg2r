# ============================================================================
# Tests: LIML / Fuller / k-class estimation (Ticket H1; M-15 fixture re-base)
#
# The M-15 fixture re-base moved the LIML/Fuller/k-class Stata-parity cells off Card/simulated fixtures and onto klein/abdata fixtures (klein for the native K1=2 multi-endogenous + panel-lag case, abdata for the K1=3 multi-endogenous cluster case). Card is retained only for fixture-free behavioral, error, display, and identity tests.
# ============================================================================

data(card, package = "ivreg2r")
data(klein, package = "ivreg2r")
data(abdata, package = "ivreg2r")

# --- Helper: compare AR LIML overid against an already-read diagnostics data frame. All call sites in this file are overidentified; the fixture-free just-identified test in section 18 owns exact-id (df = 0) behavior. ---
compare_ar_liml_overid <- function(fit, dx,
                                   tol_stat = stata_tol$stat,
                                   tol_pval = stata_tol$pval) {
  aro <- fit$diagnostics$anderson_rubin_overid

  expect_false(is.null(aro), info = "anderson_rubin_overid should not be NULL")

  expect_equal(aro$df, as.integer(dx$arubindf),
               info = "AR LIML overid df mismatch")
  expect_equal(aro$lr_stat, dx$arubin, tolerance = tol_stat,
               info = "AR LR stat mismatch")
  expect_equal(aro$lr_p, dx$arubinp, tolerance = tol_pval,
               info = "AR LR p-value mismatch")
  expect_equal(aro$lin_stat, dx$arubin_lin, tolerance = tol_stat,
               info = "AR lin stat mismatch")
  expect_equal(aro$lin_p, dx$arubin_linp, tolerance = tol_pval,
               info = "AR lin p-value mismatch")
}

# --- Helper: the RSS/sigma/R-squared/model-F core shared by every Stata parity cell in this file, promoted at the M-15 review to replace its triplicated inline block (sections 9, 10, 15-iid). ---
expect_liml_core_diagnostics <- function(fit, dx) {
  expect_equal(fit$rss, dx$rss, tolerance = stata_tol$coef, info = "RSS mismatch")
  expect_equal(fit$sigma, dx$rmse, tolerance = stata_tol$coef, info = "sigma mismatch")
  expect_equal(fit$r.squared, dx$r2, tolerance = stata_tol$coef,
               info = "R-squared mismatch")
  expect_equal(fit$model_f, dx$F_stat, tolerance = stata_tol$stat,
               info = "Model F mismatch")
}


# ============================================================================
# 1. LIML exactly-identified — lambda=1, matches 2SLS (fixture-free)
# ============================================================================

test_that("LIML exactly-identified: lambda=1, matches 2SLS", {
  fit_liml <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, method = "liml")
  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card)

  # lambda must be exactly 1 for just-identified

  expect_equal(fit_liml$lambda, 1.0)
  expect_equal(fit_liml$kclass_value, 1.0)

  # Bit-identical to 2SLS: plain just-id LIML short-circuits in .fit_kclass()
  expect_equal(coef(fit_liml), coef(fit_2sls))
  expect_equal(vcov(fit_liml), vcov(fit_2sls))
  expect_equal(fit_liml$sigma, fit_2sls$sigma)
})


# ============================================================================
# 2. kclass(1) matches 2SLS exactly; kclass(0) matches OLS exactly
# ============================================================================

test_that("kclass(1) matches 2SLS exactly", {
  fit_k1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                   data = card, kclass = 1)
  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card)

  # Near machine-precision agreement (cross-product vs QR paths differ ~1e-11 on the retired card_data.csv export; ~1e-10 on the bundled card dataset now used here, still far tighter than the 1e-6 Stata-parity standard).
  expect_equal(coef(fit_k1), coef(fit_2sls), tolerance = 1e-9)
  expect_equal(vcov(fit_k1), vcov(fit_2sls), tolerance = 1e-9)
  expect_equal(fit_k1$sigma, fit_2sls$sigma, tolerance = 1e-9)
  expect_equal(fit_k1$rss, fit_2sls$rss, tolerance = 1e-9)
})

test_that("kclass(0) matches OLS exactly", {
  fit_k0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                   data = card, kclass = 0)
  fit_ols <- ivreg2(lwage ~ exper + expersq + black + south + educ, data = card)

  # OLS with all regressors (including educ)
  expect_equal(coef(fit_k0)[names(coef(fit_ols))], coef(fit_ols), tolerance = 1e-10)
  expect_equal(fit_k0$sigma, fit_ols$sigma, tolerance = 1e-10)
})


# ============================================================================
# 3. Return object fields (fixture-free)
# ============================================================================

test_that("LIML return object has correct fields", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")

  expect_equal(fit$method, "liml")
  expect_true(is.numeric(fit$lambda))
  expect_true(!is.na(fit$lambda))
  expect_true(fit$lambda > 1)  # overidentified: lambda > 1
  expect_equal(fit$kclass_value, fit$lambda)  # pure LIML: k = lambda
  expect_equal(fit$fuller_parameter, 0)
})

test_that("Fuller return object has correct fields", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 1)

  expect_equal(fit$method, "liml")
  expect_true(!is.na(fit$lambda))
  expect_equal(fit$fuller_parameter, 1)
  # k = lambda - fuller/(N-L)
  expected_k <- fit$lambda - 1 / (nobs(fit) - 7)  # L=7 for this model
  expect_equal(fit$kclass_value, expected_k, tolerance = 1e-12)
})

test_that("kclass return object has correct fields", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, kclass = 0.5)

  expect_equal(fit$method, "kclass")
  expect_true(is.na(fit$lambda))  # no eigenvalue computed
  expect_equal(fit$kclass_value, 0.5)
  expect_equal(fit$fuller_parameter, 0)
})

test_that("2SLS return object has correct method field", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)

  expect_equal(fit$method, "2sls")
  expect_true(is.na(fit$lambda))
  expect_true(is.na(fit$kclass_value))
  expect_equal(fit$fuller_parameter, 0)
})

test_that("OLS return object has correct method field", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)

  expect_equal(fit$method, "ols")
  expect_true(is.na(fit$lambda))
  expect_true(is.na(fit$kclass_value))
  expect_equal(fit$fuller_parameter, 0)
})


# ============================================================================
# 4. Error tests (fixture-free)
# ============================================================================

test_that("coviv must be logical", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "liml", coviv = "yes"),
    "must be TRUE or FALSE"
  )
})

test_that("coviv is silently ignored for 2SLS", {
  expect_no_warning(
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                  data = card, coviv = TRUE)
  )
  expect_false(fit$coviv)
})

test_that("LIML + kclass is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "liml", kclass = 0.5),
    "Cannot specify `kclass` with"
  )
})

test_that("fuller + kclass is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, fuller = 1, kclass = 0.5),
    "Cannot specify both"
  )
})

test_that("fuller < 0 is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, fuller = -1),
    "non-negative"
  )
})

test_that("kclass < 0 is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, kclass = -0.5),
    "non-negative"
  )
})

test_that("LIML without IV is rejected", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, method = "liml"),
    "requires an IV model"
  )
})

test_that("kclass without IV is rejected", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, kclass = 0.5),
    "requires an IV model"
  )
})

test_that("invalid method is rejected", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, method = "gmm"),
    'must be one of'
  )
})

test_that("method='kclass' without kclass value is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "kclass"),
    'requires a numeric'
  )
})

test_that("kclass = Inf is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, kclass = Inf),
    'finite'
  )
})

test_that("fuller = Inf is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, fuller = Inf),
    'finite'
  )
})

test_that("kclass = NaN is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, kclass = NaN),
    'finite'
  )
})

test_that("fuller auto-promotes method to liml", {
  # fuller with default method="2sls" should auto-promote
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 1)
  expect_equal(fit$method, "liml")
})

test_that("kclass auto-promotes method to kclass", {
  # kclass with default method="2sls" should auto-promote
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, kclass = 0.5)
  expect_equal(fit$method, "kclass")
})

test_that("LIML with collinear instruments drops them gracefully", {
  set.seed(42)
  n <- 100
  z1 <- rnorm(n)
  z2 <- z1  # exact collinearity
  x <- z1 + rnorm(n)
  y <- x + rnorm(n)
  d <- data.frame(y = y, x = x, z1 = z1, z2 = z2)
  # Collinear instruments: z2 = z1. Should be dropped with warning,
  # reducing to exact identification (lambda = 1). No cryptic LAPACK error.
  suppressWarnings(
    expect_warning(
      fit <- ivreg2(y ~ 1 | x | z1 + z2, data = d, method = "liml"),
      "collinear"
    )
  )
  expect_equal(fit$lambda, 1.0)
})


# ============================================================================
# 5. Display and broom methods (fixture-free)
# ============================================================================

test_that("print.ivreg2 shows LIML estimation type", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")
  output <- capture.output(print(fit))
  expect_true(any(grepl("LIML Estimation", output)))
})

test_that("print.ivreg2 shows Fuller LIML estimation type", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 1)
  output <- capture.output(print(fit))
  expect_true(any(grepl("Fuller LIML Estimation", output)))
})

test_that("print.ivreg2 shows k-class estimation type", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, kclass = 0.5)
  output <- capture.output(print(fit))
  expect_true(any(grepl("k-class Estimation", output)))
})

test_that("summary.ivreg2 shows lambda and kclass", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")
  output <- capture.output(print(summary(fit)))
  expect_true(any(grepl("lambda:", output)))
  expect_true(any(grepl("kclass:", output)))
})

test_that("glance includes LIML columns", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")
  g <- glance(fit)

  expect_true("method" %in% names(g))
  expect_true("lambda" %in% names(g))
  expect_true("kclass_value" %in% names(g))
  expect_true("fuller_parameter" %in% names(g))
  expect_equal(g$method, "liml")
  expect_false(is.na(g$lambda))
  expect_equal(g$fuller_parameter, 0)
})

test_that("glance includes Fuller columns", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, fuller = 1)
  g <- glance(fit)

  expect_equal(g$method, "liml")
  expect_equal(g$fuller_parameter, 1)
  expect_false(is.na(g$kclass_value))
})


# ============================================================================
# 6. Existing functionality unchanged (fixture: hf_mroz H31, M-10 re-base)
# 2SLS-unchanged guard now anchored on the hf H31 overid fixture
# (help.txt:1274) instead of the retired card_overid cells.
# ============================================================================

test_that("2SLS results unchanged after LIML addition", {
  data(mroz, package = "ivreg2r")
  mroz_overid_formula <- lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6

  fit <- ivreg2(mroz_overid_formula, data = mroz)
  expect_coef_fixture(fit, "hf_mroz_coef_H31.csv")
  expect_vcov_fixture(fit, "hf_mroz_vcov_H31.csv")
})

test_that("OLS results unchanged after LIML addition", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$method, "ols")
  expect_true(is.na(fit$lambda))
})


# ============================================================================
# 7. COVIV behavioral tests (fixture-free)
# ============================================================================

test_that("COVIV changes VCV but not coefficients", {
  fit_no_coviv <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                         data = card, method = "liml")
  fit_coviv <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                      data = card, method = "liml", coviv = TRUE)

  # Coefficients must be identical

  expect_equal(coef(fit_no_coviv), coef(fit_coviv))
  # VCV must differ (LIML overid: k-class bread != 2SLS bread)
  expect_false(isTRUE(all.equal(vcov(fit_no_coviv), vcov(fit_coviv))))
})

test_that("COVIV silently ignored for 2SLS: VCV equals non-COVIV", {
  fit_coviv <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                      data = card, coviv = TRUE)
  fit_plain <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                      data = card)
  expect_false(fit_coviv$coviv)
  expect_equal(vcov(fit_coviv), vcov(fit_plain))
})

test_that("COVIV appears in summary output", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml", coviv = TRUE)
  output <- capture.output(print(summary(fit)))
  expect_true(any(grepl("coviv:", output)))
})

test_that("COVIV appears in glance output", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml", coviv = TRUE)
  g <- glance(fit)
  expect_true("coviv" %in% names(g))
  expect_true(g$coviv)

  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                 data = card, method = "liml")
  g2 <- glance(fit2)
  expect_false(g2$coviv)
})


# ============================================================================
# 8. LIML just-identified — robust/cluster VCV equals 2SLS (identity)
# ============================================================================
# Just-identified LIML has lambda=1 and k=1, so all breads coincide. The robust half stays on Card (fixture-free). The cluster half moved to a just-identified abdata spec at M-15 (K1=3, L1=3), replacing the retired card clusters=~smsa66 cell (M=2 cluster anti-pattern).

test_that("LIML justid robust VCV equals 2SLS robust VCV", {
  # HC0
  fit_liml <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, method = "liml", vcov = "robust")
  fit_2sls <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                     data = card, vcov = "robust")
  expect_equal(vcov(fit_liml), vcov(fit_2sls))

  # Cluster: just-identified abdata spec (K1 = L1 = 3), cluster(id)
  ab_justid <- n ~ 1 | w + k + ys | d(w, 1) + d(k, 1) + d(ys, 1)
  fit_liml_cl <- ivreg2(ab_justid, data = abdata, tvar = "year", ivar = "id",
                        method = "liml", clusters = ~id)
  fit_2sls_cl <- ivreg2(ab_justid, data = abdata, tvar = "year", ivar = "id",
                        clusters = ~id)
  expect_equal(vcov(fit_liml_cl), vcov(fit_2sls_cl))
})


# ============================================================================
# 9. Klein LIML — Stata parity (small=TRUE; IID reuses tsop_klein fixtures)
# ============================================================================

test_that("Klein LIML (iid, small=TRUE): coef, vcov, diagnostics match Stata", {
  skip_if(!file.exists(fixture_path("klein_liml_coef_iid_small.csv")),
          "Klein LIML small fixture not found")

  fit <- ivreg2(klein_formula, data = klein, tvar = "yr", method = "liml",
                small = TRUE)
  expect_coef_fixture(fit, "klein_liml_coef_iid_small.csv")
  expect_vcov_fixture(fit, "klein_liml_vcov_iid_small.csv")

  dx <- read_diagnostics(fixture_path("klein_liml_diagnostics_iid_small.csv"))
  expect_liml_core_diagnostics(fit, dx)
  expect_equal(fit$lambda, dx$lambda, tolerance = stata_tol$coef)
  expect_equal(fit$kclass_value, dx$kclass, tolerance = stata_tol$coef)

  compare_ar_liml_overid(fit, dx)
})


# ============================================================================
# 10. Klein LIML — Stata parity (robust)
# ============================================================================

test_that("Klein LIML (robust): coef, vcov, diagnostics match Stata", {
  skip_if(!file.exists(fixture_path("klein_liml_coef_hc1.csv")),
          "Klein LIML robust fixture not found")

  fit <- ivreg2(klein_formula, data = klein, tvar = "yr", method = "liml",
                vcov = "robust")
  expect_coef_fixture(fit, "klein_liml_coef_hc1.csv")
  expect_vcov_fixture(fit, "klein_liml_vcov_hc1.csv")

  dx <- read_diagnostics(fixture_path("klein_liml_diagnostics_hc1.csv"))
  expect_liml_core_diagnostics(fit, dx)

  # AR overid is IID-only
  expect_null(fit$diagnostics$anderson_rubin_overid)
})


# ============================================================================
# 11. Klein LIML — Stata parity (robust, small=TRUE)
# ============================================================================

test_that("Klein LIML (robust, small=TRUE): coef, vcov match Stata", {
  skip_if(!file.exists(fixture_path("klein_liml_coef_hc1_small.csv")),
          "Klein LIML robust-small fixture not found")

  fit <- ivreg2(klein_formula, data = klein, tvar = "yr", method = "liml",
                vcov = "robust", small = TRUE)
  expect_coef_fixture(fit, "klein_liml_coef_hc1_small.csv")
  expect_vcov_fixture(fit, "klein_liml_vcov_hc1_small.csv")

  # Small-sample model F is the quantity that genuinely varies under small
  dx <- read_diagnostics(fixture_path("klein_liml_diagnostics_hc1_small.csv"))
  expect_equal(fit$model_f, dx$F_stat, tolerance = stata_tol$stat,
               info = "Model F mismatch")
})


# ============================================================================
# 12. Klein LIML COVIV — Stata parity (robust)
# ============================================================================

test_that("Klein LIML COVIV (robust): coef, vcov match Stata; coefs unchanged", {
  skip_if(!file.exists(fixture_path("klein_liml_coviv_coef_hc1.csv")),
          "Klein LIML COVIV fixture not found")

  fit <- ivreg2(klein_formula, data = klein, tvar = "yr", method = "liml",
                vcov = "robust", coviv = TRUE)
  expect_coef_fixture(fit, "klein_liml_coviv_coef_hc1.csv")
  expect_vcov_fixture(fit, "klein_liml_coviv_vcov_hc1.csv")

  # Model F is computed from the COVIV VCV, so it differs from plain robust
  dx <- read_diagnostics(fixture_path("klein_liml_coviv_diagnostics_hc1.csv"))
  expect_equal(fit$model_f, dx$F_stat, tolerance = stata_tol$stat,
               info = "Model F mismatch")

  # coviv changes the VCV bread only, not the point estimates
  fit_no_coviv <- ivreg2(klein_formula, data = klein, tvar = "yr",
                         method = "liml", vcov = "robust")
  expect_equal(coef(fit), coef(fit_no_coviv))
})


# ============================================================================
# 13. Klein Fuller(1) — Stata parity (robust)
# ============================================================================

test_that("Klein Fuller(1) (robust): coef, vcov match Stata; kclass matches", {
  skip_if(!file.exists(fixture_path("klein_fuller1_coef_hc1.csv")),
          "Klein Fuller(1) fixture not found")

  fit <- ivreg2(klein_formula, data = klein, tvar = "yr", method = "liml",
                fuller = 1, vcov = "robust")
  expect_coef_fixture(fit, "klein_fuller1_coef_hc1.csv")
  expect_vcov_fixture(fit, "klein_fuller1_vcov_hc1.csv")

  dx <- read_diagnostics(fixture_path("klein_fuller1_diagnostics_hc1.csv"))
  expect_equal(fit$kclass_value, dx$kclass, tolerance = stata_tol$coef)
})


# ============================================================================
# 14. Klein kclass(1.19) — Stata parity (robust)
# ============================================================================

test_that("Klein kclass(1.19) (robust): coef, vcov match Stata", {
  skip_if(!file.exists(fixture_path("klein_kclass119_coef_hc1.csv")),
          "Klein kclass(1.19) fixture not found")

  fit <- ivreg2(klein_formula, data = klein, tvar = "yr", kclass = 1.19,
                vcov = "robust")
  expect_coef_fixture(fit, "klein_kclass119_coef_hc1.csv")
  expect_vcov_fixture(fit, "klein_kclass119_vcov_hc1.csv")

  dx <- read_diagnostics(fixture_path("klein_kclass119_diagnostics_hc1.csv"))
  expect_equal(fit$kclass_value, dx$kclass, tolerance = stata_tol$coef)
  expect_equal(fit$model_f, dx$F_stat, tolerance = stata_tol$stat,
               info = "Model F mismatch")
})


# ============================================================================
# 15. Klein weighted LIML — Stata parity (analytic weights)
# ============================================================================

test_that("Klein weighted LIML (iid): coef, vcov, diagnostics, AR overid match Stata", {
  skip_if(!file.exists(fixture_path("klein_liml_aw_coef_iid.csv")),
          "Klein weighted LIML fixture not found")

  fit <- ivreg2(klein_formula, data = klein_awt, tvar = "yr", weights = awt,
                method = "liml")
  expect_coef_fixture(fit, "klein_liml_aw_coef_iid.csv")
  expect_vcov_fixture(fit, "klein_liml_aw_vcov_iid.csv")

  dx <- read_diagnostics(fixture_path("klein_liml_aw_diagnostics_iid.csv"))
  expect_liml_core_diagnostics(fit, dx)
  expect_equal(fit$lambda, dx$lambda, tolerance = stata_tol$coef)
  expect_equal(fit$kclass_value, dx$kclass, tolerance = stata_tol$coef)

  compare_ar_liml_overid(fit, dx)
})

test_that("Klein weighted LIML (robust): coef, vcov match Stata", {
  skip_if(!file.exists(fixture_path("klein_liml_aw_coef_hc1.csv")),
          "Klein weighted LIML robust fixture not found")

  fit <- ivreg2(klein_formula, data = klein_awt, tvar = "yr", weights = awt,
                method = "liml", vcov = "robust")
  expect_coef_fixture(fit, "klein_liml_aw_coef_hc1.csv")
  expect_vcov_fixture(fit, "klein_liml_aw_vcov_hc1.csv")

  dx <- read_diagnostics(fixture_path("klein_liml_aw_diagnostics_hc1.csv"))
  expect_equal(fit$model_f, dx$F_stat, tolerance = stata_tol$stat,
               info = "Model F mismatch")
})


# ============================================================================
# 16. abdata LIML — Stata parity (cluster)
# ============================================================================
# This cell carries the LIML x cluster and K1=3 multi-endogenous intersections that the retired sim_multi_endo and card cluster (M=2, clusters=~smsa66) cells used to cover.

for (cell in list(list(suffix = "cl", small = FALSE),
                  list(suffix = "cl_small", small = TRUE))) {
  test_that(paste0("abdata LIML (", cell$suffix,
                   "): coef, vcov, diagnostics match Stata"), {
    coef_file <- paste0("ab_liml_coef_", cell$suffix, ".csv")
    skip_if(!file.exists(fixture_path(coef_file)),
            "abdata LIML cluster fixture not found")

    fit <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                  method = "liml", clusters = ~id, small = cell$small)
    expect_coef_fixture(fit, coef_file)
    expect_vcov_fixture(fit, paste0("ab_liml_vcov_", cell$suffix, ".csv"))

    dx <- read_diagnostics(
      fixture_path(paste0("ab_liml_diagnostics_", cell$suffix, ".csv")))
    expect_equal(nobs(fit), dx$N)
    expect_equal(fit$n_clusters, dx$N_clust)
    expect_equal(fit$lambda, dx$lambda, tolerance = stata_tol$coef)
    expect_equal(fit$diagnostics$overid$stat, dx$j, tolerance = stata_tol$stat)
    expect_equal(fit$diagnostics$underid$stat, dx$idstat,
                 tolerance = stata_tol$stat)
  })
}

test_that("small correction factor is method-uniform across the k-class family", {
  # The small-sample factor is applied in the shared .vcov_from_omega with no method branch, so the retired fuller x small/kclass x small Stata cells are pinned by this uniformity plus the klein_liml small cells; a future method-specific small branch breaks this test.
  ratio_for <- function(method_args) {
    fit <- do.call(ivreg2, c(list(klein_formula, data = klein, tvar = "yr",
                                  vcov = "robust"), method_args))
    fit_small <- do.call(ivreg2, c(list(klein_formula, data = klein,
                                        tvar = "yr", vcov = "robust",
                                        small = TRUE), method_args))
    vcov(fit_small) / vcov(fit)
  }
  r_liml <- ratio_for(list(method = "liml"))
  r_fuller <- ratio_for(list(method = "liml", fuller = 1))
  r_kclass <- ratio_for(list(kclass = 1.19))
  for (r in list(r_liml, r_fuller, r_kclass)) {
    expect_lt(max(r) - min(r), 1e-12)
  }
  expect_equal(mean(r_liml), mean(r_fuller), tolerance = 1e-12)
  expect_equal(mean(r_liml), mean(r_kclass), tolerance = 1e-12)
})


# ============================================================================
# 17. Klein LIML (iid) — extras beyond H72's coverage (reused tsop fixtures)
# ============================================================================
# H72 in test-ts-operators.R asserts coef/vcov/N/lambda/kclass/AR-lr against the same tsop fixtures but fits klein from tsop_klein_data.csv, whereas this test fits the bundled data(klein) (CSV/bundled coherence is pinned by test-bundled-data.R), so it asserts the full AR-overid set on the bundled-data fit plus the Sargan stat.

test_that("Klein LIML (iid): Sargan and AR overid extras match Stata", {
  skip_if(!file.exists(fixture_path("tsop_klein_diagnostics_liml.csv")),
          "tsop_klein LIML fixture not found")

  fit <- ivreg2(klein_formula, data = klein, tvar = "yr", method = "liml")
  dx <- read_diagnostics(fixture_path("tsop_klein_diagnostics_liml.csv"))

  # Sargan, computed from LIML residuals
  expect_equal(fit$diagnostics$overid$stat, dx$sargan, tolerance = stata_tol$stat)

  aro <- fit$diagnostics$anderson_rubin_overid
  expect_equal(aro$lr_stat, dx$arubin, tolerance = stata_tol$stat)
  expect_equal(aro$lr_p, dx$arubinp, tolerance = stata_tol$pval)
  expect_equal(aro$lin_stat, dx$arubin_lin, tolerance = stata_tol$stat)
  expect_equal(aro$lin_p, dx$arubin_linp, tolerance = stata_tol$pval)
})


# ============================================================================
# 18. AR LIML overidentification statistics (Ticket H3)
# ============================================================================

test_that("AR LIML overid: just-identified LIML returns zero-stat placeholder", {
  # Fixture-free: exact identification always yields df = 0, stats = 0, p = NA
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, method = "liml")
  aro <- fit$diagnostics$anderson_rubin_overid
  expect_equal(aro$df, 0L)
  expect_equal(aro$lr_stat, 0)
  expect_equal(aro$lin_stat, 0)
  expect_true(is.na(aro$lr_p))
  expect_true(is.na(aro$lin_p))
})

test_that("AR LIML overid: Klein Fuller(1) overid IID matches Stata", {
  skip_if(!file.exists(fixture_path("tsop_klein_diagnostics_fuller1.csv")),
          "tsop_klein Fuller(1) fixture not found")

  fit <- ivreg2(klein_formula, data = klein, tvar = "yr", method = "liml",
                fuller = 1)
  dx <- read_diagnostics(fixture_path("tsop_klein_diagnostics_fuller1.csv"))
  compare_ar_liml_overid(fit, dx)
})

# Klein is natively multi-endogenous (K1 = 2: profits, wagetot), so the AR overid path for K1 > 1 is already exercised by the Klein LIML/Fuller(1) tests above; the retired Fuller(4) and sim_multi_endo AR tests have no replacement cell.

# "AR LIML overid: NULL for robust LIML" is already covered by section 10's Klein LIML (robust) test on this exact fit; not repeated here.

test_that("AR LIML overid: NULL for cluster LIML", {
  fit <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                method = "liml", clusters = ~id)
  expect_null(fit$diagnostics$anderson_rubin_overid)
})

test_that("AR LIML overid: NULL for kclass", {
  # IID is the meaningful gate here: AR overid must stay NULL because the method is kclass, not because the VCE is robust.
  fit <- ivreg2(klein_formula, data = klein, tvar = "yr", kclass = 1.19)
  expect_null(fit$diagnostics$anderson_rubin_overid)
})

test_that("AR LIML overid: NULL for 2SLS", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card)
  expect_null(fit$diagnostics$anderson_rubin_overid)
})

test_that("AR LIML overid: NULL for OLS", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south, data = card)
  expect_true(is.null(fit$diagnostics) ||
              is.null(fit$diagnostics$anderson_rubin_overid))
})

test_that("AR LIML overid appears in glance()", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, method = "liml")
  gl <- generics::glance(fit)
  expect_true("ar_overid_lr_stat" %in% names(gl))
  expect_true("ar_overid_lr_p" %in% names(gl))
  expect_true("ar_overid_lin_stat" %in% names(gl))
  expect_true("ar_overid_lin_p" %in% names(gl))
  expect_true("ar_overid_df" %in% names(gl))
  expect_false(is.na(gl$ar_overid_lr_stat))
  expect_false(is.na(gl$ar_overid_lin_stat))
})

test_that("AR LIML overid: glance() returns NA for 2SLS", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card)
  gl <- generics::glance(fit)
  expect_true(is.na(gl$ar_overid_lr_stat))
  expect_true(is.na(gl$ar_overid_lin_stat))
})
