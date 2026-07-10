# ============================================================================
# Tests: Importance weights (iweight)
# ============================================================================
#
# Importance weights on the mroz base (H31 model, help.txt:1274; Stata
# documents the `iweight` option for `ivreg2` but provides no worked example,
# so cells are option-variations on the H31 base). Weights are synthetic and
# deterministic. The anti-pattern of passing the non-integer `wage` variable
# as an iweight was retired 2026-07-04.

# read_coef_fixture comes from helper-fixtures.R

read_diagnostics_fixture <- function(path) {
  read.csv(path)
}

# --- Load mroz data (H31 canonical base) ---
data(mroz, package = "ivreg2r")
mroz_lw <- mroz[!is.na(mroz$lwage), ]

# Deterministic weights mirroring the Stata generator
# (`gen int iwint = mod(age,5)+1`; `iwfrac` adds 0.7 so
# sum(iwfrac) = 1556.6 is non-integral and floor() binds).
mroz_lw$iwint <- (mroz_lw$age %% 5) + 1
mroz_lw$iwfrac <- (mroz_lw$age %% 5) + 1.7


# ============================================================================
# Section 1: Validation tests
# ============================================================================

test_that("iweight is accepted as a valid weight_type", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  fit <- ivreg2(mpg ~ wt + hp, data = d, weights = w, weight_type = "iweight")
  expect_s3_class(fit, "ivreg2")
  expect_equal(fit$weight_type, "iweight")
})

test_that("iweight blocks robust VCE", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w,
           weight_type = "iweight", vcov = "robust"),
    "iweights not allowed with robust"
  )
})

test_that("iweight blocks cluster VCE", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w,
           weight_type = "iweight", clusters = ~cyl),
    "iweights not allowed with cluster"
  )
})

test_that("iweight blocks kernel-based VCE", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w,
           weight_type = "iweight", kernel = "bartlett", bw = 3,
           tvar = "gear"),
    "iweights not allowed"
  )
})

test_that("iweight blocks GMM estimation", {
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
           data = mroz_lw, weights = iwint, weight_type = "iweight",
           method = "gmm2s"),
    "iweights not allowed with GMM"
  )
})

test_that("iweight blocks GMM via mixed-case method", {
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
           data = mroz_lw, weights = iwint, weight_type = "iweight",
           method = "GMM2S"),
    "iweights not allowed with GMM"
  )
})

test_that("iweight blocks bw without kernel", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w,
           weight_type = "iweight", bw = 3, tvar = "gear"),
    "iweights not allowed with kernel"
  )
})

test_that("iweight blocks smatrix", {
  # H31 model: L = 6 (exog exper, expersq, intercept + excluded age,
  # kidslt6, kidsge6).
  S <- diag(6)
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
           data = mroz_lw, weights = iwint, weight_type = "iweight",
           smatrix = S),
    "iweights not allowed.*smatrix"
  )
})

test_that("iweight blocks wmatrix", {
  W <- diag(6)
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
           data = mroz_lw, weights = iwint, weight_type = "iweight",
           wmatrix = W),
    "iweights not allowed.*wmatrix"
  )
})

test_that("iweight revalidates fuller against weighted N", {
  # With small fractional weights, sum(w) can be much less than nrow.
  # fuller must be validated against the weighted N, not the physical row count.
  d <- data.frame(y = rnorm(50), x = rnorm(50), z = rnorm(50),
                  endo = rnorm(50), w = rep(0.1, 50))
  # sum(w) = 5, L = 3 (intercept + x + z), so N - L = 2. fuller = 4 is invalid.
  # But physical N = 50, so N - L = 47 would pass the pre-weight check.
  expect_error(
    ivreg2(y ~ x | endo | z, data = d, weights = w,
           weight_type = "iweight", method = "liml", fuller = 4),
    "fuller.*must be less than N - L"
  )
})

test_that("iweight fuller uses float N, not floor(N)", {
  # Discriminating test: sum(w) = 5.5 (non-integer), L = 3.
  # float N - L = 2.5; floor(N) - L = 2.
  # fuller = 2 is valid under float N (2 < 2.5) but would be
  # spuriously rejected under floor(N) (2 >= 2).
  set.seed(42)
  d <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100),
                  endo = rnorm(100), w = rep(0.055, 100))
  fit <- ivreg2(y ~ x | endo | z, data = d, weights = w,
                weight_type = "iweight", method = "liml", fuller = 2)
  expect_s3_class(fit, "ivreg2")
})

test_that("iweight does not require integer weights", {
  d <- mtcars
  d$w <- runif(nrow(d), 0.5, 5.5)  # non-integer
  fit <- ivreg2(mpg ~ wt + hp, data = d, weights = w, weight_type = "iweight")
  expect_s3_class(fit, "ivreg2")
})

test_that("iweight N = floor(sum(weights))", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 10)
  fit <- ivreg2(mpg ~ wt + hp, data = d, weights = w, weight_type = "iweight")
  expect_equal(nobs(fit), as.integer(floor(sum(d$w))))
})

test_that("iweight preserves n_physical", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 10)
  fit <- ivreg2(mpg ~ wt + hp, data = d, weights = w, weight_type = "iweight")
  expect_equal(fit$n_physical, nrow(d))
})


# ============================================================================
# Section 2: Stata parity — OLS with iweight
# ============================================================================

# --- OLS, IID, small=FALSE ---
test_that("iweight OLS IID matches Stata", {
  coef_path <- fixture_path("mroz_iweight_ols_coef_iid.csv")
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq, data = mroz_lw,
                weights = iwint, weight_type = "iweight")
  stata_coef <- read_coef_fixture(coef_path)
  stata_diag <- read_diagnostics_fixture(
    fixture_path("mroz_iweight_ols_diagnostics_iid.csv"))
  stata_vcov <- read_vcov_fixture(
    fixture_path("mroz_iweight_ols_vcov_iid.csv"))

  # N
  expect_equal(nobs(fit), stata_diag$N)

  # Coefficients
  for (nm in names(stata_coef$estimate)) {
    expect_equal(coef(fit)[nm], stata_coef$estimate[nm],
                 tolerance = stata_tol$coef, info = paste("coef:", nm))
  }

  # Standard errors
  se_r <- sqrt(diag(vcov(fit)))
  for (nm in names(stata_coef$std_error)) {
    expect_equal(se_r[nm], stata_coef$std_error[nm],
                 tolerance = stata_tol$se, info = paste("SE:", nm))
  }

  # VCV
  expect_vcov_equal(vcov(fit), stata_vcov)

  # Summary stats
  expect_equal(fit$sigma^2, stata_diag$sigmasq, tolerance = stata_tol$coef)
  expect_equal(fit$r.squared, stata_diag$r2, tolerance = stata_tol$coef)
  expect_equal(fit$adj.r.squared, stata_diag$r2_a, tolerance = stata_tol$coef)
  expect_equal(fit$rss, stata_diag$rss, tolerance = stata_tol$coef)
})

# --- OLS, IID, small=TRUE ---
test_that("iweight OLS IID small matches Stata", {
  coef_path <- fixture_path("mroz_iweight_ols_coef_iid_small.csv")
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq, data = mroz_lw,
                weights = iwint, weight_type = "iweight", small = TRUE)
  stata_coef <- read_coef_fixture(coef_path)
  stata_diag <- read_diagnostics_fixture(
    fixture_path("mroz_iweight_ols_diagnostics_iid_small.csv"))
  stata_vcov <- read_vcov_fixture(
    fixture_path("mroz_iweight_ols_vcov_iid_small.csv"))

  # N
  expect_equal(nobs(fit), stata_diag$N)

  # Coefficients and SEs
  for (nm in names(stata_coef$estimate)) {
    expect_equal(coef(fit)[nm], stata_coef$estimate[nm],
                 tolerance = stata_tol$coef, info = paste("coef:", nm))
  }
  se_r <- sqrt(diag(vcov(fit)))
  for (nm in names(stata_coef$std_error)) {
    expect_equal(se_r[nm], stata_coef$std_error[nm],
                 tolerance = stata_tol$se, info = paste("SE:", nm))
  }

  # VCV
  expect_vcov_equal(vcov(fit), stata_vcov)

  # Summary stats
  expect_equal(fit$sigma^2, stata_diag$sigmasq, tolerance = stata_tol$coef)
  expect_equal(fit$r.squared, stata_diag$r2, tolerance = stata_tol$coef)
  expect_equal(fit$adj.r.squared, stata_diag$r2_a, tolerance = stata_tol$coef)

  # F-stat
  expect_equal(fit$model_f, stata_diag$F_stat, tolerance = stata_tol$stat)
  expect_equal(fit$model_f_df1, stata_diag$F_df1)
  expect_equal(fit$model_f_df2, stata_diag$F_df2)
})


# ============================================================================
# Section 3: Stata parity — 2SLS overidentified with iweight
# ============================================================================

test_that("iweight 2SLS overid IID matches Stata", {
  coef_path <- fixture_path("mroz_iweight_overid_coef_iid.csv")
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, weights = iwint, weight_type = "iweight")
  stata_coef <- read_coef_fixture(coef_path)
  stata_diag <- read_diagnostics_fixture(
    fixture_path("mroz_iweight_overid_diagnostics_iid.csv"))
  stata_vcov <- read_vcov_fixture(
    fixture_path("mroz_iweight_overid_vcov_iid.csv"))

  expect_equal(nobs(fit), stata_diag$N)

  for (nm in names(stata_coef$estimate)) {
    expect_equal(coef(fit)[nm], stata_coef$estimate[nm],
                 tolerance = stata_tol$coef, info = paste("coef:", nm))
  }
  se_r <- sqrt(diag(vcov(fit)))
  for (nm in names(stata_coef$std_error)) {
    expect_equal(se_r[nm], stata_coef$std_error[nm],
                 tolerance = stata_tol$se, info = paste("SE:", nm))
  }
  expect_vcov_equal(vcov(fit), stata_vcov)

  expect_equal(fit$sigma^2, stata_diag$sigmasq, tolerance = stata_tol$coef)
  expect_equal(fit$r.squared, stata_diag$r2, tolerance = stata_tol$coef)
  expect_equal(fit$adj.r.squared, stata_diag$r2_a, tolerance = stata_tol$coef)
  expect_equal(fit$rss, stata_diag$rss, tolerance = stata_tol$coef)

  # Diagnostics: Sargan (IID overid)
  if (!is.na(stata_diag$sargan)) {
    expect_equal(fit$diagnostics$overid$stat, stata_diag$sargan,
                 tolerance = stata_tol$stat, info = "Sargan stat")
    expect_equal(fit$diagnostics$overid$p, stata_diag$sarganp,
                 tolerance = stata_tol$pval, info = "Sargan p-value")
  }
})

test_that("iweight 2SLS overid IID small matches Stata", {
  coef_path <- fixture_path("mroz_iweight_overid_coef_iid_small.csv")
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                data = mroz_lw, weights = iwint, weight_type = "iweight",
                small = TRUE)
  stata_coef <- read_coef_fixture(coef_path)
  stata_diag <- read_diagnostics_fixture(
    fixture_path("mroz_iweight_overid_diagnostics_iid_small.csv"))
  stata_vcov <- read_vcov_fixture(
    fixture_path("mroz_iweight_overid_vcov_iid_small.csv"))

  expect_equal(nobs(fit), stata_diag$N)

  for (nm in names(stata_coef$estimate)) {
    expect_equal(coef(fit)[nm], stata_coef$estimate[nm],
                 tolerance = stata_tol$coef, info = paste("coef:", nm))
  }
  se_r <- sqrt(diag(vcov(fit)))
  for (nm in names(stata_coef$std_error)) {
    expect_equal(se_r[nm], stata_coef$std_error[nm],
                 tolerance = stata_tol$se, info = paste("SE:", nm))
  }
  expect_vcov_equal(vcov(fit), stata_vcov)

  expect_equal(fit$sigma^2, stata_diag$sigmasq, tolerance = stata_tol$coef)
  expect_equal(fit$model_f, stata_diag$F_stat, tolerance = stata_tol$stat)
  expect_equal(fit$model_f_df1, stata_diag$F_df1)
  expect_equal(fit$model_f_df2, stata_diag$F_df2)
})


# ============================================================================
# Section 4: Stata parity — fractional-weight N convention
# ============================================================================

test_that("iweight fractional weights pin N = floor(sum(w))", {
  coef_path <- fixture_path("mroz_iweight_frac_coef_iid.csv")
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq, data = mroz_lw,
                weights = iwfrac, weight_type = "iweight")
  stata_coef <- read_coef_fixture(coef_path)
  stata_diag <- read_diagnostics_fixture(
    fixture_path("mroz_iweight_frac_diagnostics_iid.csv"))
  stata_vcov <- read_vcov_fixture(
    fixture_path("mroz_iweight_frac_vcov_iid.csv"))

  # N = floor(sum(w)); previously pinned by the retired wage fixtures.
  # sum(iwfrac) = 1556.6 -> N = 1556.
  expect_equal(nobs(fit), stata_diag$N)
  expect_equal(nobs(fit), as.integer(floor(sum(mroz_lw$iwfrac))))

  for (nm in names(stata_coef$estimate)) {
    expect_equal(coef(fit)[nm], stata_coef$estimate[nm],
                 tolerance = stata_tol$coef, info = paste("coef:", nm))
  }
  se_r <- sqrt(diag(vcov(fit)))
  for (nm in names(stata_coef$std_error)) {
    expect_equal(se_r[nm], stata_coef$std_error[nm],
                 tolerance = stata_tol$se, info = paste("SE:", nm))
  }
  expect_vcov_equal(vcov(fit), stata_vcov)
})


# ============================================================================
# Section 5: Self-verifying identities
# ============================================================================
# The Stata generator asserts the same identities at generation time
# (see generate-iweight-fixtures.do self-check block).

test_that("integral iweight is equivalent to fweight", {
  fit_iw <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                    data = mroz_lw, weights = iwint, weight_type = "iweight")
  fit_fw <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                    data = mroz_lw, weights = iwint, weight_type = "fweight")

  expect_equal(coef(fit_iw), coef(fit_fw), tolerance = 1e-12)
  expect_equal(vcov(fit_iw), vcov(fit_fw), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(nobs(fit_iw), nobs(fit_fw))
  expect_equal(fit_iw$rss, fit_fw$rss, tolerance = 1e-12)
})

test_that("integral iweight is equivalent to physically duplicated rows", {
  dup <- mroz_lw[rep(seq_len(nrow(mroz_lw)), mroz_lw$iwint), ]

  fit_iw  <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                     data = mroz_lw, weights = iwint, weight_type = "iweight")
  fit_dup <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                     data = dup)

  expect_equal(nobs(fit_iw), 1257L)
  expect_equal(nobs(fit_dup), 1257L)

  expect_equal(coef(fit_iw), coef(fit_dup), tolerance = 1e-12)
  expect_equal(vcov(fit_iw), vcov(fit_dup), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(nobs(fit_iw), nobs(fit_dup))
  expect_equal(fit_iw$rss, fit_dup$rss, tolerance = 1e-12)
})
