# ============================================================================
# Tests: Importance weights (iweight) — Ticket Q2
# ============================================================================

read_coef_fixture <- function(path) {
  d <- read.csv(path)
  nms <- ifelse(d$term == "_cons", "(Intercept)", d$term)
  list(
    estimate  = setNames(d$estimate, nms),
    std_error = setNames(d$std_error, nms)
  )
}

read_diagnostics_fixture <- function(path) {
  read.csv(path)
}

# --- Load Card data ---
card_path <- fixture_path("card_data.csv")
if (file.exists(card_path)) {
  card <- read.csv(card_path)
}


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
  skip_if(!file.exists(card_path))
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc4, data = card,
           weights = wage, weight_type = "iweight", method = "gmm2s"),
    "iweights not allowed with GMM"
  )
})

test_that("iweight blocks GMM via mixed-case method", {
  skip_if(!file.exists(card_path))
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc4, data = card,
           weights = wage, weight_type = "iweight", method = "GMM2S"),
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
  skip_if(!file.exists(card_path))
  S <- diag(5)
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc4, data = card,
           weights = wage, weight_type = "iweight", smatrix = S),
    "iweights not allowed.*smatrix"
  )
})

test_that("iweight blocks wmatrix", {
  skip_if(!file.exists(card_path))
  W <- diag(5)
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc4, data = card,
           weights = wage, weight_type = "iweight", wmatrix = W),
    "iweights not allowed.*wmatrix"
  )
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
  skip_if(!file.exists(card_path))
  coef_path <- fixture_path("card_iweight_ols_coef_iid.csv")
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq, data = card,
                weights = wage, weight_type = "iweight")
  stata_coef <- read_coef_fixture(coef_path)
  stata_diag <- read_diagnostics_fixture(
    fixture_path("card_iweight_ols_diagnostics_iid.csv"))
  stata_vcov <- read_vcov_fixture(
    fixture_path("card_iweight_ols_vcov_iid.csv"))

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
  skip_if(!file.exists(card_path))
  coef_path <- fixture_path("card_iweight_ols_coef_iid_small.csv")
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq, data = card,
                weights = wage, weight_type = "iweight", small = TRUE)
  stata_coef <- read_coef_fixture(coef_path)
  stata_diag <- read_diagnostics_fixture(
    fixture_path("card_iweight_ols_diagnostics_iid_small.csv"))

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
# Section 3: Stata parity — 2SLS just identified with iweight
# ============================================================================

test_that("iweight 2SLS just-id IID matches Stata", {
  skip_if(!file.exists(card_path))
  coef_path <- fixture_path("card_iweight_just_id_coef_iid.csv")
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4, data = card,
                weights = wage, weight_type = "iweight")
  stata_coef <- read_coef_fixture(coef_path)
  stata_diag <- read_diagnostics_fixture(
    fixture_path("card_iweight_just_id_diagnostics_iid.csv"))
  stata_vcov <- read_vcov_fixture(
    fixture_path("card_iweight_just_id_vcov_iid.csv"))

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
})

test_that("iweight 2SLS just-id IID small matches Stata", {
  skip_if(!file.exists(card_path))
  coef_path <- fixture_path("card_iweight_just_id_coef_iid_small.csv")
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4, data = card,
                weights = wage, weight_type = "iweight", small = TRUE)
  stata_coef <- read_coef_fixture(coef_path)
  stata_diag <- read_diagnostics_fixture(
    fixture_path("card_iweight_just_id_diagnostics_iid_small.csv"))

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

  expect_equal(fit$sigma^2, stata_diag$sigmasq, tolerance = stata_tol$coef)
  expect_equal(fit$model_f, stata_diag$F_stat, tolerance = stata_tol$stat)
  expect_equal(fit$model_f_df1, stata_diag$F_df1)
  expect_equal(fit$model_f_df2, stata_diag$F_df2)
})


# ============================================================================
# Section 4: Stata parity — 2SLS overidentified with iweight
# ============================================================================

test_that("iweight 2SLS overid IID matches Stata", {
  skip_if(!file.exists(card_path))
  coef_path <- fixture_path("card_iweight_overid_coef_iid.csv")
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4, data = card,
                weights = wage, weight_type = "iweight")
  stata_coef <- read_coef_fixture(coef_path)
  stata_diag <- read_diagnostics_fixture(
    fixture_path("card_iweight_overid_diagnostics_iid.csv"))
  stata_vcov <- read_vcov_fixture(
    fixture_path("card_iweight_overid_vcov_iid.csv"))

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
  skip_if(!file.exists(card_path))
  coef_path <- fixture_path("card_iweight_overid_coef_iid_small.csv")
  skip_if(!file.exists(coef_path))

  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc2 + nearc4, data = card,
                weights = wage, weight_type = "iweight", small = TRUE)
  stata_coef <- read_coef_fixture(coef_path)
  stata_diag <- read_diagnostics_fixture(
    fixture_path("card_iweight_overid_diagnostics_iid_small.csv"))

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

  expect_equal(fit$sigma^2, stata_diag$sigmasq, tolerance = stata_tol$coef)
  expect_equal(fit$model_f, stata_diag$F_stat, tolerance = stata_tol$stat)
  expect_equal(fit$model_f_df1, stata_diag$F_df1)
  expect_equal(fit$model_f_df2, stata_diag$F_df2)
})
