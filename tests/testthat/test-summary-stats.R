# ============================================================================
# Tests: Summary statistics — R², adj R², sigma, RSS (Ticket F2)
# ============================================================================
#
# Fixture-based tests verify r2, r2_a, r2u, r2c, rss, rmse, sigmasq against
# Stata's ivreg2 output. Structural tests verify invariance properties.

# --- Load datasets ---
sim_multi_path <- fixture_path("sim_multi_endo_data.csv")
if (file.exists(sim_multi_path)) {
  sim_multi <- read.csv(sim_multi_path)
}

sim_noconst_path <- fixture_path("sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

data(abdata, package = "ivreg2r")
data(mroz, package = "ivreg2r")


# ============================================================================
# Structural tests (no fixtures)
# ============================================================================

test_that("OLS r.squared matches summary(lm(...))$r.squared", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$r.squared, summary(lm_fit)$r.squared, tolerance = 1e-10)
})

test_that("OLS adj.r.squared matches summary(lm(...))$adj.r.squared", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  lm_fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$adj.r.squared, summary(lm_fit)$adj.r.squared,
               tolerance = 1e-10)
})

test_that("R-squared is invariant to small", {
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = TRUE)
  expect_identical(fit1$r.squared, fit2$r.squared)
  expect_identical(fit1$adj.r.squared, fit2$adj.r.squared)
  expect_identical(fit1$r2u, fit2$r2u)
  expect_identical(fit1$r2c, fit2$r2c)
})

test_that("R-squared is invariant to VCE type", {
  fit_iid <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, vcov = "iid")
  fit_hc  <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz,
                    vcov = "robust", small = TRUE)
  expect_identical(fit_iid$r.squared, fit_hc$r.squared)
  expect_identical(fit_iid$adj.r.squared, fit_hc$adj.r.squared)
  expect_identical(fit_iid$r2u, fit_hc$r2u)
  expect_identical(fit_iid$r2c, fit_hc$r2c)
})

test_that("Sigma changes with small: sigma_small^2 * (N-K) == sigma_large^2 * N", {
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = TRUE)
  expect_equal(fit2$sigma^2 * (fit2$nobs - fit2$rank),
               fit1$sigma^2 * fit1$nobs,
               tolerance = 1e-12)
})

test_that("Negative R-squared for IV is allowed (no error/warning)", {
  # exper/expersq as excluded instruments for educ, controlling only for age,
  # is a deliberately poor specification (weak/irrelevant instruments) chosen
  # to exercise negative-R2 IV behavior; it is not a canonical model.
  expect_no_warning(
    fit <- ivreg2(lwage ~ age | educ | exper + expersq, data = mroz)
  )
  expect_true(fit$r.squared < 0)
})

test_that("mss = tss - rss identity holds", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz)
  tss <- fit$rss / (1 - fit$r.squared)
  expect_equal(fit$mss, tss - fit$rss, tolerance = 1e-12)
})

test_that("With intercept: r2 == r2c; without: r2 == r2u", {
  # With intercept
  fit_int <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_identical(fit_int$r.squared, fit_int$r2c)

  # Without intercept
  fit_noint <- ivreg2(mpg ~ 0 + wt + hp, data = mtcars)
  expect_identical(fit_noint$r.squared, fit_noint$r2u)
})


# ============================================================================
# Helper: run summary stats fixture comparison for one spec/VCE combination
# ============================================================================

check_summary_stats <- function(fit, fixture_path, label) {
  fixture <- read_diagnostics(fixture_path)

  test_that(paste("r.squared matches Stata", label), {
    expect_equal(fit$r.squared, fixture$r2,
                 tolerance = stata_tol$coef)
  })

  test_that(paste("adj.r.squared matches Stata", label), {
    expect_equal(fit$adj.r.squared, fixture$r2_a,
                 tolerance = stata_tol$coef)
  })

  test_that(paste("r2u matches Stata", label), {
    expect_equal(fit$r2u, fixture$r2u,
                 tolerance = stata_tol$coef)
  })

  test_that(paste("r2c matches Stata", label), {
    expect_equal(fit$r2c, fixture$r2c,
                 tolerance = stata_tol$coef)
  })

  test_that(paste("rss matches Stata", label), {
    expect_equal(fit$rss, fixture$rss,
                 tolerance = stata_tol$coef)
  })

  test_that(paste("sigma matches Stata rmse", label), {
    expect_equal(fit$sigma, fixture$rmse,
                 tolerance = stata_tol$coef)
  })

  # mroz fixtures (family M-10) carry no sigmasq column (schema note in
  # planning); skip this assertion there while ab_cl and card_fw (whose
  # fixtures DO carry sigmasq) keep it.
  if (!is.null(fixture$sigmasq)) {
    test_that(paste("sigma^2 matches Stata sigmasq", label), {
      expect_equal(fit$sigma^2, fixture$sigmasq,
                   tolerance = stata_tol$coef)
    })
  }
}


# ============================================================================
# Helper: small-fit summary-stats check (sigma is NOT small-invariant)
# ============================================================================
#
# `small` rescales sigma by the N/(N-K) finite-sample factor, so the small
# fit cannot be checked against the fixture's rmse/sigmasq the way the
# non-small fit is. r2/r2_a/r2u/r2c/rss ARE small-invariant (byte-diff
# evidence, 2026-07-06, that retired the card *_small summary fixtures), so
# those are checked against both the fixture and the non-small fit; sigma is
# instead pinned by the fixture-free N/(N-K) identity.

check_summary_stats_small <- function(fit, fit_small, fixture_path, label) {
  fixture <- read_diagnostics(fixture_path)
  K <- length(coef(fit))

  test_that(paste("small fit r2/r2_a/r2u/r2c/rss match Stata and non-small fit",
                   label), {
    expect_equal(fit_small$r.squared, fixture$r2, tolerance = stata_tol$coef)
    expect_equal(fit_small$adj.r.squared, fixture$r2_a, tolerance = stata_tol$coef)
    expect_equal(fit_small$r2u, fixture$r2u, tolerance = stata_tol$coef)
    expect_equal(fit_small$r2c, fixture$r2c, tolerance = stata_tol$coef)
    expect_equal(fit_small$rss, fixture$rss, tolerance = stata_tol$coef)

    expect_equal(fit_small$r.squared, fit$r.squared)
    expect_equal(fit_small$adj.r.squared, fit$adj.r.squared)
    expect_equal(fit_small$r2u, fit$r2u)
    expect_equal(fit_small$r2c, fit$r2c)
    expect_equal(fit_small$rss, fit$rss)
  })

  test_that(paste("small-fit sigma satisfies the N/(N-K) identity", label), {
    expect_equal(fit_small$sigma^2 * (nobs(fit) - K), fit$sigma^2 * nobs(fit),
                 tolerance = 1e-12)
  })
}


# ============================================================================
# mroz_justid: lwage ~ exper+expersq | educ | age (D5a disclosed
# single-instrument restriction; family M-10 re-base, replacing card_just_id)
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",    suffix = "iid"),
  list(vcov = "robust", suffix = "robust")
)) {
  fixture_file <- fixture_path(
    paste0("mroz_justid_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("mroz_justid", vce_combo$suffix)

  if (file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq | educ | age,
                  data = mroz, vcov = vce_combo$vcov)
    fit_small <- ivreg2(lwage ~ exper + expersq | educ | age,
                         data = mroz, vcov = vce_combo$vcov, small = TRUE)
    check_summary_stats(fit, fixture_file, label)
    check_summary_stats_small(fit, fit_small, fixture_file, paste(label, "small"))
  }
}


# ============================================================================
# mroz overid: lwage ~ exper+expersq | educ | age+kidslt6+kidsge6 (the
# H31/H41 help-file model, help.txt:1274 / help.txt:1325; family M-10
# re-base, replacing card_overid)
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",    suffix = "H31"),
  list(vcov = "robust", suffix = "H41")
)) {
  fixture_file <- fixture_path(
    paste0("hf_mroz_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("mroz_overid", vce_combo$suffix)

  if (file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                  data = mroz, vcov = vce_combo$vcov)
    fit_small <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                         data = mroz, vcov = vce_combo$vcov, small = TRUE)
    check_summary_stats(fit, fixture_file, label)
    check_summary_stats_small(fit, fit_small, fixture_file, paste(label, "small"))
  }
}


# ============================================================================
# sim_multi_endo: y ~ x1+x2 | endo1+endo2 | z1+z2+z3+z4
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "robust",  small = FALSE, suffix = "hc1"),
  list(vcov = "robust",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- fixture_path(
    paste0("sim_multi_endo_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_multi_endo", vce_combo$suffix)

  if (file.exists(sim_multi_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                  data = sim_multi, vcov = vce_combo$vcov,
                  small = vce_combo$small)
    check_summary_stats(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_no_constant: y ~ 0+x1 | endo1 | z1+z2
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "iid",  small = TRUE,  suffix = "iid_small"),
  list(vcov = "robust",  small = FALSE, suffix = "hc1"),
  list(vcov = "robust",  small = TRUE,  suffix = "hc1_small")
)) {
  fixture_file <- fixture_path(
    paste0("sim_no_constant_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_no_constant", vce_combo$suffix)

  if (file.exists(sim_noconst_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2,
                  data = sim_noconst, vcov = vce_combo$vcov,
                  small = vce_combo$small)
    check_summary_stats(fit, fixture_file, label)
  }
}


# ============================================================================
# ab_cl: abdata H88-minus-gmm2s 2SLS, cluster(id) (family M-13, re-based;
# the iid/hc1 rows were retired -- unweighted iid/hc1 parity for these
# statistics is owned by the card blocks above, family M-10)
# ============================================================================

for (vce_combo in ab_cl_vce_cells) {
  fixture_file <- fixture_path(
    paste0("ab_cl_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("ab_cl", vce_combo$suffix)

  if (file.exists(fixture_file)) {
    fit <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                  clusters = ~id, small = vce_combo$small, endog = "w")
    check_summary_stats(fit, fixture_file, label)
  }
}


# ============================================================================
# card_fw: lwage ~ exper+expersq+black+south | educ | nearc4+nearc2,
# fweight fwt = mod(age,5)+1, cluster(region) M=9 (family M-11, re-based;
# card_wt from helper-fixtures.R)
# ============================================================================

for (vce_combo in card_fw_vce_cells) {
  fixture_file <- fixture_path(paste0("card_fw_diagnostics_", vce_combo$suffix, ".csv"))
  label <- paste("card_fw", vce_combo$suffix)

  if (file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                  data = card_wt, weights = fwt, weight_type = "fweight",
                  vcov = vce_combo$vcov, small = vce_combo$small,
                  clusters = vce_combo$clusters)
    check_summary_stats(fit, fixture_file, label)
  }
}
