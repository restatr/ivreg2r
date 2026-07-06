# ============================================================================
# Tests: Anderson-Rubin Test (Ticket E3)
# ============================================================================
#
# Weak-instrument-robust test of H0: all endogenous coefficients = 0.
# Verifies against Stata ivreg2 fixtures.

# --- Load datasets ---
sim_multi_path <- fixture_path("sim_multi_endo_data.csv")
if (file.exists(sim_multi_path)) {
  sim_multi <- read.csv(sim_multi_path)
}

data(abdata, package = "ivreg2r")
data(mroz, package = "ivreg2r")


# ============================================================================
# Structural tests (no fixtures)
# ============================================================================

test_that("anderson_rubin is NULL for OLS models", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$diagnostics$anderson_rubin)
})

test_that("anderson_rubin has all expected fields", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz)
  ar <- fit$diagnostics$anderson_rubin
  expect_type(ar, "list")
  expect_named(ar, c("f_stat", "f_p", "f_df1", "f_df2",
                      "chi2_stat", "chi2_p", "chi2_df"),
               ignore.order = TRUE)
})

test_that("small=TRUE and small=FALSE produce identical AR stats (IID)", {
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = TRUE)
  ar1 <- fit1$diagnostics$anderson_rubin
  ar2 <- fit2$diagnostics$anderson_rubin
  expect_equal(ar1$f_stat, ar2$f_stat)
  expect_equal(ar1$chi2_stat, ar2$chi2_stat)
  expect_equal(ar1$f_p, ar2$f_p)
  expect_equal(ar1$chi2_p, ar2$chi2_p)
})

test_that("small=TRUE and small=FALSE produce identical AR stats (HC1)", {
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz,
                 vcov = "robust", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz,
                 vcov = "robust", small = TRUE)
  ar1 <- fit1$diagnostics$anderson_rubin
  ar2 <- fit2$diagnostics$anderson_rubin
  expect_equal(ar1$f_stat, ar2$f_stat)
  expect_equal(ar1$chi2_stat, ar2$chi2_stat)
})

# HC0 and HC1 are numerically identical for AR (the HC1 N/(N-K) correction
# doesn't apply to the AR meat). This test guards against someone accidentally
# introducing a correction inside the AR function.
test_that("HC0 and HC1 produce identical AR stats", {
  fit0 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, vcov = "robust")
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, vcov = "robust")
  ar0 <- fit0$diagnostics$anderson_rubin
  ar1 <- fit1$diagnostics$anderson_rubin
  expect_equal(ar0$f_stat, ar1$f_stat)
  expect_equal(ar0$chi2_stat, ar1$chi2_stat)
  expect_equal(ar0$f_p, ar1$f_p)
  expect_equal(ar0$chi2_p, ar1$chi2_p)
})


# ============================================================================
# Helper: run fixture comparison for one spec/VCE combination
# ============================================================================

check_anderson_rubin <- function(fit, fixture_path, label) {
  fixture <- read_diagnostics(fixture_path)
  ar <- fit$diagnostics$anderson_rubin

  test_that(paste("AR F-stat matches Stata", label), {
    expect_equal(ar$f_stat, fixture$arf,
                 tolerance = stata_tol$stat)
  })

  test_that(paste("AR F p-value matches Stata", label), {
    expect_equal(ar$f_p, fixture$arfp,
                 tolerance = stata_tol$pval)
  })

  test_that(paste("AR chi2 matches Stata", label), {
    expect_equal(ar$chi2_stat, fixture$archi2,
                 tolerance = stata_tol$stat)
  })

  test_that(paste("AR chi2 p-value matches Stata", label), {
    expect_equal(ar$chi2_p, fixture$archi2p,
                 tolerance = stata_tol$pval)
  })

  test_that(paste("AR F df1 matches Stata", label), {
    expect_identical(ar$f_df1, as.integer(fixture$ardf))
  })

  test_that(paste("AR F df2 matches Stata", label), {
    expect_identical(ar$f_df2, as.integer(fixture$ardf_r))
  })

  test_that(paste("AR chi2 df matches Stata", label), {
    expect_identical(ar$chi2_df, as.integer(fixture$ardf))
  })
}


# ============================================================================
# mroz_justid: lwage ~ exper+expersq | educ | age (D5a disclosed
# single-instrument restriction; family M-10 re-base, replacing card_just_id)
# K1=1, L1=1, just-identified
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
    check_anderson_rubin(fit, fixture_file, label)
    check_anderson_rubin(fit_small, fixture_file, paste(label, "small"))

    test_that(paste("Anderson-Rubin direct equality: small vs non-small", label), {
      ar <- fit$diagnostics$anderson_rubin
      ar_small <- fit_small$diagnostics$anderson_rubin
      expect_equal(ar_small$f_stat, ar$f_stat)
      expect_equal(ar_small$f_p, ar$f_p)
      expect_equal(ar_small$chi2_stat, ar$chi2_stat)
      expect_equal(ar_small$chi2_p, ar$chi2_p)
      expect_identical(ar_small$f_df1, ar$f_df1)
      expect_identical(ar_small$f_df2, ar$f_df2)
      expect_identical(ar_small$chi2_df, ar$chi2_df)
    })
  }
}


# ============================================================================
# mroz overid: lwage ~ exper+expersq | educ | age+kidslt6+kidsge6. iid
# anchor = H36 (help.txt:1305, same fit as H31 but with arf/arfp populated,
# unlike H31); robust anchor = H41 (help.txt:1325). Family M-10 re-base,
# replacing card_overid.
# K1=1, L1=3, overidentified
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",    suffix = "H36"),
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
    check_anderson_rubin(fit, fixture_file, label)
    check_anderson_rubin(fit_small, fixture_file, paste(label, "small"))

    test_that(paste("Anderson-Rubin direct equality: small vs non-small", label), {
      ar <- fit$diagnostics$anderson_rubin
      ar_small <- fit_small$diagnostics$anderson_rubin
      expect_equal(ar_small$f_stat, ar$f_stat)
      expect_equal(ar_small$f_p, ar$f_p)
      expect_equal(ar_small$chi2_stat, ar$chi2_stat)
      expect_equal(ar_small$chi2_p, ar$chi2_p)
      expect_identical(ar_small$f_df1, ar$f_df1)
      expect_identical(ar_small$f_df2, ar$f_df2)
      expect_identical(ar_small$chi2_df, ar$chi2_df)
    })
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
    check_anderson_rubin(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_multi_endo: y ~ x1+x2 | endo1+endo2 | z1+z2+z3+z4
# K1=2, L1=4
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
    check_anderson_rubin(fit, fixture_file, label)
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
    check_anderson_rubin(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_no_constant: y ~ 0 + x1 | endo1 | z1+z2 (noconstant)
# ============================================================================

sim_noconst_path <- fixture_path("sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

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
    check_anderson_rubin(fit, fixture_file, label)
  }
}
