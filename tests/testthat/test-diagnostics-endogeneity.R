# ============================================================================
# Tests: Endogeneity Test / C-Statistic (Ticket E4)
# ============================================================================
#
# Tests H0: specified endogenous regressors are actually exogenous.
# Difference-in-Sargan (IID) or difference-in-J (robust/cluster) statistic.
# Verifies against Stata ivreg2 fixtures (e(estat), e(estatp), e(estatdf)).

# --- Load datasets ---
sim_multi_path <- fixture_path("sim_multi_endo_data.csv")
if (file.exists(sim_multi_path)) {
  sim_multi <- read.csv(sim_multi_path)
}

sim_nc_path <- fixture_path("sim_no_constant_data.csv")
if (file.exists(sim_nc_path)) {
  sim_nc <- read.csv(sim_nc_path)
}

data(abdata, package = "ivreg2r")
data(mroz, package = "ivreg2r")


# ============================================================================
# Structural tests (no fixtures)
# ============================================================================

test_that("endogeneity is NULL for OLS models", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$diagnostics$endogeneity)
})

test_that("endogeneity has all expected fields", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz)
  endog <- fit$diagnostics$endogeneity
  expect_type(endog, "list")
  expect_named(endog, c("stat", "p", "df", "test_name", "tested_vars"),
               ignore.order = TRUE)
  expect_equal(endog$tested_vars, "educ")
  expect_identical(endog$df, 1L)
})

test_that("small=TRUE and small=FALSE produce identical endogeneity stats (IID)", {
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = TRUE)
  expect_equal(fit1$diagnostics$endogeneity$stat,
               fit2$diagnostics$endogeneity$stat)
  expect_equal(fit1$diagnostics$endogeneity$p,
               fit2$diagnostics$endogeneity$p)
})

test_that("small=TRUE and small=FALSE produce identical endogeneity stats (HC1)", {
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz,
                 vcov = "robust", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz,
                 vcov = "robust", small = TRUE)
  expect_equal(fit1$diagnostics$endogeneity$stat,
               fit2$diagnostics$endogeneity$stat)
  expect_equal(fit1$diagnostics$endogeneity$p,
               fit2$diagnostics$endogeneity$p)
})

test_that("HC0 and HC1 produce identical endogeneity stat", {
  fit0 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, vcov = "robust")
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, vcov = "robust")
  expect_equal(fit0$diagnostics$endogeneity$stat,
               fit1$diagnostics$endogeneity$stat)
})

test_that("invalid endog variable name produces error", {
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | age,
           data = mroz, endog = "not_a_var"),
    "not in the endogenous list"
  )
})

test_that("endog = endo_names equals default endog = NULL", {
  fit_default <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz)
  fit_explicit <- ivreg2(lwage ~ exper + expersq | educ | age,
                         data = mroz, endog = "educ")
  expect_equal(fit_default$diagnostics$endogeneity$stat,
               fit_explicit$diagnostics$endogeneity$stat)
  expect_equal(fit_default$diagnostics$endogeneity$p,
               fit_explicit$diagnostics$endogeneity$p)
})

test_that("endog must be character or NULL", {
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | age,
           data = mroz, endog = 1),
    "character vector or NULL"
  )
})

test_that("duplicate endog entries are deduplicated with warning", {
  fit_single <- ivreg2(lwage ~ exper + expersq | educ | age,
                       data = mroz, endog = "educ")
  expect_warning(
    fit_dup <- ivreg2(lwage ~ exper + expersq | educ | age,
                      data = mroz, endog = c("educ", "educ")),
    "duplicate entries"
  )
  expect_equal(fit_dup$diagnostics$endogeneity$stat,
               fit_single$diagnostics$endogeneity$stat)
  expect_equal(fit_dup$diagnostics$endogeneity$df,
               fit_single$diagnostics$endogeneity$df)
})

# ============================================================================
# Estimation invariance: endog() is diagnostics-only (M-10 anchor for
# mroz_endog_diagnostics_robust.csv, whose overid + robust fit is generated
# with `endog(educ)` but whose coef/vcov are identical to the plain H41 fit)
# ============================================================================

test_that("endog = 'educ' does not change coef/vcov of the overid robust fit", {
  fit_plain <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                      data = mroz, vcov = "robust")
  fit_endog <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                      data = mroz, vcov = "robust", endog = "educ")
  expect_identical(coef(fit_plain), coef(fit_endog))
  expect_identical(vcov(fit_plain), vcov(fit_endog))
})


# ============================================================================
# Helper: run fixture comparison for one spec/VCE combination
# ============================================================================

check_endogeneity <- function(fit, fixture_path, label) {
  fixture <- read_diagnostics(fixture_path)
  endog <- fit$diagnostics$endogeneity

  # Skip if Stata didn't compute the stat (e.g., too few clusters)
  if (is.na(fixture$estat)) return(invisible())

  test_that(paste("Endogeneity stat matches Stata", label), {
    expect_equal(endog$stat, fixture$estat,
                 tolerance = stata_tol$stat)
  })

  test_that(paste("Endogeneity p-value matches Stata", label), {
    expect_equal(endog$p, fixture$estatp,
                 tolerance = stata_tol$pval)
  })

  test_that(paste("Endogeneity df matches Stata", label), {
    expect_identical(endog$df, as.integer(fixture$estatdf))
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
                  data = mroz, vcov = vce_combo$vcov, endog = "educ")
    fit_small <- ivreg2(lwage ~ exper + expersq | educ | age,
                         data = mroz, vcov = vce_combo$vcov, small = TRUE,
                         endog = "educ")
    check_endogeneity(fit, fixture_file, label)
    check_endogeneity(fit_small, fixture_file, paste(label, "small"))

    test_that(paste("Endogeneity direct equality: small vs non-small", label), {
      expect_equal(fit_small$diagnostics$endogeneity$stat,
                   fit$diagnostics$endogeneity$stat)
      expect_equal(fit_small$diagnostics$endogeneity$p,
                   fit$diagnostics$endogeneity$p)
      expect_identical(fit_small$diagnostics$endogeneity$df,
                        fit$diagnostics$endogeneity$df)
    })
  }
}


# ============================================================================
# mroz overid: lwage ~ exper+expersq | educ | age+kidslt6+kidsge6. iid
# anchor = H33 (help.txt:1283, H31's fit + endog(educ)); robust anchor =
# mroz_endog_diagnostics_robust.csv (H41's fit + endog(educ)). Family M-10
# re-base, replacing card_overid.
# K1=1, L1=3, overidentified
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",    fixture = "hf_mroz_diagnostics_H33.csv"),
  list(vcov = "robust", fixture = "mroz_endog_diagnostics_robust.csv")
)) {
  fixture_file <- fixture_path(vce_combo$fixture)
  label <- paste("mroz_overid", vce_combo$vcov)

  if (file.exists(fixture_file)) {
    fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                  data = mroz, vcov = vce_combo$vcov, endog = "educ")
    fit_small <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                         data = mroz, vcov = vce_combo$vcov, small = TRUE,
                         endog = "educ")
    check_endogeneity(fit, fixture_file, label)
    check_endogeneity(fit_small, fixture_file, paste(label, "small"))

    test_that(paste("Endogeneity direct equality: small vs non-small", label), {
      expect_equal(fit_small$diagnostics$endogeneity$stat,
                   fit$diagnostics$endogeneity$stat)
      expect_equal(fit_small$diagnostics$endogeneity$p,
                   fit$diagnostics$endogeneity$p)
      expect_identical(fit_small$diagnostics$endogeneity$df,
                        fit$diagnostics$endogeneity$df)
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
    check_endogeneity(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_multi_endo: y ~ x1+x2 | endo1+endo2 | z1+z2+z3+z4
# K1=2, L1=4, tests both endogenous
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
    check_endogeneity(fit, fixture_file, label)
  }
}


# ============================================================================
# sim_no_constant: y ~ 0 + x1 | endo1 | z1+z2 (noconstant)
# ============================================================================

for (vce_combo in list(
  list(vcov = "iid",  small = FALSE, suffix = "iid"),
  list(vcov = "robust",  small = FALSE, suffix = "hc1")
)) {
  fixture_file <- fixture_path(
    paste0("sim_no_constant_diagnostics_", vce_combo$suffix, ".csv")
  )
  label <- paste("sim_no_constant", vce_combo$suffix)

  if (file.exists(sim_nc_path) && file.exists(fixture_file)) {
    fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2,
                  data = sim_nc, vcov = vce_combo$vcov,
                  small = vce_combo$small)
    check_endogeneity(fit, fixture_file, label)
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
    check_endogeneity(fit, fixture_file, label)
  }
}
