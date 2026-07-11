# ============================================================================
# Tests: Identification Diagnostics (Ticket D2)
# ============================================================================
#
# Anderson LM / Cragg-Donald F (IID) and Kleibergen-Paap rk LM / rk Wald F
# (robust/cluster) underidentification and weak identification tests.
# Verifies against Stata ivreg2 fixtures.

# --- Load datasets ---
sim_multi_path <- fixture_path("sim_multi_endo_data.csv")
if (file.exists(sim_multi_path)) {
  sim_multi <- read.csv(sim_multi_path)
}

data(abdata, package = "ivreg2r")
data(mroz, package = "ivreg2r")

# Shared fits (M-17 hoisting precedent, see test-vcov-cluster.R): the ab_cl
# cluster fit is reused across the four tests below instead of being refit.
fit_ab_cl <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                     clusters = ~id, endog = "w")
fit_ab_cl_small <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                          clusters = ~id, small = TRUE, endog = "w")

# Shared mroz_justid fits (M-17 hoisting precedent): the plain iid and robust
# fits are each reused across several tests below instead of being refit.
fit_mroz_justid        <- ivreg2(mroz_justid_formula, data = mroz)
fit_mroz_justid_robust <- ivreg2(mroz_justid_formula, data = mroz, vcov = "robust")

# overid Anderson LM / CD F / KP rk LM / KP Wald F Stata parity is hf-owned
# (compare_hf asserts idstat/idp/iddf/cdf/widstat for H31/H41 in
# test-helpfile-examples.R; M-13 reuse precedent) -- the justid cells below
# are this family's own anchors.


# ============================================================================
# Anderson LM — IID
# ============================================================================

test_that("Anderson LM matches Stata mroz_justid iid fixture", {
  diag_path <- fixture_path("mroz_justid_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- fit_mroz_justid
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$test_name,
               "Anderson canon. corr. LM statistic")
  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("Anderson LM matches Stata sim_multi_endo iid fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- fixture_path("sim_multi_endo_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})


# ============================================================================
# Cragg-Donald F — IID
# ============================================================================

test_that("Cragg-Donald F matches Stata mroz_justid iid fixture", {
  diag_path <- fixture_path("mroz_justid_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- fit_mroz_justid
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$test_name,
               "Cragg-Donald Wald F statistic")
  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})

test_that("Cragg-Donald F matches Stata sim_multi_endo iid fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- fixture_path("sim_multi_endo_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})


# ============================================================================
# KP rk LM — HC1
# ============================================================================

test_that("KP rk LM matches Stata mroz_justid robust fixture", {
  diag_path <- fixture_path("mroz_justid_diagnostics_robust.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- fit_mroz_justid_robust
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$test_name,
               "Kleibergen-Paap rk LM statistic")
  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("KP rk LM matches Stata sim_multi_endo hc1 fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- fixture_path("sim_multi_endo_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi, vcov = "robust")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})


# ============================================================================
# KP rk Wald F — HC1
# ============================================================================

test_that("KP rk Wald F matches Stata mroz_justid robust fixture", {
  diag_path <- fixture_path("mroz_justid_diagnostics_robust.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- fit_mroz_justid_robust
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id_robust$test_name,
               "Kleibergen-Paap rk Wald F statistic")
  expect_equal(fit$diagnostics$weak_id_robust$stat, fixture$widstat,
               tolerance = stata_tol$stat)
})

test_that("KP rk Wald F matches Stata sim_multi_endo hc1 fixture", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")
  diag_path <- fixture_path("sim_multi_endo_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi, vcov = "robust")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id_robust$stat, fixture$widstat,
               tolerance = stata_tol$stat)
})


# ============================================================================
# KP rk LM — Cluster (family M-13, re-based: abdata H88-minus-gmm2s 2SLS,
# cluster(id))
# ============================================================================

test_that("KP rk LM matches Stata ab_cl cl fixture", {
  diag_path <- fixture_path("ab_cl_diagnostics_cl.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- fit_ab_cl
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$test_name,
               "Kleibergen-Paap rk LM statistic")
  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})


# ============================================================================
# KP rk Wald F — Cluster (family M-13, re-based: abdata ab_cl)
# ============================================================================

test_that("KP rk Wald F matches Stata ab_cl cl fixture", {
  diag_path <- fixture_path("ab_cl_diagnostics_cl.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- fit_ab_cl
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id_robust$stat, fixture$widstat,
               tolerance = stata_tol$stat)
})


# ============================================================================
# CD F always present even with robust/cluster
# ============================================================================

test_that("Cragg-Donald F present alongside KP stats (HC1)", {
  diag_path <- fixture_path("mroz_justid_diagnostics_robust.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- fit_mroz_justid_robust
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$weak_id$test_name,
               "Cragg-Donald Wald F statistic")
})

test_that("Cragg-Donald F present alongside KP stats (cluster)", {
  diag_path <- fixture_path("ab_cl_diagnostics_cl.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- fit_ab_cl
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})


# ============================================================================
# Test name dispatch: IID vs robust
# ============================================================================

test_that("IID reports Anderson LM, robust reports KP rk LM", {
  fit_iid <- fit_mroz_justid
  fit_hc1 <- fit_mroz_justid_robust

  expect_equal(fit_iid$diagnostics$underid$test_name,
               "Anderson canon. corr. LM statistic")
  expect_equal(fit_hc1$diagnostics$underid$test_name,
               "Kleibergen-Paap rk LM statistic")
})


# ============================================================================
# weak_id_robust is NULL for IID
# ============================================================================

test_that("weak_id_robust is NULL for IID models", {
  fit <- fit_mroz_justid
  expect_null(fit$diagnostics$weak_id_robust)
})


# ============================================================================
# small flag invariance: stats identical for small=TRUE vs FALSE
# ============================================================================

test_that("small does not change Anderson LM statistic", {
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz, small = TRUE)
  expect_equal(fit1$diagnostics$underid$stat, fit2$diagnostics$underid$stat)
  expect_equal(fit1$diagnostics$weak_id$stat, fit2$diagnostics$weak_id$stat)
})

test_that("small does not change KP rk LM statistic (HC1)", {
  fit1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz,
                 vcov = "robust", small = FALSE)
  fit2 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz,
                 vcov = "robust", small = TRUE)
  expect_equal(fit1$diagnostics$underid$stat, fit2$diagnostics$underid$stat)
  expect_equal(fit1$diagnostics$weak_id_robust$stat,
               fit2$diagnostics$weak_id_robust$stat)
})

test_that("small does not change KP stats (cluster)", {
  fit1 <- fit_ab_cl
  fit2 <- fit_ab_cl_small
  expect_equal(fit1$diagnostics$underid$stat, fit2$diagnostics$underid$stat)
  expect_equal(fit1$diagnostics$weak_id_robust$stat,
               fit2$diagnostics$weak_id_robust$stat)
})


# ============================================================================
# HC0 = HC1 equivalence: KP stats identical
# ============================================================================

test_that("HC0 and HC1 produce identical KP rk LM and Wald F", {
  fit_hc0 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz,
                     vcov = "robust", small = FALSE)
  fit_hc1 <- ivreg2(lwage ~ exper + expersq | educ | age, data = mroz,
                     vcov = "robust", small = TRUE)

  expect_equal(fit_hc0$diagnostics$underid$stat,
               fit_hc1$diagnostics$underid$stat)
  expect_equal(fit_hc0$diagnostics$weak_id_robust$stat,
               fit_hc1$diagnostics$weak_id_robust$stat)
})


# ============================================================================
# OLS: identification diagnostics are NULL
# ============================================================================

test_that("OLS model has NULL identification diagnostics", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$diagnostics$underid)
  expect_null(fit$diagnostics$weak_id)
  expect_null(fit$diagnostics$weak_id_robust)
})


# ============================================================================
# CD F is VCE-invariant
# ============================================================================

test_that("Cragg-Donald F is identical across VCE types", {
  fit_iid <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                    data = mroz)
  fit_hc1 <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                    data = mroz, vcov = "robust")

  expect_equal(fit_iid$diagnostics$weak_id$stat,
               fit_hc1$diagnostics$weak_id$stat)
})


# ============================================================================
# sim_no_constant: y ~ 0 + x1 | endo1 | z1+z2 (noconstant)
# ============================================================================

sim_noconst_path <- fixture_path("sim_no_constant_data.csv")
if (file.exists(sim_noconst_path)) {
  sim_noconst <- read.csv(sim_noconst_path)
}

test_that("Anderson LM matches Stata sim_no_constant iid fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- fixture_path("sim_no_constant_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("Cragg-Donald F matches Stata sim_no_constant iid fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- fixture_path("sim_no_constant_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst)
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})

test_that("KP rk LM matches Stata sim_no_constant hc1 fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- fixture_path("sim_no_constant_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, vcov = "robust")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$underid$stat, fixture$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$underid$p, fixture$idp,
               tolerance = stata_tol$pval)
  expect_identical(fit$diagnostics$underid$df, as.integer(fixture$iddf))
})

test_that("KP rk Wald F matches Stata sim_no_constant hc1 fixture", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- fixture_path("sim_no_constant_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, vcov = "robust")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id_robust$stat, fixture$widstat,
               tolerance = stata_tol$stat)
})

test_that("small does not change id stats (sim_no_constant)", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")

  fit1 <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, small = FALSE)
  fit2 <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, small = TRUE)
  expect_equal(fit1$diagnostics$underid$stat, fit2$diagnostics$underid$stat)
  expect_equal(fit1$diagnostics$weak_id$stat, fit2$diagnostics$weak_id$stat)
})

test_that("Cragg-Donald F present alongside KP stats sim_no_constant (HC1)", {
  skip_if(!file.exists(sim_noconst_path), "sim_no_constant data not found")
  diag_path <- fixture_path("sim_no_constant_diagnostics_hc1.csv")
  skip_if(!file.exists(diag_path), "diagnostics fixture not found")

  fit <- ivreg2(y ~ 0 + x1 | endo1 | z1 + z2, data = sim_noconst, vcov = "robust")
  fixture <- read_diagnostics(diag_path)

  expect_equal(fit$diagnostics$weak_id$stat, fixture$cdf,
               tolerance = stata_tol$stat)
})


test_that("center does not change the KP identification statistics", {
  # Stata's ranktest never receives centering for the underid or weak-id
  # calls (ivreg2.ado:1639-1650, 1699-1711), so the KP omega behind both
  # statistics is always uncentered -- center = TRUE must be a complete
  # no-op here, exactly as for the redundancy statistic
  # (test-diagnostics-redundancy.R). Tripwire: fails if center forwarding
  # is ever reintroduced at the .compute_id_tests() call site.
  data(card)
  f <- lwage ~ exper + expersq + black + smsa + south | educ | nearc4 + nearc2
  fit_no <- ivreg2(f, data = card, vcov = "robust", center = FALSE)
  fit_yes <- ivreg2(f, data = card, vcov = "robust", center = TRUE)
  expect_equal(fit_yes$diagnostics$underid$stat,
               fit_no$diagnostics$underid$stat, tolerance = 1e-12)
  expect_equal(fit_yes$diagnostics$underid$p,
               fit_no$diagnostics$underid$p, tolerance = 1e-12)
  expect_equal(fit_yes$diagnostics$weak_id_robust$stat,
               fit_no$diagnostics$weak_id_robust$stat, tolerance = 1e-12)
})
