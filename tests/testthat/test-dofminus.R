# ============================================================================
# Tests: dofminus / sdofminus (Ticket K1, family M-29 re-base)
# ============================================================================
#
# M-29 re-base (2026-07-06): canonical bases per planning/22-spec-matrix.md
# row M-29. All old card_just_id, card_overid, and sim_cluster dofminus
# fixtures (including the M=2 smsa66 cluster cells and their 1e-5 tolerance
# exception) are retired; nothing below references them.
# New base: mroz iid/robust (D5a option-variation on hf H31/H41), ab cluster
# (diagnostics-only, endog(w)), and grunfeld within/fe (small, dofminus=9,
# self-verified against `xtreg, fe` at reldif < 1e-10).
#
# Invariance/identity evidence (byte-diffs and exact ratios measured on the
# retired fixtures, 2026-07-06 -- these justify the fixture-free blocks below):
#   1. Coefficients are invariant to both `small` and `dofminus`/`sdofminus`.
#   2. Under `small`, diagnostics are byte-identical except rmse/sigmasq and
#      F_stat floating-point noise (F ratio exactly 1.0, df identical).
#   3. Exact sigma identity:
#      sigma_small^2 * (N - K - dofminus - sdofminus) == sigma^2 * (N - dofminus)
#   4. iid AND robust small VCV factor:
#      V_small == V * (N - dofminus) / (N - K - dofminus - sdofminus)
#   5. Cluster (non-small) coef and vcov are byte-identical under dofminus
#      (Stata) and `identical()` in R.
#   6. Cluster small VCV factor:
#      V_cl_small == V_cl * (N - 1) / (N - K - sdofminus) * M / (M - 1)
#      -- dofminus does NOT enter the cluster small factor.

data(mroz, package = "ivreg2r")
data(abdata, package = "ivreg2r")
data(grunfeld, package = "ivreg2r")


# ============================================================================
# Block A: mroz iid cell (Stata parity)
# ============================================================================

fit_iid <- ivreg2(mroz_overid_formula, data = mroz,
                   dofminus = 1L, sdofminus = 1L, endog = "educ",
                   vcov = "iid", small = FALSE)

test_that("mroz iid dofminus coef/VCV match Stata (D5a option-variation on H31/H41)", {
  expect_coef_fixture(fit_iid, "mroz_dofminus_coef_iid.csv")
  expect_vcov_fixture(fit_iid, "mroz_dofminus_vcov_iid.csv")
})

test_that("mroz iid dofminus diagnostics match Stata", {
  fx <- read_diagnostics(fixture_path("mroz_dofminus_diagnostics_iid.csv"))
  d <- fit_iid$diagnostics

  expect_equal(fit_iid$sigma, sqrt(fx$sigmasq), tolerance = stata_tol$coef)
  expect_identical(fit_iid$df.residual, as.integer(fx$F_df2))

  expect_equal(fit_iid$r.squared, fx$r2, tolerance = stata_tol$coef)
  expect_equal(fit_iid$adj.r.squared, fx$r2_a, tolerance = stata_tol$coef)

  # iid posts Sargan for overid.
  expect_equal(d$overid$stat, fx$sargan, tolerance = stata_tol$stat)
  expect_equal(d$overid$p, fx$sarganp, tolerance = stata_tol$pval)

  expect_equal(d$underid$stat, fx$idstat, tolerance = stata_tol$stat)
  expect_equal(d$underid$p, fx$idp, tolerance = stata_tol$pval)

  # iid: weak_id is the Cragg-Donald F (no weak_id_robust).
  expect_equal(d$weak_id$stat, fx$cdf, tolerance = stata_tol$stat)

  expect_equal(d$anderson_rubin$f_stat, fx$arf, tolerance = stata_tol$stat)
  expect_equal(d$anderson_rubin$f_p, fx$arfp, tolerance = stata_tol$pval)
  expect_equal(d$anderson_rubin$chi2_stat, fx$archi2, tolerance = stata_tol$stat)
  expect_equal(d$anderson_rubin$chi2_p, fx$archi2p, tolerance = stata_tol$pval)

  expect_equal(d$endogeneity$stat, fx$estat, tolerance = stata_tol$stat)
  expect_equal(d$endogeneity$p, fx$estatp, tolerance = stata_tol$pval)

  expect_equal(fit_iid$model_f, fx$F_stat, tolerance = stata_tol$stat)
  expect_equal(fit_iid$model_f_p, fx$F_p, tolerance = stata_tol$pval)
})


# ============================================================================
# Block B: mroz iid small variant (fixture-free identities)
# ============================================================================

fit_iid_small <- ivreg2(mroz_overid_formula, data = mroz,
                         dofminus = 1L, sdofminus = 1L, endog = "educ",
                         vcov = "iid", small = TRUE)

test_that("mroz iid small: coefficients and VCV invariant identities (facts 1, 3, 4)", {
  N <- nobs(fit_iid)
  K <- length(coef(fit_iid))
  dm <- 1L
  sdm <- 1L

  expect_identical(coef(fit_iid_small), coef(fit_iid))

  expect_equal(fit_iid_small$vcov,
               fit_iid$vcov * (N - dm) / (N - K - dm - sdm),
               tolerance = 1e-12)

  expect_equal(fit_iid_small$sigma^2 * (N - K - dm - sdm),
               fit_iid$sigma^2 * (N - dm),
               tolerance = 1e-12)
})

test_that("mroz iid small: diagnostics invariant to small (fact 2)", {
  d <- fit_iid$diagnostics
  ds <- fit_iid_small$diagnostics

  expect_equal(ds$overid$stat, d$overid$stat)
  expect_equal(ds$overid$p, d$overid$p)
  expect_identical(ds$overid$df, d$overid$df)

  expect_equal(ds$underid$stat, d$underid$stat)
  expect_equal(ds$underid$p, d$underid$p)

  expect_equal(ds$weak_id$stat, d$weak_id$stat)

  expect_equal(ds$anderson_rubin$f_stat, d$anderson_rubin$f_stat)
  expect_equal(ds$anderson_rubin$f_p, d$anderson_rubin$f_p)
  expect_identical(ds$anderson_rubin$f_df1, d$anderson_rubin$f_df1)
  expect_equal(ds$anderson_rubin$chi2_stat, d$anderson_rubin$chi2_stat)
  expect_equal(ds$anderson_rubin$chi2_p, d$anderson_rubin$chi2_p)

  expect_equal(ds$endogeneity$stat, d$endogeneity$stat)
  expect_equal(ds$endogeneity$p, d$endogeneity$p)

  expect_equal(fit_iid_small$model_f, fit_iid$model_f)
})


# ============================================================================
# Block C: mroz robust cell + small variant
# ============================================================================

fit_rob <- ivreg2(mroz_overid_formula, data = mroz,
                   dofminus = 1L, sdofminus = 1L, endog = "educ",
                   vcov = "robust", small = FALSE)
fit_rob_small <- ivreg2(mroz_overid_formula, data = mroz,
                         dofminus = 1L, sdofminus = 1L, endog = "educ",
                         vcov = "robust", small = TRUE)

test_that("mroz robust dofminus coef/VCV match Stata", {
  expect_coef_fixture(fit_rob, "mroz_dofminus_coef_robust.csv")
  expect_vcov_fixture(fit_rob, "mroz_dofminus_vcov_robust.csv")
})

test_that("mroz robust dofminus diagnostics match Stata", {
  fx <- read_diagnostics(fixture_path("mroz_dofminus_diagnostics_robust.csv"))
  d <- fit_rob$diagnostics

  expect_equal(fit_rob$sigma, sqrt(fx$sigmasq), tolerance = stata_tol$coef)
  expect_identical(fit_rob$df.residual, as.integer(fx$F_df2))

  expect_equal(fit_rob$r.squared, fx$r2, tolerance = stata_tol$coef)
  expect_equal(fit_rob$adj.r.squared, fx$r2_a, tolerance = stata_tol$coef)

  # robust posts Hansen J for overid.
  expect_equal(d$overid$stat, fx$j, tolerance = stata_tol$stat)
  expect_equal(d$overid$p, fx$jp, tolerance = stata_tol$pval)

  expect_equal(d$underid$stat, fx$idstat, tolerance = stata_tol$stat)
  expect_equal(d$underid$p, fx$idp, tolerance = stata_tol$pval)

  expect_equal(d$weak_id$stat, fx$cdf, tolerance = stata_tol$stat)
  expect_equal(d$weak_id_robust$stat, fx$widstat, tolerance = stata_tol$stat)

  expect_equal(d$anderson_rubin$f_stat, fx$arf, tolerance = stata_tol$stat)
  expect_equal(d$anderson_rubin$f_p, fx$arfp, tolerance = stata_tol$pval)
  expect_equal(d$anderson_rubin$chi2_stat, fx$archi2, tolerance = stata_tol$stat)
  expect_equal(d$anderson_rubin$chi2_p, fx$archi2p, tolerance = stata_tol$pval)

  expect_equal(d$endogeneity$stat, fx$estat, tolerance = stata_tol$stat)
  expect_equal(d$endogeneity$p, fx$estatp, tolerance = stata_tol$pval)

  expect_equal(fit_rob$model_f, fx$F_stat, tolerance = stata_tol$stat)
  expect_equal(fit_rob$model_f_p, fx$F_p, tolerance = stata_tol$pval)
})

test_that("mroz robust small: coefficients and VCV invariant identities (facts 1, 3, 4)", {
  # Same shared VCV factor as the iid case (fact 4): the factor does not
  # depend on the vcov type, only on N, K, dofminus, sdofminus.
  N <- nobs(fit_rob)
  K <- length(coef(fit_rob))
  dm <- 1L
  sdm <- 1L

  expect_identical(coef(fit_rob_small), coef(fit_rob))

  expect_equal(fit_rob_small$vcov,
               fit_rob$vcov * (N - dm) / (N - K - dm - sdm),
               tolerance = 1e-12)

  expect_equal(fit_rob_small$sigma^2 * (N - K - dm - sdm),
               fit_rob$sigma^2 * (N - dm),
               tolerance = 1e-12)
})

test_that("mroz robust small: diagnostics invariant to small (fact 2)", {
  d <- fit_rob$diagnostics
  ds <- fit_rob_small$diagnostics

  expect_equal(ds$overid$stat, d$overid$stat)
  expect_equal(ds$overid$p, d$overid$p)
  expect_identical(ds$overid$df, d$overid$df)

  expect_equal(ds$underid$stat, d$underid$stat)
  expect_equal(ds$underid$p, d$underid$p)

  expect_equal(ds$weak_id$stat, d$weak_id$stat)
  expect_equal(ds$weak_id_robust$stat, d$weak_id_robust$stat)

  expect_equal(ds$anderson_rubin$f_stat, d$anderson_rubin$f_stat)
  expect_equal(ds$anderson_rubin$f_p, d$anderson_rubin$f_p)
  expect_equal(ds$anderson_rubin$chi2_stat, d$anderson_rubin$chi2_stat)
  expect_equal(ds$anderson_rubin$chi2_p, d$anderson_rubin$chi2_p)

  expect_equal(ds$endogeneity$stat, d$endogeneity$stat)
  expect_equal(ds$endogeneity$p, d$endogeneity$p)

  expect_equal(fit_rob_small$model_f, fit_rob$model_f)
})


# ============================================================================
# Block D: ab cluster cell
# ============================================================================

fit_ab <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                  clusters = ~id, endog = "w",
                  dofminus = 1L, sdofminus = 1L, small = FALSE)
fit_ab_plain <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                        clusters = ~id, endog = "w")

test_that("ab cluster dofminus is invariant to plain fit (fact 5)", {
  # fit_ab_plain (no dofminus) is itself fixture-anchored in
  # test-vcov-cluster.R (ab_cl cells); this identity transitively anchors the
  # estimation surface here, which is why no coef/vcov fixture is shipped for
  # this cell.
  expect_identical(coef(fit_ab), coef(fit_ab_plain))
  expect_identical(fit_ab$vcov, fit_ab_plain$vcov)
})

test_that("ab cluster dofminus diagnostics match Stata", {
  fx <- read_diagnostics(fixture_path("ab_cl_dofminus_diagnostics_cl.csv"))
  d <- fit_ab$diagnostics

  expect_equal(fit_ab$sigma, sqrt(fx$sigmasq), tolerance = stata_tol$coef)
  expect_equal(d$overid$stat, fx$j, tolerance = stata_tol$stat)
  expect_equal(d$overid$p, fx$jp, tolerance = stata_tol$pval)
  expect_equal(d$underid$stat, fx$idstat, tolerance = stata_tol$stat)
  expect_equal(d$underid$p, fx$idp, tolerance = stata_tol$pval)
  expect_equal(d$weak_id$stat, fx$cdf, tolerance = stata_tol$stat)
  expect_equal(d$weak_id_robust$stat, fx$widstat, tolerance = stata_tol$stat)
  expect_equal(d$anderson_rubin$f_stat, fx$arf, tolerance = stata_tol$stat)
  expect_equal(d$anderson_rubin$f_p, fx$arfp, tolerance = stata_tol$pval)
  expect_equal(d$anderson_rubin$chi2_stat, fx$archi2, tolerance = stata_tol$stat)
  expect_equal(d$anderson_rubin$chi2_p, fx$archi2p, tolerance = stata_tol$pval)
  expect_equal(d$endogeneity$stat, fx$estat, tolerance = stata_tol$stat)
  expect_equal(d$endogeneity$p, fx$estatp, tolerance = stata_tol$pval)
  expect_equal(fit_ab$model_f, fx$F_stat, tolerance = stata_tol$stat)
  expect_equal(fit_ab$model_f_p, fx$F_p, tolerance = stata_tol$pval)
  expect_equal(fit_ab$r.squared, fx$r2, tolerance = stata_tol$coef)
  expect_equal(fit_ab$adj.r.squared, fx$r2_a, tolerance = stata_tol$coef)

  # This fixture's purpose: cluster + dofminus DOES move these diagnostics
  # (unlike the byte-identical coef/vcov asserted above).
  expect_identical(fit_ab$n_clusters, as.integer(fx$N_clust))
})

test_that("ab cluster dofminus df.residual follows the cluster rule (M - 1)", {
  expect_identical(fit_ab$df.residual, as.integer(fit_ab$n_clusters - 1L))
})


# ============================================================================
# Block E: ab cluster small variant (fixture-free identities)
# ============================================================================

fit_ab_small <- ivreg2(ab_formula, data = abdata, tvar = "year", ivar = "id",
                        clusters = ~id, endog = "w",
                        dofminus = 1L, sdofminus = 1L, small = TRUE)

test_that("ab cluster small: coefficients and VCV invariant identities (facts 1, 3, 6)", {
  N <- nobs(fit_ab)
  K <- length(coef(fit_ab))
  M <- fit_ab$n_clusters
  sdm <- 1L
  dm <- 1L

  expect_identical(coef(fit_ab_small), coef(fit_ab))

  # dofminus does not enter the cluster small factor (measured on the
  # retired sim and card fixtures).
  expect_equal(fit_ab_small$vcov,
               fit_ab$vcov * (N - 1) / (N - K - sdm) * M / (M - 1),
               tolerance = 1e-12)

  expect_equal(fit_ab_small$sigma^2 * (N - K - dm - sdm),
               fit_ab$sigma^2 * (N - dm),
               tolerance = 1e-12)
})

test_that("ab cluster small: diagnostics invariant to small (fact 2)", {
  d <- fit_ab$diagnostics
  ds <- fit_ab_small$diagnostics

  expect_equal(ds$overid$stat, d$overid$stat)
  expect_equal(ds$overid$p, d$overid$p)
  expect_identical(ds$overid$df, d$overid$df)

  expect_equal(ds$underid$stat, d$underid$stat)
  expect_equal(ds$underid$p, d$underid$p)

  expect_equal(ds$weak_id$stat, d$weak_id$stat)
  expect_equal(ds$weak_id_robust$stat, d$weak_id_robust$stat)

  expect_equal(ds$anderson_rubin$f_stat, d$anderson_rubin$f_stat)
  expect_equal(ds$anderson_rubin$f_p, d$anderson_rubin$f_p)
  expect_equal(ds$anderson_rubin$chi2_stat, d$anderson_rubin$chi2_stat)
  expect_equal(ds$anderson_rubin$chi2_p, d$anderson_rubin$chi2_p)

  expect_equal(ds$endogeneity$stat, d$endogeneity$stat)
  expect_equal(ds$endogeneity$p, d$endogeneity$p)

  expect_equal(fit_ab_small$model_f, fit_ab$model_f)
})


# ============================================================================
# Block F: grunfeld within/fe cell
# ============================================================================

# Mirrors the .do construction exactly: per-company demean with the grand
# mean added back. The Stata side self-verified this cell against
# `xtreg, fe` (slope coefs and SEs at reldif < 1e-10), so the fixture
# anchors the documented xtdata-fe semantics of dofminus, not just recorded
# output.
gr <- grunfeld
for (v in c("invest", "mvalue", "kstock")) {
  gr[[paste0(v, "_w")]] <- gr[[v]] - ave(gr[[v]], gr$company) + mean(gr[[v]])
}
fit_gr <- ivreg2(invest_w ~ mvalue_w + kstock_w, data = gr,
                  small = TRUE, dofminus = 9L)

test_that("grunfeld within/fe dofminus matches Stata (xtreg fe self-verified)", {
  expect_coef_fixture(fit_gr, "grun_fe_dofminus_coef_small.csv")
  expect_vcov_fixture(fit_gr, "grun_fe_dofminus_vcov_small.csv")

  fx <- read_diagnostics(fixture_path("grun_fe_dofminus_diagnostics_small.csv"))
  expect_equal(fit_gr$sigma, sqrt(fx$sigmasq), tolerance = stata_tol$coef)
  expect_equal(fit_gr$r.squared, fx$r2, tolerance = stata_tol$coef)
  expect_equal(fit_gr$adj.r.squared, fx$r2_a, tolerance = stata_tol$coef)
  expect_equal(fit_gr$model_f, fx$F_stat, tolerance = stata_tol$stat)
  expect_equal(fit_gr$model_f_p, fx$F_p, tolerance = stata_tol$pval)

  # N - K - dofminus = 200 - 3 - 9 = 188 = e(df_r) of xtreg fe
  expect_identical(fit_gr$df.residual, 188L)
})


# ============================================================================
# Block G: pure-R behavior tests (ported from the retired card/sim version)
# ============================================================================

test_that("dofminus=0 sdofminus=0 is identical to default", {
  fit_default <- ivreg2(mroz_overid_formula, data = mroz, endog = "educ")
  fit_zero <- ivreg2(mroz_overid_formula, data = mroz, endog = "educ",
                      dofminus = 0L, sdofminus = 0L)

  expect_identical(coef(fit_default), coef(fit_zero))
  expect_identical(fit_default$vcov, fit_zero$vcov)
  expect_identical(fit_default$sigma, fit_zero$sigma)
  expect_identical(fit_default$df.residual, fit_zero$df.residual)
  expect_identical(fit_default$diagnostics$overid$stat,
                    fit_zero$diagnostics$overid$stat)
})

test_that("df.residual = N - K - dofminus - sdofminus (non-cluster)", {
  fit <- ivreg2(mroz_overid_formula, data = mroz, endog = "educ",
                dofminus = 1L, sdofminus = 1L)

  expected_df <- nobs(fit) - length(coef(fit)) - 1L - 1L
  expect_identical(fit$df.residual, expected_df)
})

test_that("dofminus and sdofminus are stored in the object", {
  fit <- ivreg2(mroz_overid_formula, data = mroz, endog = "educ",
                dofminus = 3L, sdofminus = 2L)

  expect_identical(fit$dofminus, 3L)
  expect_identical(fit$sdofminus, 2L)
})

test_that("default dofminus/sdofminus are 0L", {
  fit <- ivreg2(mroz_overid_formula, data = mroz, endog = "educ")

  expect_identical(fit$dofminus, 0L)
  expect_identical(fit$sdofminus, 0L)
})

test_that("dofminus rejects invalid inputs", {
  # Negative values
  expect_error(ivreg2(lwage ~ exper, data = mroz, dofminus = -1L),
               "non-negative")
  expect_error(ivreg2(lwage ~ exper, data = mroz, sdofminus = -1L),
               "non-negative")
  # NA
  expect_error(ivreg2(lwage ~ exper, data = mroz, dofminus = NA_integer_),
               "non-negative")
  # Fractional values (should be rejected, not silently truncated)
  expect_error(ivreg2(lwage ~ exper, data = mroz, dofminus = 1.5),
               "non-negative integer")
  expect_error(ivreg2(lwage ~ exper, data = mroz, sdofminus = 2.9),
               "non-negative integer")
  # Non-numeric
  expect_error(ivreg2(lwage ~ exper, data = mroz, dofminus = "1"),
               "non-negative integer")
  # Inf and overflow doubles
  expect_error(ivreg2(lwage ~ exper, data = mroz, dofminus = Inf),
               "non-negative integer")
  expect_error(ivreg2(lwage ~ exper, data = mroz, sdofminus = Inf),
               "non-negative integer")
  expect_error(ivreg2(lwage ~ exper, data = mroz, dofminus = 1e20),
               "non-negative integer")
})

test_that("dofminus rejects values too large for model dimensions", {
  # mroz has N = 428 complete cases for these formulas.

  # dofminus >= N
  expect_error(ivreg2(lwage ~ exper, data = mroz, dofminus = 5000L),
               "must be less than N")
  # N - K - dofminus - sdofminus <= 0 (OLS: K = 2, intercept + exper)
  expect_error(ivreg2(lwage ~ exper, data = mroz, dofminus = 400L, sdofminus = 30L),
               "too large")
  # IV model: N - K - dofminus - sdofminus <= 0 (K = 4: intercept, exper,
  # expersq, educ)
  expect_error(
    ivreg2(mroz_overid_formula, data = mroz, endog = "educ",
           dofminus = 400L, sdofminus = 30L),
    "too large"
  )
})
