# ============================================================================
# Tests: Empty-endogenous / HOLS syntax — Ticket F3
# ============================================================================
# Stata's `(=z1 z2)` form: no endogenous regressors, excluded instruments as
# surplus moment conditions. R formula form: y ~ exog | 0 | excluded_iv.
# Estimation is exact OLS (Stata posts e(model)="ols", ivreg2.ado:2101-2106);
# with method="gmm2s" it is Cragg's (1983) HOLS. The surplus conditions drive
# the Sargan/Hansen J, orthog() C-tests, and the J-equivalence identities.

data(mroz, package = "ivreg2r")
mroz_lw <- mroz[!is.na(mroz$lwage), ]   # the lwage estimation sample (N = 428)

hols_fml <- lwage ~ exper + expersq + educ | 0 | age + kidslt6 + kidsge6
ols_fml  <- lwage ~ exper + expersq + educ


# ============================================================================
# Estimation: K1 = 0 is exact OLS
# ============================================================================

test_that("empty-endogenous fit is exact OLS (iid)", {
  f_hols <- ivreg2(hols_fml, data = mroz)
  f_ols  <- ivreg2(ols_fml, data = mroz)

  # Same .fit_ols code path: bitwise identical
  expect_identical(coef(f_hols), coef(f_ols))
  expect_identical(vcov(f_hols), vcov(f_ols))
  expect_identical(f_hols$sigma, f_ols$sigma)
  expect_equal(f_hols$method, "ols")
  expect_equal(nobs(f_hols), 428L)

  # But the instrument set is real: K1 = 0, L1 = 3
  expect_equal(length(f_hols$endogenous), 0L)
  expect_equal(f_hols$instruments, c("age", "kidslt6", "kidsge6"))
})

test_that("empty-endogenous fit is exact OLS (robust and cluster)", {
  f_hols <- ivreg2(hols_fml, data = mroz, vcov = "robust")
  f_ols  <- ivreg2(ols_fml, data = mroz, vcov = "robust")
  expect_identical(coef(f_hols), coef(f_ols))
  expect_identical(vcov(f_hols), vcov(f_ols))

  f_hols_cl <- ivreg2(hols_fml, data = mroz, clusters = ~hushrs)
  f_ols_cl  <- ivreg2(ols_fml, data = mroz, clusters = ~hushrs)
  expect_identical(vcov(f_hols_cl), vcov(f_ols_cl))
})

test_that("HOLS via gmm2s differs from OLS and is efficient GMM", {
  f_hols <- suppressWarnings(
    ivreg2(hols_fml, data = mroz, method = "gmm2s", vcov = "robust")
  )
  f_ols <- ivreg2(ols_fml, data = mroz, vcov = "robust")
  expect_equal(f_hols$method, "gmm2s")
  # The surplus moment conditions shift the point estimates
  expect_gt(max(abs(coef(f_hols) - coef(f_ols)[names(coef(f_hols))])), 1e-6)
  expect_true(is.finite(f_hols$diagnostics$overid$stat))
})


# ============================================================================
# Diagnostics gating (Stata master gate: endo1_ct > 0, ivreg2.ado:1625)
# ============================================================================

test_that("K1 = 0: overid computed with df = L1; endo diagnostics skipped", {
  fit <- ivreg2(hols_fml, data = mroz)

  expect_equal(fit$diagnostics$overid$test_name, "Sargan")
  expect_equal(fit$diagnostics$overid$df, 3L)
  expect_true(is.finite(fit$diagnostics$overid$stat))

  # Identification / endogenous-regressor diagnostics: all skipped
  expect_null(fit$diagnostics$underid)
  expect_null(fit$diagnostics$weak_id)
  expect_null(fit$diagnostics$weak_id_robust)
  expect_null(fit$diagnostics$weak_id_sy)
  expect_null(fit$diagnostics$anderson_rubin)
  expect_null(fit$diagnostics$stock_wright)
  expect_null(fit$diagnostics$endogeneity)
  expect_null(fit$diagnostics$redundancy)
  expect_null(fit$first_stage)
})

test_that("orthog() of an included exogenous regressor works with K1 = 0", {
  # Help-file H34: C-test of educ exogeneity in the OLS model
  fit <- ivreg2(hols_fml, data = mroz, orthog = "educ")
  expect_true(is.finite(fit$diagnostics$orthog$stat))
  expect_equal(fit$diagnostics$orthog$df, 1L)
})

test_that("stored S and W cover the full instrument set", {
  fit <- ivreg2(hols_fml, data = mroz)
  # L = 4 included (incl intercept) + 3 excluded = 7
  expect_equal(dim(fit$S), c(7L, 7L))
  expect_equal(dim(fit$W), c(7L, 7L))
  expect_true(all(c("age", "kidslt6", "kidsge6") %in% colnames(fit$S)))
})


# ============================================================================
# J-equivalence identities (help-file H48 / H52 demos)
# ============================================================================

test_that("J of the empty-endo first-stage equation equals KP idstat (H48)", {
  f_iv <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                 data = mroz_lw, vcov = "robust")
  f_j <- ivreg2(educ ~ exper + expersq | 0 | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust")
  expect_equal(f_j$diagnostics$overid$stat, f_iv$diagnostics$underid$stat,
               tolerance = 1e-8)
})

test_that("J of the empty-endo equation equals the redundancy LM stat (H52)", {
  f_iv <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                 data = mroz_lw, vcov = "robust", redundant = "age")
  f_j <- ivreg2(educ ~ exper + expersq + kidslt6 + kidsge6 | 0 | age,
                data = mroz_lw, vcov = "robust")
  expect_equal(f_j$diagnostics$overid$stat, f_iv$diagnostics$redundancy$stat,
               tolerance = 1e-8)
})


# ============================================================================
# Guards
# ============================================================================

test_that("liml / kclass error with no endogenous regressors", {
  expect_error(
    ivreg2(hols_fml, data = mroz, method = "liml"),
    "requires at least one endogenous regressor"
  )
  expect_error(
    ivreg2(hols_fml, data = mroz, kclass = 0.5),
    "requires at least one endogenous regressor"
  )
  expect_error(
    ivreg2(hols_fml, data = mroz, fuller = 1),
    "requires at least one endogenous regressor"
  )
})

test_that("endog() errors with no endogenous regressors (Stata :519-523)", {
  expect_error(
    ivreg2(hols_fml, data = mroz, endog = "educ"),
    "not in the endogenous list"
  )
})

test_that("redundant() warns and is skipped with K1 = 0", {
  expect_warning(
    fit <- ivreg2(hols_fml, data = mroz, redundant = "age"),
    "no endogenous regressors"
  )
  expect_null(fit$diagnostics$redundancy)
})

test_that("both parts empty errors; reduced_form silently skipped", {
  expect_error(
    ivreg2(lwage ~ exper | 0 | 0, data = mroz),
    "both empty"
  )
  fit <- ivreg2(hols_fml, data = mroz, reduced_form = "rf")
  expect_null(fit$reduced_form)
})

test_that("all excluded instruments collinear degrades to plain OLS", {
  d <- mroz
  d$exper2 <- 2 * d$exper   # collinear with the included exper
  expect_warning(
    expect_warning(
      fit <- ivreg2(lwage ~ exper + expersq | 0 | exper2, data = d),
      "Dropped 1 collinear instrument"
    ),
    "model is now plain OLS"
  )
  expect_null(fit$diagnostics)
  expect_equal(length(fit$instruments), 0L)
  expect_equal(fit$method, "ols")
})


# ============================================================================
# Reclassification parity: all endo reclassified keeps Z (latent F3 fix)
# ============================================================================

test_that("all-endo-reclassified model keeps instruments and posts J", {
  d <- mroz
  d$exper_c <- 3 * d$exper + 1   # endogenous slot, exactly collinear with exper
  ws <- character()
  fit <- withCallingHandlers(
    ivreg2(lwage ~ exper | exper_c | age + kidslt6, data = d),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("Now treated as exogenous: exper_c", ws)))
  expect_true(any(grepl("Dropped 1 collinear regressor: exper_c", ws)))
  expect_equal(fit$method, "ols")
  expect_equal(length(fit$endogenous), 0L)
  expect_equal(fit$instruments, c("age", "kidslt6"))
  expect_equal(fit$diagnostics$overid$df, 2L)
  expect_true(is.finite(fit$diagnostics$overid$stat))
})


# ============================================================================
# S3 methods smoke tests
# ============================================================================

test_that("summary/print/glance/predict work on a K1 = 0 fit", {
  fit <- ivreg2(hols_fml, data = mroz)

  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("OLS Estimation", out)))
  expect_true(any(grepl("Excluded instruments:.*age", out)))
  expect_true(any(grepl("Overidentification test \\(Sargan\\)", out)))
  expect_false(any(grepl("Instrumented:", out)))
  expect_false(any(grepl("Stock-Yogo", out)))

  gl <- glance(fit)
  expect_equal(gl$method, "ols")

  pr <- predict(fit, newdata = mroz[1:5, ])
  expect_equal(length(pr), 5L)

  Z <- model.matrix(fit, "instruments")
  expect_equal(ncol(Z), 7L)
  Xp <- model.matrix(fit, "projected")
  expect_equal(ncol(Xp), 4L)
})
