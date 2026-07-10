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

test_that("HOLS via gmm2s + robust differs from OLS; bare gmm2s equals OLS", {
  # With robust VCE the surplus moment conditions buy efficiency and shift
  # the point estimates (Cragg's 1983 HOLS) ...
  f_hols <- ivreg2(hols_fml, data = mroz, method = "gmm2s", vcov = "robust")
  f_ols <- ivreg2(ols_fml, data = mroz, vcov = "robust")
  expect_equal(f_hols$method, "gmm2s")
  expect_gt(max(abs(coef(f_hols) - coef(f_ols)[names(coef(f_hols))])), 1e-6)
  expect_true(is.finite(f_hols$diagnostics$overid$stat))

  # ... but with the default iid VCE the 2-step weighting matrix is
  # proportional to (Z'Z)^-1, so gmm2s = 2SLS = OLS (Stata help.txt:344-358:
  # bare gmm2s is "IV/2SLS, SEs consistent under homoskedasticity")
  f_iid <- ivreg2(hols_fml, data = mroz, method = "gmm2s")
  expect_equal(coef(f_iid), coef(ivreg2(ols_fml, data = mroz)),
               tolerance = 1e-10)
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

test_that("liml / kclass with K1 = 0 equal OLS exactly (any k)", {
  # X is contained in Z, so X'M_Z = 0 and the k-class normal equations
  # collapse to X'X for every k: liml, Fuller, and arbitrary kclass all
  # reproduce OLS. Stata runs this case (e(model)="liml"); verified.
  f_ols <- ivreg2(ols_fml, data = mroz)
  for (fit in list(ivreg2(hols_fml, data = mroz, method = "liml"),
                   ivreg2(hols_fml, data = mroz, fuller = 1),
                   ivreg2(hols_fml, data = mroz, kclass = 0.5))) {
    expect_equal(coef(fit), coef(f_ols), tolerance = 1e-10)
    expect_equal(unname(vcov(fit)), unname(vcov(f_ols)), tolerance = 1e-10)
  }
  f_liml <- ivreg2(hols_fml, data = mroz, method = "liml")
  # lambda = RSS(y ~ X) / RSS(y ~ Z) >= 1, the 1x1 concentration eigenvalue
  expect_gt(f_liml$lambda, 1)
  # Anderson-Rubin LR overid posted (df = L1)
  expect_equal(f_liml$diagnostics$anderson_rubin_overid$df, 3L)
  expect_true(is.finite(f_liml$diagnostics$anderson_rubin_overid$lr_stat))
})

test_that("null model (no regressors at all) errors like Stata", {
  # Stata's CheckMisc: "no regressors specified" (ivreg2.ado:4246-4249)
  expect_error(ivreg2(lwage ~ 0 | 0 | age, data = mroz),
               "No regressors specified")
  expect_error(ivreg2(lwage ~ 0, data = mroz),
               "No regressors specified")
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

test_that("both parts empty errors; reduced_form skipped with a warning", {
  expect_error(
    ivreg2(lwage ~ exper | 0 | 0, data = mroz),
    "both empty"
  )
  # No reduced form is computed or stored (matching Stata, which skips the
  # whole RF block when the endogenous list is empty), but the ignored
  # explicit request is warned about rather than dropped silently.
  expect_warning(
    fit <- ivreg2(hols_fml, data = mroz, reduced_form = "rf"),
    "`reduced_form` ignored: the model has no endogenous regressors",
    fixed = TRUE
  )
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
# Stata fixture comparisons (generate-hols-fixtures.do)
# ============================================================================

# Compare a HOLS fit against its Stata fixture set. e(S) is compared only
# for the canonical-base configs (compare_s = TRUE), which share the
# mroz_hols_inames.csv instrument set; the H48/H52 companion fits have
# different instrument lists.
check_hols_fixture <- function(fit, suffix, compare_s = FALSE) {
  coef_path <- fixture_path(paste0("mroz_hols_coef_", suffix, ".csv"))
  diag_path <- fixture_path(paste0("mroz_hols_diagnostics_", suffix, ".csv"))
  skip_if(!file.exists(coef_path), "Fixture not found")

  stata_coef <- read.csv(coef_path, stringsAsFactors = FALSE)
  stata_coef$term[stata_coef$term == "_cons"] <- "(Intercept)"
  r_names <- names(coef(fit))
  m <- match(r_names, stata_coef$term)
  expect_equal(unname(coef(fit)), stata_coef$estimate[m],
               tolerance = stata_tol$coef, info = paste(suffix, "coef"))
  expect_equal(unname(sqrt(diag(fit$vcov))), stata_coef$std_error[m],
               tolerance = stata_tol$se, info = paste(suffix, "se"))

  vcov_path <- fixture_path(paste0("mroz_hols_vcov_", suffix, ".csv"))
  V_stata <- as.matrix(read.csv(vcov_path))
  stata_to_r <- match(stata_coef$term, r_names)
  expect_vcov_match(unname(fit$vcov[stata_to_r, stata_to_r]), V_stata,
                    label = suffix)

  if (compare_s) {
    S_stata <- read_stata_matrix(
      fixture_path(paste0("mroz_hols_eS_", suffix, ".csv")),
      fixture_path("mroz_hols_inames.csv")
    )
    expect_vcov_match(fit$S[rownames(S_stata), colnames(S_stata)], S_stata,
                      label = paste(suffix, "e(S)"))
  }

  stata_diag <- read.csv(diag_path, na.strings = c("", "."))
  diag <- fit$diagnostics
  expect_identical(fit$nobs, as.integer(stata_diag$N))

  # Stata posts NO idstat/widstat for (=z) fits — and neither do we
  if (is.na(stata_diag$underid_stat)) {
    expect_null(diag$underid)
  } else {
    expect_equal(diag$underid$stat, stata_diag$underid_stat,
                 tolerance = stata_tol$stat, info = paste(suffix, "idstat"))
  }

  if (!is.na(stata_diag$overid_stat)) {
    expect_equal(diag$overid$stat, stata_diag$overid_stat,
                 tolerance = stata_tol$stat, info = paste(suffix, "J"))
    expect_equal(diag$overid$p, stata_diag$overid_p,
                 tolerance = stata_tol$pval, info = paste(suffix, "J p"))
    expect_identical(diag$overid$df, as.integer(stata_diag$overid_df))
  }
  if (!is.na(stata_diag$cstat) && !is.null(diag$orthog)) {
    expect_equal(diag$orthog$stat, stata_diag$cstat,
                 tolerance = stata_tol$stat, info = paste(suffix, "C-stat"))
    expect_identical(diag$orthog$df, as.integer(stata_diag$cstat_df))
  }
  if (!is.na(stata_diag$redstat) && !is.null(diag$redundancy)) {
    expect_equal(diag$redundancy$stat, stata_diag$redstat,
                 tolerance = stata_tol$stat, info = paste(suffix, "redstat"))
  }
  if (!is.na(stata_diag$model_f)) {
    expect_equal(fit$model_f, stata_diag$model_f,
                 tolerance = stata_tol$stat, info = paste(suffix, "model F"))
  }
}

test_that("HOLS base iid matches Stata (incl. Sargan and e(S))", {
  skip_if(!file.exists(fixture_path("mroz_hols_coef_iid.csv")),
          "Fixture not found")
  check_hols_fixture(ivreg2(hols_fml, data = mroz), "iid", compare_s = TRUE)
})

test_that("HOLS base robust matches Stata (incl. Hansen J)", {
  skip_if(!file.exists(fixture_path("mroz_hols_coef_robust.csv")),
          "Fixture not found")
  check_hols_fixture(ivreg2(hols_fml, data = mroz, vcov = "robust"),
                     "robust", compare_s = TRUE)
})

test_that("H34: orthog(educ) C-test matches Stata", {
  skip_if(!file.exists(fixture_path("mroz_hols_coef_orthog.csv")),
          "Fixture not found")
  check_hols_fixture(ivreg2(hols_fml, data = mroz, orthog = "educ"),
                     "orthog", compare_s = TRUE)
})

test_that("H35 exact command (bare gmm2s) matches Stata and equals OLS", {
  skip_if(!file.exists(fixture_path("mroz_hols_coef_gmm2s_iid.csv")),
          "Fixture not found")
  fit <- ivreg2(hols_fml, data = mroz, method = "gmm2s")
  check_hols_fixture(fit, "gmm2s_iid", compare_s = TRUE)
})

test_that("HOLS proper (gmm2s + robust) matches Stata", {
  skip_if(!file.exists(fixture_path("mroz_hols_coef_gmm2s.csv")),
          "Fixture not found")
  fit <- ivreg2(hols_fml, data = mroz, method = "gmm2s", vcov = "robust")
  check_hols_fixture(fit, "gmm2s", compare_s = TRUE)
})

test_that("intercept-only partial = \"_all\" is a no-op, matching Stata", {
  # Stata's `_all` expands to the included exogenous regressors EXCLUDING
  # the constant (ivreg2.ado:3708-3710); for y ~ 1 the list is empty, the
  # partial block is skipped, and the intercept fit proceeds (live-verified
  # rc=0, e(partial_ct)=0). Locks the parity so a future "fix" doesn't
  # break it.
  fit <- ivreg2(lwage ~ 1, data = mroz, partial = "_all")
  expect_equal(length(coef(fit)), 1L)
  expect_equal(fit$partial_ct %||% 0L, 0L)
})

test_that("partial = \"_all\" with no surviving regressors errors like Stata", {
  # Stata: CheckMisc recounts rhs1_ct AFTER partialling (ivreg2.ado:633/641,
  # 4246-4249) and exits 102. Reachable for empty-endogenous and 1-part OLS
  # models, whose regressors are all exogenous.
  expect_error(
    ivreg2(hols_fml, data = mroz, partial = "_all"),
    "No regressors remain after partialling"
  )
  expect_error(
    ivreg2(lwage ~ exper, data = mroz, partial = "_all"),
    "No regressors remain after partialling"
  )
})

test_that("LIML with K1 = 0 matches Stata (coef = OLS, lambda, arubin)", {
  skip_if(!file.exists(fixture_path("mroz_hols_coef_liml.csv")),
          "Fixture not found")
  fit <- ivreg2(hols_fml, data = mroz, method = "liml")
  check_hols_fixture(fit, "liml", compare_s = TRUE)
  stata_diag <- read.csv(fixture_path("mroz_hols_diagnostics_liml.csv"),
                         na.strings = c("", "."))
  expect_equal(fit$lambda, stata_diag$lambda, tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$anderson_rubin_overid$lr_stat,
               stata_diag$arubin, tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$anderson_rubin_overid$lr_p,
               stata_diag$arubin_p, tolerance = stata_tol$pval)
})

test_that("H48: companion IV fit and J-equation both match Stata", {
  skip_if(!file.exists(fixture_path("mroz_hols_coef_h48_iv.csv")),
          "Fixture not found")
  f_iv <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                 data = mroz_lw, vcov = "robust")
  check_hols_fixture(f_iv, "h48_iv")
  f_j <- ivreg2(educ ~ exper + expersq | 0 | age + kidslt6 + kidsge6,
                data = mroz_lw, vcov = "robust")
  check_hols_fixture(f_j, "h48_j")
})

test_that("H52: companion redundant(age) fit and J-equation match Stata", {
  skip_if(!file.exists(fixture_path("mroz_hols_coef_h52_iv.csv")),
          "Fixture not found")
  f_iv <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                 data = mroz_lw, vcov = "robust", redundant = "age")
  check_hols_fixture(f_iv, "h52_iv")
  f_j <- ivreg2(educ ~ exper + expersq + kidslt6 + kidsge6 | 0 | age,
                data = mroz_lw, vcov = "robust")
  check_hols_fixture(f_j, "h52_j")
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

  expect_equal(fit$method, "ols")

  pr <- predict(fit, newdata = mroz[1:5, ])
  expect_equal(length(pr), 5L)

  Z <- model.matrix(fit, "instruments")
  expect_equal(ncol(Z), 7L)
  Xp <- model.matrix(fit, "projected")
  expect_equal(ncol(Xp), 4L)
})
