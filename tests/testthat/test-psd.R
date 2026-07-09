# ============================================================================
# Tests: PSD Corrections — Ticket P2
# ============================================================================

# --- Helpers ---
data(card)

# ============================================================================
# Unit tests for .psd_correct()
# ============================================================================

test_that(".psd_correct returns input unchanged when psd is NULL", {
  mat <- matrix(c(4, 2, 2, 3), 2, 2)
  expect_identical(.psd_correct(mat, NULL), mat)
})

test_that(".psd_correct returns PSD matrix unchanged", {
  mat <- matrix(c(4, 2, 2, 3), 2, 2)  # eigenvalues: 5, 2 (both positive)
  expect_equal(.psd_correct(mat, "psd0"), mat)
  expect_equal(.psd_correct(mat, "psda"), mat)
})

test_that(".psd_correct with psd0 zeros negative eigenvalues", {
  # Construct matrix with known eigenvalues: one negative
  # V = [1/sqrt(2), 1/sqrt(2); 1/sqrt(2), -1/sqrt(2)]
  # D = diag(3, -1)
  # A = V %*% D %*% t(V)
  V <- matrix(c(1, 1, 1, -1) / sqrt(2), 2, 2)
  D <- diag(c(3, -1))
  mat <- V %*% D %*% t(V)

  result <- suppressWarnings(.psd_correct(mat, "psd0"))

  # Expected: eigenvalue -1 → 0, so D_corrected = diag(3, 0)
  expected <- V %*% diag(c(3, 0)) %*% t(V)
  expected <- (expected + t(expected)) / 2
  expect_equal(result, expected)

  # Verify PSD: all eigenvalues >= 0
  eig <- eigen(result, symmetric = TRUE)$values
  expect_true(all(eig >= 0))
})

test_that(".psd_correct with psda takes absolute value of negative eigenvalues", {
  V <- matrix(c(1, 1, 1, -1) / sqrt(2), 2, 2)
  D <- diag(c(3, -1))
  mat <- V %*% D %*% t(V)

  result <- suppressWarnings(.psd_correct(mat, "psda"))

  # Expected: eigenvalue -1 → 1, so D_corrected = diag(3, 1)
  expected <- V %*% diag(c(3, 1)) %*% t(V)
  expected <- (expected + t(expected)) / 2
  expect_equal(result, expected)

  # Verify PSD
  eig <- eigen(result, symmetric = TRUE)$values
  expect_true(all(eig >= 0))
})

test_that(".psd_correct handles multiple negative eigenvalues", {
  V <- diag(3)  # identity eigenvectors for simplicity
  D <- diag(c(5, -2, -0.5))
  mat <- V %*% D %*% t(V)

  result_psd0 <- suppressWarnings(.psd_correct(mat, "psd0"))
  result_psda <- suppressWarnings(.psd_correct(mat, "psda"))

  # psd0: negative → 0
  expect_equal(eigen(result_psd0, symmetric = TRUE)$values, c(5, 0, 0))
  # psda: negative → abs
  expect_equal(eigen(result_psda, symmetric = TRUE)$values, c(5, 2, 0.5))
})

test_that(".psd_correct result is symmetric", {
  V <- matrix(c(0.6, 0.8, 0.8, -0.6), 2, 2)
  D <- diag(c(2, -0.3))
  mat <- V %*% D %*% t(V)

  result <- suppressWarnings(.psd_correct(mat, "psd0"))
  expect_true(isSymmetric(result))

  result <- suppressWarnings(.psd_correct(mat, "psda"))
  expect_true(isSymmetric(result))
})

test_that(".psd_correct emits warning when correction is applied", {
  V <- matrix(c(1, 1, 1, -1) / sqrt(2), 2, 2)
  D <- diag(c(3, -1))
  mat <- V %*% D %*% t(V)

  expect_warning(.psd_correct(mat, "psd0"),
                 "1 negative eigenvalue corrected via the 'psd0' method")
  expect_warning(.psd_correct(mat, "psda"),
                 "1 negative eigenvalue corrected via the 'psda' method")
})

test_that(".psd_correct emits warning with correct plural for multiple eigenvalues", {
  D <- diag(c(5, -2, -0.5))
  mat <- D  # diagonal = eigenvalues

  expect_warning(.psd_correct(mat, "psd0"),
                 "2 negative eigenvalues corrected via the 'psd0' method")
})

test_that(".psd_correct emits no warning when no correction needed", {
  mat <- matrix(c(4, 2, 2, 3), 2, 2)
  expect_no_warning(.psd_correct(mat, "psd0"))
  expect_no_warning(.psd_correct(mat, "psda"))
})

test_that(".psd_correct ignores floating-point noise near zero", {
  # Construct a PSD matrix, then perturb one eigenvalue to be
 # very slightly negative (floating-point noise level)
  V <- matrix(c(1, 1, 1, -1) / sqrt(2), 2, 2)
  D <- diag(c(3, 1e-15))  # nearly zero but positive
  mat <- V %*% D %*% t(V)
  # Manually flip the tiny eigenvalue to a tiny negative
  D_neg <- diag(c(3, -1e-15))
  mat_neg <- V %*% D_neg %*% t(V)
  # Should NOT warn — the negative eigenvalue is floating-point noise
  expect_no_warning(.psd_correct(mat_neg, "psd0"))
  expect_no_warning(.psd_correct(mat_neg, "psda"))
})


# ============================================================================
# Integration tests with ivreg2()
# ============================================================================

test_that("psd parameter validation works", {
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc4, data = card,
           vcov = "robust", psd = "invalid"),
    "should be one of"
  )
})

test_that("psd = NULL is the default and changes nothing", {
  fit_default <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                         data = card, vcov = "robust")
  fit_null <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                      data = card, vcov = "robust", psd = NULL)
  expect_equal(coef(fit_default), coef(fit_null))
  expect_equal(vcov(fit_default), vcov(fit_null))
  expect_null(fit_default$psd)
  expect_null(fit_null$psd)
})

test_that("psd is stored in the returned object", {
  fit0 <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                  data = card, vcov = "robust", psd = "psd0")
  fita <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                  data = card, vcov = "robust", psd = "psda")
  expect_equal(fit0$psd, "psd0")
  expect_equal(fita$psd, "psda")
})

test_that("psd is stored on the fitted object", {
  fit0 <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                  data = card, vcov = "robust", psd = "psd0")
  expect_equal(fit0$psd, "psd0")

  fit_none <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                      data = card, vcov = "robust")
  expect_null(fit_none$psd)
})

test_that("psd appears in summary footer", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                 data = card, vcov = "robust", psd = "psda")
  out <- capture.output(print(summary(fit)))
  vcv_line <- grep("VCV type:", out, value = TRUE)
  expect_true(any(grepl("psda", vcv_line)))
})

test_that("non-pathological model: psd has no effect (VCV already PSD)", {
  # Card dataset with robust VCE: VCV should already be PSD
  fit_plain <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                       data = card, vcov = "robust")
  fit_psd0 <- suppressWarnings(
    ivreg2(lwage ~ exper + expersq | educ | nearc4,
            data = card, vcov = "robust", psd = "psd0")
  )
  fit_psda <- suppressWarnings(
    ivreg2(lwage ~ exper + expersq | educ | nearc4,
            data = card, vcov = "robust", psd = "psda")
  )
  # Coefficients unchanged
  expect_equal(coef(fit_plain), coef(fit_psd0))
  expect_equal(coef(fit_plain), coef(fit_psda))
  # VCV unchanged (no negative eigenvalues to correct)
  expect_equal(vcov(fit_plain), vcov(fit_psd0))
  expect_equal(vcov(fit_plain), vcov(fit_psda))
})

test_that("psd works with clustered VCE", {
  # Should not error; for well-conditioned data, results are the same
  # M=2 clusters → expected rank-deficient diagnostics
  fit <- muffle_rank_warnings(
    ivreg2(lwage ~ exper + expersq | educ | nearc4,
           data = card, clusters = ~smsa, psd = "psd0")
  )
  expect_equal(fit$psd, "psd0")
  expect_true(all(is.finite(coef(fit))))
})

test_that("psd works with iid VCE", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                 data = card, vcov = "iid", psd = "psda")
  expect_equal(fit$psd, "psda")
  expect_true(all(is.finite(coef(fit))))
})

test_that("psd works with GMM2S", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4 + nearc2,
                 data = card, method = "gmm2s", vcov = "robust", psd = "psd0")
  expect_equal(fit$psd, "psd0")
  expect_true(all(is.finite(coef(fit))))
})

test_that("psd works with OLS model", {
  # OLS with psd should not error (no-op for IID VCV)
  fit <- ivreg2(lwage ~ exper + expersq, data = card, psd = "psd0")
  expect_equal(fit$psd, "psd0")
  expect_true(all(is.finite(coef(fit))))
})


# ============================================================================
# F2: psd correction site — S-level correction, VCV assembled from S
# ============================================================================
# Stata applies psd0/psda to the L x L moment covariance S inside m_omega
# (livreg2.do:607-617) and assembles the plain non-GMM VCV from the corrected
# S by congruence with the first-stage map A (s_iegmm, ivreg2.ado:5556).
# The binding configuration below (Driscoll-Kraay + truncated kernel on the
# wagepan panel) produces an indefinite S, so the correction actually fires.

data(wagepan)

# Capture-all-warnings fitter. The binding DK config legitimately warns from
# several distinct sites (one psd correction per corrected matrix, plus the
# pre-existing singular Hansen J / Stock-Wright NA warnings), so tests assert
# on the captured set instead of nesting expect_warning() calls.
fit_capturing_warnings <- function(expr) {
  ws <- character()
  fit <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(fit = fit, warnings = ws)
}

dk_binding_fit <- function(psd = NULL, small = FALSE) {
  fit_capturing_warnings(ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wagepan, dkraay = 2, kernel = "truncated",
    tvar = "year", ivar = "nr", psd = psd, small = small
  ))
}

# Rebuild the plain-path VCV from fit$S: V = N * bread * (A' S A) * bread
# * small-factor. A and the bread are reconstructed from the model matrices
# (X_hat = Z A), independently of the internal code path.
vcov_from_stored_s <- function(fit, small_factor = 1) {
  Z <- model.matrix(fit, "instruments")
  X <- model.matrix(fit, "regressors")
  w <- fit$weights
  A <- if (is.null(w)) {
    qr.coef(qr(Z), X)
  } else {
    qr.coef(qr(sqrt(w) * Z), sqrt(w) * X)
  }
  A[is.na(A)] <- 0
  X_hat <- Z %*% A
  XtX <- if (is.null(w)) crossprod(X_hat) else crossprod(sqrt(w) * X_hat)
  bread <- chol2inv(chol(XtX))
  N <- nobs(fit)
  V <- N * (bread %*% crossprod(A, fit$S %*% A) %*% bread) * small_factor
  (V + t(V)) / 2
}

test_that("binding case: psd correction fires at the S level", {
  res_null <- dk_binding_fit(psd = NULL)
  res_psd0 <- dk_binding_fit(psd = "psd0")
  res_psda <- dk_binding_fit(psd = "psda")

  # S is indefinite without correction; corrected S is PSD
  eig_null <- eigen(res_null$fit$S, symmetric = TRUE)$values
  expect_true(min(eig_null) < 0)
  tol_eig <- .Machine$double.eps^0.5 * max(abs(eig_null))
  expect_true(min(eigen(res_psd0$fit$S, symmetric = TRUE)$values) >= -tol_eig)
  expect_true(min(eigen(res_psda$fit$S, symmetric = TRUE)$values) >= -tol_eig)

  # The correction warns; the uncorrected fit does not warn about psd
  expect_true(any(grepl("corrected via the 'psd0' method",res_psd0$warnings)))
  expect_true(any(grepl("corrected via the 'psda' method",res_psda$warnings)))
  expect_false(any(grepl("corrected via", res_null$warnings)))

  # Coefficients are untouched by psd on plain (non-GMM) paths
  expect_equal(coef(res_psd0$fit), coef(res_null$fit))
  expect_equal(coef(res_psda$fit), coef(res_null$fit))

  # The VCV changes when the correction binds, and psd0 != psda
  expect_gt(max(abs(vcov(res_psd0$fit) - vcov(res_null$fit))), 1e-4)
  expect_gt(max(abs(vcov(res_psd0$fit) - vcov(res_psda$fit))), 1e-6)
})

test_that("binding case: correction site is S, not the final VCV", {
  res_null <- dk_binding_fit(psd = NULL)
  res_psd0 <- dk_binding_fit(psd = "psd0")

  # Correcting the final VCV (the pre-F2 site) must NOT reproduce the
  # S-level correction: B psd(S) B' != psd(B S B').
  v_oldsite <- suppressWarnings(.psd_correct(vcov(res_null$fit), "psd0"))
  expect_gt(max(abs(v_oldsite - vcov(res_psd0$fit))), 1e-6)
})

test_that("plain-path VCV equals the sandwich rebuilt from fit$S", {
  # Binding DK case (cluster family, small = FALSE: factor 1)
  res_psd0 <- dk_binding_fit(psd = "psd0")
  v_rebuilt <- vcov_from_stored_s(res_psd0$fit)
  expect_equal(unname(vcov(res_psd0$fit)), v_rebuilt, tolerance = 1e-8)

  # Same with psda
  res_psda <- dk_binding_fit(psd = "psda")
  v_rebuilt_a <- vcov_from_stored_s(res_psda$fit)
  expect_equal(unname(vcov(res_psda$fit)), v_rebuilt_a, tolerance = 1e-8)

  # Binding DK case with small = TRUE: cluster-family factor
  res_small <- dk_binding_fit(psd = "psd0", small = TRUE)
  fit_s <- res_small$fit
  N <- nobs(fit_s)
  K <- length(coef(fit_s))
  M <- fit_s$n_clusters
  factor_cl <- ((N - 1) / (N - K)) * (M / (M - 1))
  v_rebuilt_s <- vcov_from_stored_s(fit_s, small_factor = factor_cl)
  expect_equal(unname(vcov(fit_s)), v_rebuilt_s, tolerance = 1e-8)
})

test_that("non-binding robust path: psd inert; VCV-from-S identity holds", {
  fit_null <- ivreg2(lwage ~ exper + expersq | educ | nearc4 + nearc2,
                     data = card, vcov = "robust")
  fit_psd <- ivreg2(lwage ~ exper + expersq | educ | nearc4 + nearc2,
                    data = card, vcov = "robust", psd = "psd0")

  # Same numbers up to assembly-route floating-point noise
  expect_equal(vcov(fit_null), vcov(fit_psd), tolerance = 1e-9)

  # V = N * bread * A' S A * bread (robust family, small = FALSE)
  v_rebuilt <- vcov_from_stored_s(fit_psd)
  expect_equal(unname(vcov(fit_psd)), v_rebuilt, tolerance = 1e-8)
})

test_that("Stock-Watson path: psd routes through the L x L SW omega", {
  fit_sw <- muffle_rank_warnings(
    ivreg2(lwage ~ exper + expersq + married + union | hours | educ + black,
           data = wagepan, sw = TRUE, ivar = "nr", psd = "psd0")
  )
  # SW omega is per-observation scale; same uniform identity applies
  v_rebuilt <- vcov_from_stored_s(fit_sw)
  expect_equal(unname(vcov(fit_sw)), v_rebuilt, tolerance = 1e-8)
})

test_that("KP identification stats are never psd-corrected (ranktest parity)", {
  # Stata's ranktest does not receive the psd option (ivreg2.ado:1639-1650),
  # so underid/weak-id statistics must be identical across psd settings even
  # when the VCV correction binds.
  res_null <- dk_binding_fit(psd = NULL)
  res_psd0 <- dk_binding_fit(psd = "psd0")
  res_psda <- dk_binding_fit(psd = "psda")

  expect_identical(res_psd0$fit$diagnostics$underid$stat,
                   res_null$fit$diagnostics$underid$stat)
  expect_identical(res_psda$fit$diagnostics$underid$stat,
                   res_null$fit$diagnostics$underid$stat)
  expect_identical(res_psd0$fit$diagnostics$weak_id_robust$stat,
                   res_null$fit$diagnostics$weak_id_robust$stat)
  expect_identical(res_psda$fit$diagnostics$weak_id_robust$stat,
                   res_null$fit$diagnostics$weak_id_robust$stat)
})

test_that("iid VCV is never psd-corrected (Stata m_omega parity)", {
  fit_null <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                     data = card, vcov = "iid")
  fit_psd <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                    data = card, vcov = "iid", psd = "psd0")
  # Identical code path: psd is inert for the iid VCV
  expect_identical(vcov(fit_null), vcov(fit_psd))
})


# ============================================================================
# F2: Stata fixture comparisons (generate-psd-fixtures.do)
# ============================================================================

# Compare a psd fit against its Stata fixture set: coefficients, SEs, full
# VCV, e(S), and the diagnostics Stata reports. `skip_model_f = TRUE` for
# the psd0-binding configs, where the corrected VCV is exactly singular and
# the joint model F is not a well-defined Wald statistic: Stata inverts the
# near-singular slope block raw (F ~ 2.9e4 on the binding repro), while our
# conditioning guard (.is_badly_conditioned_vcov) falls back to the sweep
# inverse (F ~ 7.7). Documented intentional divergence; see
# planning/06-parity-deltas.md row 18.
check_psd_fixture <- function(fit, suffix, prefix = "wp_psd",
                              inames_file = "wp_psd_inames.csv",
                              skip_model_f = FALSE) {
  coef_path <- fixture_path(paste0(prefix, "_coef_", suffix, ".csv"))
  vcov_path <- fixture_path(paste0(prefix, "_vcov_", suffix, ".csv"))
  s_path    <- fixture_path(paste0(prefix, "_eS_", suffix, ".csv"))
  diag_path <- fixture_path(paste0(prefix, "_diagnostics_", suffix, ".csv"))
  skip_if(!file.exists(coef_path), "Fixture not found")

  stata_coef <- read.csv(coef_path, stringsAsFactors = FALSE)
  stata_coef$term[stata_coef$term == "_cons"] <- "(Intercept)"
  r_names <- names(coef(fit))

  m <- match(r_names, stata_coef$term)
  expect_equal(unname(coef(fit)), stata_coef$estimate[m],
               tolerance = stata_tol$coef,
               info = paste(suffix, "coefficients"))
  expect_equal(unname(sqrt(diag(fit$vcov))), stata_coef$std_error[m],
               tolerance = stata_tol$se,
               info = paste(suffix, "std errors"))

  # Full VCV, reordered to Stata's coefficient order
  V_stata <- as.matrix(read.csv(vcov_path))
  stata_to_r <- match(stata_coef$term, r_names)
  expect_vcov_match(unname(fit$vcov[stata_to_r, stata_to_r]), V_stata,
                    label = suffix)

  # e(S), aligned by instrument name
  S_stata <- read_stata_matrix(s_path, fixture_path(inames_file))
  expect_vcov_match(fit$S[rownames(S_stata), colnames(S_stata)], S_stata,
                    label = paste(suffix, "e(S)"))

  # Diagnostics Stata reports
  stata_diag <- read.csv(diag_path, na.strings = c("", "."))
  diag <- fit$diagnostics

  if (is.na(stata_diag$overid_stat)) {
    # Stata suppressed J (rank-deficient corrected S); we must return NA too
    expect_true(is.null(diag$overid) || is.na(diag$overid$stat),
                label = paste(suffix, "J suppressed like Stata"))
  } else {
    expect_equal(diag$overid$stat, stata_diag$overid_stat,
                 tolerance = stata_tol$stat, info = paste(suffix, "J"))
    expect_equal(diag$overid$p, stata_diag$overid_p,
                 tolerance = stata_tol$pval, info = paste(suffix, "J p"))
  }

  # KP identification stats: never psd-corrected, present in all fixtures
  expect_equal(diag$underid$stat, stata_diag$underid_stat,
               tolerance = stata_tol$stat, info = paste(suffix, "KP LM"))
  expect_equal(diag$weak_id_robust$stat, stata_diag$weak_id_kp_f,
               tolerance = stata_tol$stat, info = paste(suffix, "KP F"))
  expect_equal(diag$weak_id$stat, stata_diag$weak_id_cd_f,
               tolerance = stata_tol$stat, info = paste(suffix, "CD F"))

  if (!skip_model_f && !is.na(stata_diag$model_f)) {
    expect_equal(fit$model_f, stata_diag$model_f,
                 tolerance = stata_tol$stat, info = paste(suffix, "model F"))
  }

  expect_identical(fit$nobs, as.integer(stata_diag$N))
}

test_that("binding DK + truncated psd0 matches Stata", {
  res <- dk_binding_fit(psd = "psd0")
  expect_true(any(grepl("corrected via the 'psd0' method",res$warnings)))
  check_psd_fixture(res$fit, "dk_tru2_psd0", skip_model_f = TRUE)
})

test_that("binding DK + truncated psda matches Stata", {
  res <- dk_binding_fit(psd = "psda")
  expect_true(any(grepl("corrected via the 'psda' method",res$warnings)))
  check_psd_fixture(res$fit, "dk_tru2_psda")
})

test_that("binding DK + truncated psd0 + small matches Stata", {
  res <- dk_binding_fit(psd = "psd0", small = TRUE)
  check_psd_fixture(res$fit, "dk_tru2_psd0_small", skip_model_f = TRUE)
})

test_that("Stock-Watson + psda matches Stata", {
  fit <- ivreg2(lwage ~ exper + expersq + married + union | hours |
                  educ + black,
                data = wagepan, sw = TRUE, ivar = "nr", psd = "psda")
  check_psd_fixture(fit, "sw_psda")
})

test_that("two-way cluster + psd0 matches Stata", {
  fit <- ivreg2(lwage ~ exper + expersq + married + union | hours |
                  educ + black,
                data = wagepan, clusters = ~nr + year, psd = "psd0")
  check_psd_fixture(fit, "twoway_psd0")
})

test_that("all-zero Omega returns NA J instead of crashing", {
  # Exactly-zero residuals (y identically 0) give an all-zero moment
  # covariance; the rank check must catch it (NA + warning, matching
  # Stata's "not of full rank" suppression) rather than crashing in
  # qr.solve(). Regression test for the <= boundary in the eigenvalue
  # rank check (.compute_j_with_omega).
  set.seed(42)
  d <- data.frame(y0 = 0, x = rnorm(50), z1 = rnorm(50),
                  z2 = rnorm(50), z3 = rnorm(50))
  # The zero-residual fit degrades every moment-based diagnostic, so
  # AR and Stock-Wright also warn; capture the set and assert on it.
  res <- fit_capturing_warnings(
    ivreg2(y0 ~ 1 | x | z1 + z2 + z3, data = d, vcov = "robust")
  )
  expect_true(any(grepl("Hansen J statistic not computed", res$warnings)))
  expect_true(all(grepl(
    "Hansen J statistic not computed|Anderson-Rubin.*could not be computed|Stock-Wright: Omega is rank-deficient",
    res$warnings
  )))
  expect_true(is.na(res$fit$diagnostics$overid$stat))
})

test_that("HAC truncated + psd0 matches Stata", {
  data(phillips)
  fit <- ivreg2(cinf ~ 1 | unem | unem_1 + unem_2, data = phillips,
                vcov = "robust", kernel = "truncated", bw = 2,
                tvar = "year", psd = "psd0")
  check_psd_fixture(fit, "hac_tru2_psd0", prefix = "phil_psd",
                    inames_file = "phil_psd_inames.csv")
})
