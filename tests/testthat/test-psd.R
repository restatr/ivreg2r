# ============================================================================
# Tests: PSD Corrections — Ticket P2
# ============================================================================

# --- Helpers ---
card_path <- fixture_path("card_data.csv")

if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

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
                 "1 negative eigenvalue corrected via psd0")
  expect_warning(.psd_correct(mat, "psda"),
                 "1 negative eigenvalue corrected via psda")
})

test_that(".psd_correct emits warning with correct plural for multiple eigenvalues", {
  D <- diag(c(5, -2, -0.5))
  mat <- D  # diagonal = eigenvalues

  expect_warning(.psd_correct(mat, "psd0"),
                 "2 negative eigenvalues corrected via psd0")
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
  skip_if_not(file.exists(card_path))
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc4, data = card,
           vcov = "robust", psd = "invalid"),
    "should be one of"
  )
})

test_that("psd = NULL is the default and changes nothing", {
  skip_if_not(file.exists(card_path))
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
  skip_if_not(file.exists(card_path))
  fit0 <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                  data = card, vcov = "robust", psd = "psd0")
  fita <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                  data = card, vcov = "robust", psd = "psda")
  expect_equal(fit0$psd, "psd0")
  expect_equal(fita$psd, "psda")
})

test_that("psd appears in glance() output", {
  skip_if_not(file.exists(card_path))
  fit0 <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                  data = card, vcov = "robust", psd = "psd0")
  gl <- glance(fit0)
  expect_true("psd" %in% names(gl))
  expect_equal(gl$psd, "psd0")

  fit_none <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                      data = card, vcov = "robust")
  gl_none <- glance(fit_none)
  expect_true(is.na(gl_none$psd))
})

test_that("psd appears in summary footer", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                 data = card, vcov = "robust", psd = "psda")
  out <- capture.output(print(summary(fit)))
  vcv_line <- grep("VCV type:", out, value = TRUE)
  expect_true(any(grepl("psda", vcv_line)))
})

test_that("non-pathological model: psd has no effect (VCV already PSD)", {
  skip_if_not(file.exists(card_path))
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
  skip_if_not(file.exists(card_path))
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
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                 data = card, vcov = "iid", psd = "psda")
  expect_equal(fit$psd, "psda")
  expect_true(all(is.finite(coef(fit))))
})

test_that("psd works with GMM2S", {
  skip_if_not(file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4 + nearc2,
                 data = card, method = "gmm2s", vcov = "robust", psd = "psd0")
  expect_equal(fit$psd, "psd0")
  expect_true(all(is.finite(coef(fit))))
})

test_that("psd works with OLS model", {
  skip_if_not(file.exists(card_path))
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

wp_psd_path <- fixture_path("wp_ck_data.csv")
if (file.exists(wp_psd_path)) {
  wp_psd <- read.csv(wp_psd_path)
}

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
    data = wp_psd, dkraay = 2, kernel = "truncated",
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
  skip_if_not(file.exists(wp_psd_path))

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
  expect_true(any(grepl("corrected via psd0", res_psd0$warnings)))
  expect_true(any(grepl("corrected via psda", res_psda$warnings)))
  expect_false(any(grepl("corrected via", res_null$warnings)))

  # Coefficients are untouched by psd on plain (non-GMM) paths
  expect_equal(coef(res_psd0$fit), coef(res_null$fit))
  expect_equal(coef(res_psda$fit), coef(res_null$fit))

  # The VCV changes when the correction binds, and psd0 != psda
  expect_gt(max(abs(vcov(res_psd0$fit) - vcov(res_null$fit))), 1e-4)
  expect_gt(max(abs(vcov(res_psd0$fit) - vcov(res_psda$fit))), 1e-6)
})

test_that("binding case: correction site is S, not the final VCV", {
  skip_if_not(file.exists(wp_psd_path))

  res_null <- dk_binding_fit(psd = NULL)
  res_psd0 <- dk_binding_fit(psd = "psd0")

  # Correcting the final VCV (the pre-F2 site) must NOT reproduce the
  # S-level correction: B psd(S) B' != psd(B S B').
  v_oldsite <- suppressWarnings(.psd_correct(vcov(res_null$fit), "psd0"))
  expect_gt(max(abs(v_oldsite - vcov(res_psd0$fit))), 1e-6)
})

test_that("plain-path VCV equals the sandwich rebuilt from fit$S", {
  skip_if_not(file.exists(wp_psd_path))

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
  skip_if_not(file.exists(card_path))

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
  skip_if_not(file.exists(wp_psd_path))

  fit_sw <- muffle_rank_warnings(
    ivreg2(lwage ~ exper + expersq + married + union | hours | educ + black,
           data = wp_psd, sw = TRUE, ivar = "nr", psd = "psd0")
  )
  # SW omega is per-observation scale; same uniform identity applies
  v_rebuilt <- vcov_from_stored_s(fit_sw)
  expect_equal(unname(vcov(fit_sw)), v_rebuilt, tolerance = 1e-8)
})

test_that("KP identification stats are never psd-corrected (ranktest parity)", {
  skip_if_not(file.exists(wp_psd_path))

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
  skip_if_not(file.exists(card_path))

  fit_null <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                     data = card, vcov = "iid")
  fit_psd <- ivreg2(lwage ~ exper + expersq | educ | nearc4,
                    data = card, vcov = "iid", psd = "psd0")
  # Identical code path: psd is inert for the iid VCV
  expect_identical(vcov(fit_null), vcov(fit_psd))
})
