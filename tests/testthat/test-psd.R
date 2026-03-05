# ============================================================================
# Tests: PSD Corrections — Ticket P2
# ============================================================================

# --- Helpers ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)
card_path <- file.path(fixture_dir, "card_data.csv")

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
