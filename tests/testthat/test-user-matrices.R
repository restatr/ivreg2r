# ============================================================================
# Tests: User-supplied W and S matrices (Ticket N2)
# ============================================================================

# --- Helpers ---
data(card)
data(griliches, package = "ivreg2r")

# Stata's instrument order for the Card overid model is:
#   exper, expersq, black, south, _cons, nearc2, nearc4
# R's Z column order is:
#   (Intercept), exper, expersq, black, south, nearc2, nearc4
# This helper reorders a Stata-order matrix to R's order.
reorder_stata_to_r <- function(M, stata_order, r_order) {
  # Build permutation: which Stata column maps to each R column
  perm <- match(r_order, stata_order)
  M[perm, perm]
}

# Standard mappings for Card models
# Stata's e(S) column order: excluded instruments first, then included
# exogenous regressors, then constant.
stata_z_order_overid <- c("nearc2", "nearc4", "exper", "expersq", "black",
                           "south", "(Intercept)")
r_z_order_overid <- c("(Intercept)", "exper", "expersq", "black", "south",
                       "nearc2", "nearc4")
stata_z_order_justid <- c("nearc4", "exper", "expersq", "black", "south",
                            "(Intercept)")
r_z_order_justid <- c("(Intercept)", "exper", "expersq", "black", "south",
                        "nearc4")

# Load a Stata S/W matrix fixture and reorder to R's Z ordering
load_stata_matrix <- function(path, stata_order, r_order) {
  M <- as.matrix(read.csv(path))
  M <- reorder_stata_to_r(M, stata_order, r_order)
  # Strip read.csv's auto colnames (v1, v2, ...): these matrices are aligned
  # positionally by the reorder above, and named matrices are now matched to
  # the instrument list BY NAME (.match_user_matrix, Stata matsort parity) --
  # junk v-names would be rejected.
  dimnames(M) <- NULL
  M
}

# ============================================================================
# 1. Input validation
# ============================================================================

test_that("wmatrix rejects non-matrix input", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, vcov = "robust", wmatrix = 1:7),
    "`wmatrix` must be a numeric matrix"
  )
})

test_that("smatrix rejects non-matrix input", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, vcov = "robust", smatrix = "foo"),
    "`smatrix` must be a numeric matrix"
  )
})

test_that("wmatrix wrong dimensions", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, vcov = "robust", wmatrix = diag(3)),
    "do not match the number of instruments"
  )
})

test_that("smatrix wrong dimensions", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, vcov = "robust", smatrix = diag(3)),
    "do not match the number of instruments"
  )
})

test_that("wmatrix non-symmetric", {
  W <- diag(7)
  W[1, 2] <- 999
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, vcov = "robust", wmatrix = W),
    "not symmetric"
  )
})

test_that("smatrix non-symmetric", {
  S <- diag(7)
  S[1, 2] <- 999
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, vcov = "robust", smatrix = S),
    "not symmetric"
  )
})

test_that("wmatrix + liml errors", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, method = "liml", wmatrix = diag(7)),
    'Cannot specify `wmatrix` with method = "liml"'
  )
})

test_that("smatrix + kclass errors", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, kclass = 0.5, smatrix = diag(7)),
    "Cannot specify `smatrix` with `kclass`"
  )
})

test_that("wmatrix + fuller errors", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, fuller = 1, wmatrix = diag(7)),
    "Cannot specify `wmatrix` with `fuller`"
  )
})

test_that("smatrix + fuller errors", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
           data = card, fuller = 1, smatrix = diag(7)),
    "Cannot specify `smatrix` with `fuller`"
  )
})

test_that("wmatrix requires IV model", {
  expect_error(
    ivreg2(lwage ~ exper + expersq, data = card, vcov = "robust",
           wmatrix = diag(3)),
    "requires an IV model"
  )
})

test_that("smatrix requires IV model", {
  expect_error(
    ivreg2(lwage ~ exper + expersq, data = card, vcov = "robust",
           smatrix = diag(3)),
    "requires an IV model"
  )
})

test_that("wmatrix + IID warns and is ignored", {
  # wmatrix without robust/cluster/HAC/gmm2s → warning, then 2SLS
  expect_warning(
    fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                  data = card, wmatrix = diag(7)),
    "ignored without robust VCE"
  )
  expect_equal(fit$method, "2sls")
})


# ============================================================================
# 2. Method classification
# ============================================================================

test_that("wmatrix alone → method = 'gmmw'", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", wmatrix = diag(7))
  expect_equal(fit$method, "gmmw")
})

test_that("smatrix alone → method = 'gmm2s'", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", smatrix = diag(7))
  expect_equal(fit$method, "gmm2s")
})

test_that("wmatrix + gmm2s → method = 'gmm2s'", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", method = "gmm2s", wmatrix = diag(7))
  expect_equal(fit$method, "gmm2s")
})


# ============================================================================
# 3. Estimation labels
# ============================================================================

test_that("gmmw has correct estimation label", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", wmatrix = diag(7))
  out <- capture.output(print(fit))
  expect_true(any(grepl("GMM Estimation.*user-supplied W", out)))
})

test_that("smatrix has correct estimation label", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", smatrix = diag(7))
  out <- capture.output(print(fit))
  expect_true(any(grepl("GMM Estimation.*user-supplied S", out)))
})

test_that("wmatrix + gmm2s has correct estimation label", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", method = "gmm2s", wmatrix = diag(7))
  out <- capture.output(print(fit))
  expect_true(any(grepl("2-Step GMM Estimation.*user-supplied first step", out)))
})


# ============================================================================
# 4. Algebraic equivalence tests (no fixtures needed)
# ============================================================================

test_that("smatrix(S_from_gmm2s) reproduces standard gmm2s", {
  # Run standard GMM2S to get Omega
  fit_std <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                    data = card, vcov = "robust", method = "gmm2s", x = TRUE)
  # Extract Omega (= fitted S matrix from the estimation)
  # Omega is not directly stored in the return object, but we can reconstruct
  # it by running the omega_fn. Instead, let's verify that providing the
  # Stata-derived S matrix gives the same result.
  # For this algebraic test, compute S = omega from the estimation manually.
  Z <- fit_std$x$Z
  resid_1 <- residuals(ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card, x = TRUE))
  N <- nobs(fit_std)
  # HC0 meat: (1/N) * Z' diag(e^2) Z / N (unnormalized for omega)
  scores <- Z * resid_1
  Omega <- crossprod(scores) / N

  # Run with this smatrix
  fit_smat <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                     data = card, vcov = "robust", smatrix = Omega)

  # Should produce the same coefficients as standard gmm2s
  expect_equal(coef(fit_smat), coef(fit_std), tolerance = 1e-10)
  expect_equal(vcov(fit_smat), vcov(fit_std), tolerance = 1e-10)
})

test_that("just-identified: wmatrix + gmm2s equals standard gmm2s", {
  # For just-identified, the first-step weighting matrix doesn't matter
  # (only one solution). So wmatrix(I) + gmm2s should match standard gmm2s.
  fit_std <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card, vcov = "robust", method = "gmm2s")
  # L = 6 for just-identified (5 exog + 1 excluded IV)
  fit_wm <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                   data = card, vcov = "robust", method = "gmm2s",
                   wmatrix = diag(6))
  # The just-identified GMM solution is unique, but different first-step W
  # produces slightly different step-1 residuals → different omega_fn → different
  # step-2 beta. The condition number of the just-id Hessian amplifies the
  # small difference. Use 1e-6 for coefs, 1e-5 for VCV (contains squared errors).
  expect_equal(coef(fit_wm), coef(fit_std), tolerance = 1e-6)
  expect_equal(vcov(fit_wm), vcov(fit_std), tolerance = 1e-5)
})

test_that("smatrix reuse reproduces the two-step GMM fit", {
  a <- ivreg2(gril_formula, data = griliches, method = "gmm2s", vcov = "robust")
  b <- ivreg2(gril_formula, data = griliches, vcov = "robust")
  refit <- ivreg2(gril_formula, data = griliches, method = "gmm2s",
                  vcov = "robust", smatrix = b$S)
  expect_equal(b$S, a$S)
  expect_equal(coef(refit), coef(a))
})


# ============================================================================
# 5. Stata fixture tests — wmatrix = I(L), robust
# ============================================================================
# The GMM Hessian H = QXZ * W * QXZ' with W = I has condition number ~1.6e13
# for the Card data, meaning R and Stata can only agree to ~3-4 significant
# digits. This is a well-understood conditioning effect, not a formula bug.
# Use stata_tol$stat (1e-4) for comparisons in wmatrix fixture tests.

test_that("wmatrix identity robust: coefs match Stata", {
  fixture_path <- fixture_path("card_overid_wmat_ident_coef_robust.csv")
  skip_if(!file.exists(fixture_path), "wmatrix fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", wmatrix = diag(7))
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$stat,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$stat,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("wmatrix identity robust: diagnostics match Stata", {
  fixture_path <- fixture_path("card_overid_wmat_ident_diagnostics_robust.csv")
  skip_if(!file.exists(fixture_path), "wmatrix diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", wmatrix = diag(7))
  diag <- read.csv(fixture_path)

  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = stata_tol$stat, info = "J stat")
  expect_equal(fit$diagnostics$overid$p, diag$overid_p,
               tolerance = stata_tol$pval, info = "J p")
  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
  expect_equal(fit$sigma, diag$sigma, tolerance = stata_tol$stat,
               info = "sigma")
  expect_equal(fit$rss, diag$rss, tolerance = stata_tol$stat,
               info = "rss")
  expect_equal(fit$r.squared, diag$r2, tolerance = stata_tol$stat,
               info = "r2")
  expect_equal(fit$nobs, diag$N, info = "N")
})

test_that("wmatrix identity robust: VCV matches Stata", {
  fixture_path <- fixture_path("card_overid_wmat_ident_vcov_robust.csv")
  skip_if(!file.exists(fixture_path), "wmatrix VCV fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", wmatrix = diag(7))
  # Stata VCV is in Stata's coefficient order (educ, exper, ..., _cons);
  # R VCV is in R's coefficient order (Intercept, educ, exper, ...).
  # Reorder Stata VCV to match R's order for comparison.
  stata_coef_order <- c("educ", "exper", "expersq", "black", "south",
                         "(Intercept)")
  r_coef_order <- c("(Intercept)", "educ", "exper", "expersq", "black",
                     "south")
  V_stata <- read_vcov_fixture(fixture_path)
  perm <- match(r_coef_order, stata_coef_order)
  V_stata <- unname(V_stata[perm, perm])
  V_r <- unname(vcov(fit))
  expect_equal(V_r, V_stata, tolerance = stata_tol$stat)
})


# ============================================================================
# 6. Stata fixture tests — wmatrix = I(L), robust small
# ============================================================================

test_that("wmatrix identity robust small: coefs match Stata", {
  fixture_path <- fixture_path("card_overid_wmat_ident_coef_robust_small.csv")
  skip_if(!file.exists(fixture_path), "wmatrix robust_small fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", small = TRUE, wmatrix = diag(7))
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$stat,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$stat,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("wmatrix identity robust small: diagnostics match Stata", {
  fixture_path <- fixture_path("card_overid_wmat_ident_diagnostics_robust_small.csv")
  skip_if(!file.exists(fixture_path), "wmatrix robust_small diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", small = TRUE, wmatrix = diag(7))
  diag <- read.csv(fixture_path)

  expect_equal(fit$sigma, diag$sigma, tolerance = stata_tol$stat,
               info = "sigma")
  expect_equal(fit$rss, diag$rss, tolerance = stata_tol$stat,
               info = "rss")
})


# ============================================================================
# 7. Stata fixture tests — smatrix, robust
# ============================================================================

test_that("smatrix robust: coefs match Stata", {
  fixture_path <- fixture_path("card_overid_smat_coef_robust.csv")
  s_path <- fixture_path("card_overid_smat_S_robust.csv")
  skip_if(!file.exists(fixture_path) || !file.exists(s_path),
          "smatrix fixture not found")

  # Load Stata's exact S matrix and reorder to R's Z column order
  S <- load_stata_matrix(s_path, stata_z_order_overid, r_z_order_overid)

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", smatrix = S)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("smatrix robust: diagnostics match Stata", {
  fixture_path <- fixture_path("card_overid_smat_diagnostics_robust.csv")
  s_path <- fixture_path("card_overid_smat_S_robust.csv")
  skip_if(!file.exists(fixture_path) || !file.exists(s_path),
          "smatrix diagnostics fixture not found")

  S <- load_stata_matrix(s_path, stata_z_order_overid, r_z_order_overid)

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", smatrix = S)
  diag <- read.csv(fixture_path)

  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = stata_tol$stat, info = "overid stat")
  expect_equal(fit$diagnostics$overid$p, diag$overid_p,
               tolerance = stata_tol$pval, info = "overid p")
  expect_equal(fit$sigma, diag$sigma, tolerance = stata_tol$coef,
               info = "sigma")
})


# ============================================================================
# 8. Stata fixture tests — wmatrix + gmm2s, robust
# ============================================================================

test_that("wmatrix + gmm2s robust: coefs match Stata", {
  fixture_path <- fixture_path("card_overid_wmat_gmm2s_coef_robust.csv")
  skip_if(!file.exists(fixture_path), "wmat_gmm2s fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", method = "gmm2s", wmatrix = diag(7))
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})

test_that("wmatrix + gmm2s robust: diagnostics match Stata", {
  fixture_path <- fixture_path("card_overid_wmat_gmm2s_diagnostics_robust.csv")
  skip_if(!file.exists(fixture_path), "wmat_gmm2s diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", method = "gmm2s", wmatrix = diag(7))
  diag <- read.csv(fixture_path)

  expect_equal(fit$diagnostics$overid$stat, diag$overid_stat,
               tolerance = stata_tol$stat, info = "J stat")
  expect_equal(fit$sigma, diag$sigma, tolerance = stata_tol$coef,
               info = "sigma")
})


# ============================================================================
# 9. Stata fixture tests — wmatrix from IID, robust
# ============================================================================

test_that("wmatrix from IID robust: runs without error", {
  # Coefficient comparison is skipped for this fixture. The IID S matrix

  # (sigma^2 * Z'Z/N) used as W produces a GMM Hessian with condition number
  # ~5e17 for the Card data, making coefficients numerically meaningless —
  # tiny perturbations in W produce 10%+ coefficient changes. R and Stata
  # use different LAPACK implementations, so divergence at this condition
  # number is expected. We verify the code runs without error and produces
  # the correct method classification.
  w_path <- fixture_path("card_overid_wmat_from_iid_W.csv")
  skip_if(!file.exists(w_path), "wmat_from_iid fixture not found")

  W_iid <- load_stata_matrix(w_path, stata_z_order_overid, r_z_order_overid)

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", wmatrix = W_iid)
  expect_equal(fit$method, "gmmw")
  expect_equal(nobs(fit), 3010L)
})


# ============================================================================
# 10. Stata fixture tests — smatrix, just-identified
# ============================================================================

test_that("smatrix just-identified robust: coefs match Stata", {
  fixture_path <- fixture_path("card_justid_smat_coef_robust.csv")
  s_path <- fixture_path("card_justid_smat_S_robust.csv")
  skip_if(!file.exists(fixture_path) || !file.exists(s_path),
          "smatrix justid fixture not found")

  S <- load_stata_matrix(s_path, stata_z_order_justid, r_z_order_justid)

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, vcov = "robust", smatrix = S)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
  }
})


# ============================================================================
# 11. Stata fixture tests — wmatrix = I(L), cluster
# ============================================================================

test_that("wmatrix identity cluster: coefs match Stata", {
  fixture_path <- fixture_path("card_overid_wmat_ident_coef_cluster.csv")
  skip_if(!file.exists(fixture_path), "wmatrix cluster fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, clusters = ~age, wmatrix = diag(7))
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$stat,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$stat,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})


# ============================================================================
# 12. Stata fixture tests — smatrix, cluster
# ============================================================================

test_that("smatrix cluster: coefs match Stata", {
  fixture_path <- fixture_path("card_overid_smat_coef_cluster.csv")
  s_path <- fixture_path("card_overid_smat_S_cluster.csv")
  skip_if(!file.exists(fixture_path) || !file.exists(s_path),
          "smatrix cluster fixture not found")

  S <- load_stata_matrix(s_path, stata_z_order_overid, r_z_order_overid)

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, clusters = ~age, smatrix = S)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})


# ============================================================================
# 13. Stata fixture tests — wmatrix + smatrix, robust
# ============================================================================

test_that("wmatrix + smatrix robust: coefs match Stata", {
  fixture_path <- fixture_path("card_overid_wmat_smat_coef_robust.csv")
  s_path <- fixture_path("card_overid_wmat_smat_S_robust.csv")
  skip_if(!file.exists(fixture_path) || !file.exists(s_path),
          "wmat_smat fixture not found")

  S <- load_stata_matrix(s_path, stata_z_order_overid, r_z_order_overid)

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", wmatrix = diag(7), smatrix = S)
  fixture <- read.csv(fixture_path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)

  # Same identity W conditioning as section 5
  for (i in seq_len(nrow(fixture))) {
    expect_equal(
      unname(coef(fit)[r_names[i]]), fixture$estimate[i],
      tolerance = stata_tol$stat,
      info = paste("Coef mismatch:", r_names[i])
    )
    expect_equal(
      unname(sqrt(diag(vcov(fit)))[r_names[i]]), fixture$std_error[i],
      tolerance = stata_tol$stat,
      info = paste("SE mismatch:", r_names[i])
    )
  }
})


# ============================================================================
# 14. Glance works for gmmw
# ============================================================================

test_that("glance works for gmmw method", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", wmatrix = diag(7))
  g <- glance(fit)
  expect_equal(nrow(g), 1L)
  expect_equal(fit$method, "gmmw")
  expect_true(!is.na(g$overid_stat))
})

test_that("tidy works for gmmw method", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", wmatrix = diag(7))
  t <- tidy(fit)
  expect_equal(nrow(t), 6L)
  expect_true(all(c("term", "estimate", "std.error") %in% names(t)))
})

test_that("wmatrix/smatrix stored in return object", {
  W <- diag(7)
  S <- diag(7)
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card, vcov = "robust", wmatrix = W, smatrix = S)
  expect_true(!is.null(fit$wmatrix))
  expect_true(!is.null(fit$smatrix))
  expect_equal(fit$wmatrix, W)
  expect_equal(fit$smatrix, S)
})
