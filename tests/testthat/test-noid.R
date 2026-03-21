# ============================================================================
# Tests: noid parameter — suppress identification diagnostics (Q4)
# ============================================================================

# --- Helper: standard IV fit on Card data ---
card_path <- fixture_path("card_data.csv")
skip_if_not(file.exists(card_path), "Card fixture not found")
card <- read.csv(card_path)

base_args <- list(
  formula = lwage ~ exper + expersq + black + south + smsa + smsa66 |
    educ | nearc4 + nearc2,
  data = card
)

# ============================================================================
# Basic suppression: noid = TRUE omits underid, weak_id, weak_id_sy, redundancy
# ============================================================================

test_that("noid = TRUE suppresses identification and redundancy diagnostics", {
  fit <- do.call(ivreg2, c(base_args, list(noid = TRUE)))
  diag <- fit$diagnostics
  expect_null(diag$underid)
  expect_null(diag$weak_id)
  expect_null(diag$weak_id_robust)
  expect_null(diag$weak_id_sy)
  expect_null(diag$redundancy)
})

# ============================================================================
# Non-suppressed diagnostics still present
# ============================================================================

test_that("noid = TRUE still computes overid, AR, SW, endogeneity, first_stage", {
  fit <- do.call(ivreg2, c(base_args, list(noid = TRUE)))
  diag <- fit$diagnostics
  expect_false(is.null(diag$overid))
  expect_false(is.null(diag$anderson_rubin))
  expect_false(is.null(diag$stock_wright))
  expect_false(is.null(diag$endogeneity))
  expect_false(is.null(fit$first_stage))
})

# ============================================================================
# Model F still computed
# ============================================================================

test_that("noid = TRUE still computes model F", {
  fit <- do.call(ivreg2, c(base_args, list(noid = TRUE, small = TRUE)))
  expect_false(is.null(fit$model_f))
  expect_true(is.finite(fit$model_f))
})

# ============================================================================
# Coefficients and VCV unchanged
# ============================================================================

test_that("noid does not change coefficients or VCV", {
  fit_noid <- do.call(ivreg2, c(base_args, list(noid = TRUE)))
  fit_base <- do.call(ivreg2, c(base_args, list(noid = FALSE)))
  expect_equal(coef(fit_noid), coef(fit_base))
  expect_equal(vcov(fit_noid), vcov(fit_base))
})

# ============================================================================
# glance() handles suppressed diagnostics (NA for suppressed, non-NA for others)
# ============================================================================

test_that("glance() returns NA for suppressed diagnostics with noid = TRUE", {
  fit <- do.call(ivreg2, c(base_args, list(noid = TRUE)))
  gl <- glance(fit)
  expect_true(is.na(gl$weak_id_stat))
  expect_true(is.na(gl$weak_id_robust_stat))
  expect_true(is.na(gl$underid_stat))
  expect_true(is.na(gl$underid_p))
  # Non-suppressed should be non-NA
  expect_false(is.na(gl$overid_stat))
  expect_false(is.na(gl$endogeneity_stat))
})

# ============================================================================
# summary() prints without error
# ============================================================================

test_that("summary() works with noid = TRUE", {
  fit <- do.call(ivreg2, c(base_args, list(noid = TRUE)))
  expect_output(print(summary(fit)))
})

# ============================================================================
# Input validation
# ============================================================================

test_that("noid rejects invalid inputs", {
  expect_error(
    do.call(ivreg2, c(base_args, list(noid = "yes"))),
    "`noid` must be TRUE or FALSE"
  )
  expect_error(
    do.call(ivreg2, c(base_args, list(noid = NA))),
    "`noid` must be TRUE or FALSE"
  )
  expect_error(
    do.call(ivreg2, c(base_args, list(noid = c(TRUE, FALSE)))),
    "`noid` must be TRUE or FALSE"
  )
})

# ============================================================================
# noid + b0 interaction: both suppress identification diagnostics, no conflict
# ============================================================================

test_that("noid = TRUE + b0 works without conflict", {
  fit_base <- do.call(ivreg2, c(base_args, list(noid = FALSE)))
  b0_vec <- coef(fit_base)
  expect_no_error(
    do.call(ivreg2, c(base_args, list(noid = TRUE, b0 = b0_vec)))
  )
})

# ============================================================================
# noid + redundant: redundancy test skipped silently
# ============================================================================

test_that("noid = TRUE + redundant skips redundancy test silently", {
  fit <- do.call(ivreg2, c(base_args, list(noid = TRUE, redundant = "nearc4")))
  expect_null(fit$diagnostics$redundancy)
})

test_that("noid = TRUE + invalid redundant still errors (validation unconditional)", {
  expect_error(
    do.call(ivreg2, c(base_args, list(noid = TRUE, redundant = "not_a_var"))),
    "not in the excluded instrument list"
  )
})

# ============================================================================
# b0 + invalid argument names still error (validation outside b0 gate)
# ============================================================================

test_that("b0 + invalid redundant still errors", {
  fit_base <- do.call(ivreg2, c(base_args, list(noid = FALSE)))
  b0_vec <- coef(fit_base)
  expect_error(
    do.call(ivreg2, c(base_args, list(b0 = b0_vec, redundant = "not_a_var"))),
    "not in the excluded instrument list"
  )
})

test_that("b0 + invalid endog still errors", {
  fit_base <- do.call(ivreg2, c(base_args, list(noid = FALSE)))
  b0_vec <- coef(fit_base)
  expect_error(
    do.call(ivreg2, c(base_args, list(b0 = b0_vec, endog = "not_a_var"))),
    "not in the endogenous list"
  )
})

test_that("b0 + invalid orthog still errors", {
  fit_base <- do.call(ivreg2, c(base_args, list(noid = FALSE)))
  b0_vec <- coef(fit_base)
  expect_error(
    do.call(ivreg2, c(base_args, list(b0 = b0_vec, orthog = "not_a_var"))),
    "not in the instrument list"
  )
})

test_that("b0 + orthog + partial does not error (orthog suppressed by b0)", {
  # Regression test: partialled_out check must not fire when b0 suppresses orthog
  fit_base <- ivreg2(lwage ~ exper + expersq + black + south + smsa + smsa66 |
                       educ | nearc4 + nearc2, data = card, partial = "smsa")
  b0_vec <- coef(fit_base)
  expect_no_error(
    ivreg2(lwage ~ exper + expersq + black + south + smsa + smsa66 |
             educ | nearc4 + nearc2, data = card,
           partial = "smsa", orthog = "smsa", b0 = b0_vec)
  )
})

test_that("b0 + partial does not warn about CUE FWL invariance", {
  # b0 evaluates at fixed values — no optimization, so FWL warning doesn't apply
  fit_base <- ivreg2(lwage ~ exper + expersq + black + south + smsa + smsa66 |
                       educ | nearc4 + nearc2, data = card, partial = "smsa")
  b0_vec <- coef(fit_base)
  expect_no_warning(
    ivreg2(lwage ~ exper + expersq + black + south + smsa + smsa66 |
             educ | nearc4 + nearc2, data = card,
           partial = "smsa", b0 = b0_vec)
  )
})

# ============================================================================
# noid stored on object
# ============================================================================

test_that("noid value is stored on the fitted object", {
  fit_t <- do.call(ivreg2, c(base_args, list(noid = TRUE)))
  fit_f <- do.call(ivreg2, c(base_args, list(noid = FALSE)))
  expect_true(fit_t$noid)
  expect_false(fit_f$noid)
})
