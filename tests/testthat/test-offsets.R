# ==========================================================================
# Tests for offset() rejection
# ==========================================================================
# Stata's ivreg2 has no offset concept, so there is no parity target to
# implement against (Frank, 2026-07-10). ivreg2() must reject offset() terms
# with a clear error rather than silently ignoring them.

data(card)

offset_msg <- paste0(
  "offset\\(\\) terms are not supported by ivreg2\\(\\)\\. Include the ",
  "offset variable as a regressor, or subtract it from the response ",
  "before fitting\\."
)

test_that("offset() in an IV formula's exogenous part is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + offset(black) | educ | nearc4,
           data = card),
    offset_msg
  )
})

test_that("offset() in an OLS-style formula (no | parts) is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq + offset(black), data = card),
    offset_msg
  )
})

test_that("offset() in the excluded-instruments part is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq | educ | nearc4 + offset(black),
           data = card),
    offset_msg
  )
})

test_that("offset() in the endogenous part is rejected", {
  expect_error(
    ivreg2(lwage ~ exper + expersq | offset(black) + educ | nearc4,
           data = card),
    offset_msg
  )
})

test_that("a regressor merely named offset does not trigger the error", {
  d <- card
  d$offset <- d$black
  result <- ivreg2(lwage ~ exper + expersq + offset | educ | nearc4,
                    data = d)
  expect_s3_class(result, "ivreg2")
  expect_true("offset" %in% names(coef(result)))
})
