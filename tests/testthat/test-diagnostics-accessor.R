# ============================================================================
# Tests: diagnostics() accessor for ivreg2 objects
# ============================================================================
# No Stata fixtures needed here: every expected value is read directly from
# the fitted object's own `$diagnostics` slots, which are already fixture-
# validated elsewhere. This file only checks that diagnostics() flattens
# those slots correctly.

data(card)


# ============================================================================
# Overidentified robust fit with endog/orthog/redundant
# ============================================================================

test_that("diagnostics() lists the documented keys in order for a full fit", {
  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
    data = card, vcov = "robust",
    endog = "educ", orthog = "nearc2", redundant = "nearc2"
  )
  d <- diagnostics(fit)

  expect_s3_class(d, "tbl_df")
  expect_named(d, c("test", "test_name", "statistic", "df", "df2",
                     "p_value", "tested_vars", "note"))
  expect_type(d$test, "character")
  expect_type(d$test_name, "character")
  expect_type(d$statistic, "double")
  expect_type(d$df, "integer")
  expect_type(d$df2, "integer")
  expect_type(d$p_value, "double")
  expect_type(d$tested_vars, "character")
  expect_type(d$note, "character")
  # This fit trips none of the Stata-quirk triggers, so every note is NA.
  expect_true(all(is.na(d$note)))

  expected_keys <- c(
    "underid", "weak_id", "weak_id_robust",
    "sy_iv_size_10", "sy_iv_size_15", "sy_iv_size_20", "sy_iv_size_25",
    "overid", "anderson_rubin_f", "anderson_rubin_chi2", "stock_wright",
    "endogeneity", "orthog", "redundancy"
  )
  expect_identical(d$test, expected_keys)

  diag <- fit$diagnostics

  underid_row <- d[d$test == "underid", ]
  expect_equal(underid_row$statistic, diag$underid$stat)
  expect_equal(underid_row$df, diag$underid$df)
  expect_equal(underid_row$p_value, diag$underid$p)
  expect_equal(underid_row$test_name, diag$underid$test_name)

  weak_id_row <- d[d$test == "weak_id", ]
  expect_equal(weak_id_row$statistic, diag$weak_id$stat)
  expect_true(is.na(weak_id_row$df))
  expect_true(is.na(weak_id_row$p_value))

  weak_id_robust_row <- d[d$test == "weak_id_robust", ]
  expect_equal(weak_id_robust_row$statistic, diag$weak_id_robust$stat)

  overid_row <- d[d$test == "overid", ]
  expect_equal(overid_row$statistic, diag$overid$stat)
  expect_equal(overid_row$df, diag$overid$df)
  expect_equal(overid_row$p_value, diag$overid$p)
  expect_equal(overid_row$test_name, diag$overid$test_name)

  ar_f_row <- d[d$test == "anderson_rubin_f", ]
  expect_equal(ar_f_row$statistic, diag$anderson_rubin$f_stat)
  expect_equal(ar_f_row$df, diag$anderson_rubin$f_df1)
  expect_equal(ar_f_row$df2, diag$anderson_rubin$f_df2)
  expect_equal(ar_f_row$p_value, diag$anderson_rubin$f_p)

  ar_chi2_row <- d[d$test == "anderson_rubin_chi2", ]
  expect_equal(ar_chi2_row$statistic, diag$anderson_rubin$chi2_stat)
  expect_equal(ar_chi2_row$df, diag$anderson_rubin$chi2_df)
  expect_true(is.na(ar_chi2_row$df2))
  expect_equal(ar_chi2_row$p_value, diag$anderson_rubin$chi2_p)

  sw_row <- d[d$test == "stock_wright", ]
  expect_equal(sw_row$statistic, diag$stock_wright$stat)
  expect_equal(sw_row$df, diag$stock_wright$df)
  expect_equal(sw_row$p_value, diag$stock_wright$p)

  endog_row <- d[d$test == "endogeneity", ]
  expect_equal(endog_row$statistic, diag$endogeneity$stat)
  expect_equal(endog_row$df, diag$endogeneity$df)
  expect_equal(endog_row$p_value, diag$endogeneity$p)
  expect_equal(endog_row$test_name, diag$endogeneity$test_name)
  expect_equal(endog_row$tested_vars,
               paste(diag$endogeneity$tested_vars, collapse = ", "))

  orthog_row <- d[d$test == "orthog", ]
  expect_equal(orthog_row$statistic, diag$orthog$stat)
  expect_equal(orthog_row$df, diag$orthog$df)
  expect_equal(orthog_row$p_value, diag$orthog$p)
  expect_equal(orthog_row$tested_vars,
               paste(diag$orthog$tested_vars, collapse = ", "))

  redundancy_row <- d[d$test == "redundancy", ]
  expect_equal(redundancy_row$statistic, diag$redundancy$stat)
  expect_equal(redundancy_row$df, diag$redundancy$df)
  expect_equal(redundancy_row$p_value, diag$redundancy$p)
  expect_equal(redundancy_row$tested_vars,
               paste(diag$redundancy$tested_vars, collapse = ", "))

  # tested_vars is NA everywhere except the endog/orthog/redundancy rows.
  non_tested_rows <- d[!d$test %in% c("endogeneity", "orthog", "redundancy"), ]
  expect_true(all(is.na(non_tested_rows$tested_vars)))
})


# ============================================================================
# Consistency with glance()
# ============================================================================

test_that("diagnostics() agrees with glance() on overid/underid/weak_id", {
  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
    data = card, vcov = "robust"
  )
  d <- diagnostics(fit)
  gl <- glance(fit)

  expect_equal(d[d$test == "overid", ]$statistic, gl$overid_stat)
  expect_equal(d[d$test == "overid", ]$p_value, gl$overid_p)
  expect_equal(d[d$test == "underid", ]$statistic, gl$underid_stat)
  expect_equal(d[d$test == "underid", ]$p_value, gl$underid_p)
  expect_equal(d[d$test == "weak_id", ]$statistic, gl$weak_id_stat)
  expect_equal(d[d$test == "weak_id_robust", ]$statistic,
               gl$weak_id_robust_stat)
})


# ============================================================================
# Stock-Yogo rows
# ============================================================================

test_that("diagnostics() SY rows match the weak_id_sy data frame", {
  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
    data = card, vcov = "robust"
  )
  sy <- fit$diagnostics$weak_id_sy
  d <- diagnostics(fit)
  sy_rows <- d[grepl("^sy_", d$test), ]

  expect_equal(nrow(sy_rows), nrow(sy))

  # Literal keys, deliberately NOT derived with the implementation's own
  # slug regex, so a slug bug cannot mirror itself into the expectation.
  expect_identical(
    sy_rows$test,
    c("sy_iv_size_10", "sy_iv_size_15", "sy_iv_size_20", "sy_iv_size_25")
  )
  expect_equal(sy_rows$statistic, sy$critical_value)
  expect_true(all(is.na(sy_rows$p_value)))
  expect_true(all(is.na(sy_rows$df)))
})


# ============================================================================
# LIML iid fit
# ============================================================================

test_that("diagnostics() reports AR overid rows and omits weak_id_robust for LIML iid", {
  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
    data = card, method = "liml"
  )
  d <- diagnostics(fit)
  diag <- fit$diagnostics

  expect_true(is.null(diag$weak_id_robust))
  expect_false("weak_id_robust" %in% d$test)

  # A LIML fit gets the "LIML size" Stock-Yogo table, so this pins the slug
  # for a non-"IV size" type with a literal key.
  expect_true("sy_liml_size_10" %in% d$test)
  expect_equal(
    d[d$test == "sy_liml_size_10", ]$statistic,
    diag$weak_id_sy$critical_value[diag$weak_id_sy$threshold == "10%"]
  )

  expect_true(all(c("anderson_rubin_overid_lr", "anderson_rubin_overid_lin")
                  %in% d$test))

  lr_row <- d[d$test == "anderson_rubin_overid_lr", ]
  expect_equal(lr_row$statistic, diag$anderson_rubin_overid$lr_stat)
  expect_equal(lr_row$df, diag$anderson_rubin_overid$df)
  expect_equal(lr_row$p_value, diag$anderson_rubin_overid$lr_p)
  expect_identical(lr_row$test_name, "Anderson-Rubin overidentification (LR)")

  lin_row <- d[d$test == "anderson_rubin_overid_lin", ]
  expect_equal(lin_row$statistic, diag$anderson_rubin_overid$lin_stat)
  expect_equal(lin_row$df, diag$anderson_rubin_overid$df)
  expect_equal(lin_row$p_value, diag$anderson_rubin_overid$lin_p)
  expect_identical(lin_row$test_name,
                    "Anderson-Rubin overidentification (linearized)")
})


# ============================================================================
# OLS fit
# ============================================================================

test_that("diagnostics() returns a zero-row tibble for OLS", {
  fit <- ivreg2(lwage ~ exper + expersq, data = card)
  d <- diagnostics(fit)

  expect_null(fit$diagnostics)
  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 0L)
  expect_named(d, c("test", "test_name", "statistic", "df", "df2",
                     "p_value", "tested_vars", "note"))
  expect_type(d$test, "character")
  expect_type(d$test_name, "character")
  expect_type(d$statistic, "double")
  expect_type(d$df, "integer")
  expect_type(d$df2, "integer")
  expect_type(d$p_value, "double")
  expect_type(d$tested_vars, "character")
  expect_type(d$note, "character")
})


# ============================================================================
# b0 fit
# ============================================================================

test_that("diagnostics() reports only the overid row for a b0 fit", {
  fit_2sls <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
    data = card, vcov = "robust"
  )
  fit_b0 <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
    data = card, vcov = "robust", b0 = coef(fit_2sls)
  )
  d <- diagnostics(fit_b0)

  expect_identical(names(fit_b0$diagnostics), "overid")
  expect_identical(d$test, "overid")
  expect_equal(d$statistic, fit_b0$diagnostics$overid$stat)
  expect_equal(d$df, fit_b0$diagnostics$overid$df)
  expect_equal(d$p_value, fit_b0$diagnostics$overid$p)
})


# ============================================================================
# Just-identified fit
# ============================================================================

test_that("diagnostics() reports overid with df 0 and NA p-value when just-identified", {
  fit <- ivreg2(
    lwage ~ exper + expersq + black + south | educ | nearc4,
    data = card, vcov = "robust"
  )
  d <- diagnostics(fit)
  overid_row <- d[d$test == "overid", ]

  expect_equal(nrow(overid_row), 1L)
  expect_identical(overid_row$df, 0L)
  expect_true(is.na(overid_row$p_value))
  expect_equal(overid_row$statistic, fit$diagnostics$overid$stat)
})
