# ==========================================================================
# Tests for .parse_formula() and helpers
# ==========================================================================

# -- Helper: build a small data frame for most tests -----------------------
make_test_data <- function(n = 100) {
  set.seed(42)
  data.frame(
    y  = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n),
    endo1 = rnorm(n),
    endo2 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n),
    z3 = rnorm(n),
    w  = runif(n, 0.5, 1.5),
    grp = sample(letters[1:5], n, replace = TRUE)
  )
}

# ==========================================================================
# 1. Formula structure validation
# ==========================================================================
test_that("1-part formula (OLS) is accepted", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 + x2, data = d)
  expect_s3_class(result, "parsed_formula")
  expect_false(result$is_iv)
})

test_that("3-part formula (IV) is accepted", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 | endo1 | z1, data = d)
  expect_s3_class(result, "parsed_formula")
  expect_true(result$is_iv)
})

test_that("2-part formula gives exact error message", {
  d <- make_test_data()
  expect_error(
    .parse_formula(y ~ x1 | z1, data = d),
    "Two-part formula detected. ivreg2\\(\\) uses a three-part formula: y ~ exog \\| endo \\| instruments. Did you mean ivreg::ivreg\\(\\)\\?"
  )
})

test_that("4-part formula gives exact error message with count", {
  d <- make_test_data()
  expect_error(
    .parse_formula(y ~ x1 | x2 | x3 | z1, data = d),
    "Formula has 4 parts; ivreg2\\(\\) expects 1 part \\(OLS\\) or 3 parts: y ~ exog \\| endo \\| instruments\\."
  )
})

test_that("5-part formula error includes correct count", {
  d <- make_test_data()
  expect_error(
    .parse_formula(y ~ x1 | x2 | x3 | z1 | z2, data = d),
    "Formula has 5 parts"
  )
})

test_that("3-part formula with empty part 2 gives exact error", {
  d <- make_test_data()
  # y ~ x1 | 0 | z1  — part 2 has no intercept and no variables
  expect_error(
    .parse_formula(y ~ x1 | 0 | z1, data = d),
    "Part 2 \\(endogenous regressors\\) is empty\\. Use a one-part formula for OLS: y ~ x1 \\+ x2\\."
  )
})

test_that("3-part formula with only-intercept part 2 gives empty error", {
  d <- make_test_data()
  # y ~ x1 | 1 | z1  — part 2 has only an intercept (which we strip)
  expect_error(
    .parse_formula(y ~ x1 | 1 | z1, data = d),
    "Part 2 \\(endogenous regressors\\) is empty"
  )
})


# ==========================================================================
# 2. OLS path
# ==========================================================================
test_that("OLS: Z is NULL, is_iv FALSE, K1 = 0", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 + x2, data = d)
  expect_null(result$Z)
  expect_false(result$is_iv)
  expect_equal(result$K1, 0L)
  expect_equal(result$K2, result$K)
})

test_that("OLS: correct X and y extraction", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 + x2, data = d)
  expect_equal(as.numeric(result$y), d$y)
  expect_equal(result$y_name, "y")
  expect_equal(ncol(result$X), 3L)  # intercept + x1 + x2
  expect_true("(Intercept)" %in% colnames(result$X))
  expect_true("x1" %in% colnames(result$X))
  expect_true("x2" %in% colnames(result$X))
})

test_that("OLS: dimension counts correct", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 + x2, data = d)
  expect_equal(result$N, 100L)
  expect_equal(result$K, 3L)  # intercept + 2 regressors
  expect_equal(result$L, 3L)  # same as K for OLS
})


# ==========================================================================
# 3. IV path
# ==========================================================================
test_that("IV: correct X/Z composition", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 | endo1 | z1, data = d)
  expect_equal(colnames(result$X), c("(Intercept)", "x1", "endo1"))
  expect_equal(colnames(result$Z), c("(Intercept)", "x1", "z1"))
})

test_that("IV: exogenous regressors appear in both X and Z", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 + x2 | endo1 | z1 + z2, data = d)
  exog_in_X <- all(c("(Intercept)", "x1", "x2") %in% colnames(result$X))
  exog_in_Z <- all(c("(Intercept)", "x1", "x2") %in% colnames(result$Z))
  expect_true(exog_in_X)
  expect_true(exog_in_Z)
})

test_that("IV: intercept propagated to Z", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 | endo1 | z1, data = d)
  expect_true("(Intercept)" %in% colnames(result$Z))
})

test_that("IV: multiple endogenous variables", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 | endo1 + endo2 | z1 + z2, data = d)
  expect_equal(result$K1, 2L)
  expect_true("endo1" %in% colnames(result$X))
  expect_true("endo2" %in% colnames(result$X))
  expect_false("endo1" %in% colnames(result$Z))
  expect_false("endo2" %in% colnames(result$Z))
})

test_that("IV: exact identification (L1 == K1)", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 | endo1 | z1, data = d)
  expect_equal(result$L1, result$K1)
  expect_false(result$is_overid)
  expect_equal(result$overid_df, 0L)
})

test_that("IV: overidentification (L1 > K1)", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 | endo1 | z1 + z2, data = d)
  expect_equal(result$L1, 2L)
  expect_equal(result$K1, 1L)
  expect_true(result$is_overid)
  expect_equal(result$overid_df, 1L)
})


# ==========================================================================
# 4. Duplicate detection
# ==========================================================================
test_that("duplicate between parts 1 and 2 is an error", {
  d <- make_test_data()
  expect_error(
    .parse_formula(y ~ x1 | x1 | z1, data = d),
    "Variable\\(s\\) appear in multiple formula parts: x1"
  )
})

test_that("duplicate between parts 1 and 3 is an error", {
  d <- make_test_data()
  expect_error(
    .parse_formula(y ~ x1 | endo1 | x1, data = d),
    "Variable\\(s\\) appear in multiple formula parts: x1"
  )
})

test_that("duplicate between parts 2 and 3 is an error", {
  d <- make_test_data()
  expect_error(
    .parse_formula(y ~ x1 | endo1 | endo1, data = d),
    "Variable\\(s\\) appear in multiple formula parts: endo1"
  )
})

test_that("multiple duplicates all named in error", {
  d <- make_test_data()
  expect_error(
    .parse_formula(y ~ x1 + x2 | x1 + x2 | z1, data = d),
    "x1.*x2|x2.*x1"
  )
})


# ==========================================================================
# 5. Collinearity detection
# ==========================================================================
test_that("collinear regressor dropped with warning naming variable", {
  d <- make_test_data()
  d$x2_dup <- d$x2
  expect_warning(
    result <- .parse_formula(y ~ x1 + x2 + x2_dup, data = d),
    "Dropped 1 collinear regressor: x2_dup"
  )
  expect_false("x2_dup" %in% colnames(result$X))
  expect_equal(result$dropped_regressors, "x2_dup")
  # Name vectors updated
  expect_false("x2_dup" %in% result$exog_names)
  expect_true("x2" %in% result$exog_names)
})

test_that("collinear instrument dropped with warning", {
  d <- make_test_data()
  d$z1_dup <- d$z1
  # Two warnings: "Dropped 1 collinear instrument" + "exact identification"
  expect_warning(
    expect_warning(
      result <- .parse_formula(y ~ x1 | endo1 | z1 + z1_dup, data = d),
      "Dropped 1 collinear instrument: z1_dup"
    ),
    "exact identification"
  )
  expect_false("z1_dup" %in% colnames(result$Z))
  expect_equal(result$dropped_instruments, "z1_dup")
  # Name vectors updated
  expect_false("z1_dup" %in% result$excluded_names)
  expect_true("z1" %in% result$excluded_names)
})

test_that("no warning when no collinearity", {
  d <- make_test_data()
  expect_no_warning(
    .parse_formula(y ~ x1 + x2, data = d)
  )
})


# ==========================================================================
# 6. Under-identification after collinearity
# ==========================================================================
test_that("error when dropping instruments causes underidentification", {
  d <- make_test_data()
  d$z1_dup <- d$z1
  # 2 endogenous, 2 instruments but one is collinear -> 1 instrument < 2 endo
  expect_error(
    suppressWarnings(
      .parse_formula(y ~ x1 | endo1 + endo2 | z1 + z1_dup, data = d)
    ),
    "underidentified"
  )
})

test_that("warning when dropping reduces to exact identification", {
  d <- make_test_data()
  # 1 endo, 2 instruments where one is collinear -> 1 instr = exact id
  d$z1_dup <- d$z1
  # Two warnings: "Dropped 1 collinear instrument" + "exact identification"
  expect_warning(
    expect_warning(
      result <- .parse_formula(y ~ x1 | endo1 | z1 + z1_dup, data = d),
      "Dropped 1 collinear instrument"
    ),
    "exact identification"
  )
  expect_false(result$is_overid)
})


# ==========================================================================
# 7. Intercept handling
# ==========================================================================
test_that("intercept present by default in X", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 + x2, data = d)
  expect_true("(Intercept)" %in% colnames(result$X))
  expect_true(result$has_intercept)
})

test_that("intercept present by default in Z for IV", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 | endo1 | z1, data = d)
  expect_true("(Intercept)" %in% colnames(result$Z))
})

test_that("0 + suppresses intercept in both X and Z", {
  d <- make_test_data()
  result <- .parse_formula(y ~ 0 + x1 | endo1 | z1, data = d)
  expect_false("(Intercept)" %in% colnames(result$X))
  expect_false("(Intercept)" %in% colnames(result$Z))
  expect_false(result$has_intercept)
})

test_that("part 2 intercept spec does not leak into X or Z", {
  d <- make_test_data()
  # Part 2 written as `1 + endo1` — the intercept should be stripped
  result <- .parse_formula(y ~ x1 | endo1 | z1, data = d)
  # Only one intercept in X (from part 1)
  expect_equal(sum(colnames(result$X) == "(Intercept)"), 1L)
  # Only one intercept in Z (from part 1)
  expect_equal(sum(colnames(result$Z) == "(Intercept)"), 1L)
})

test_that("no-intercept OLS works", {
  d <- make_test_data()
  result <- .parse_formula(y ~ 0 + x1 + x2, data = d)
  expect_false("(Intercept)" %in% colnames(result$X))
  expect_false(result$has_intercept)
  expect_equal(result$K, 2L)
})


# ==========================================================================
# 8. Weights / subset / na.action
# ==========================================================================
test_that("weights are passed through", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 + x2, data = d, weights = d$w)
  expect_equal(length(result$weights), 100L)
  expect_equal(result$weights, d$w)
})

test_that("subset reduces N", {
  d <- make_test_data()
  sub <- d$x1 > 0
  result <- .parse_formula(y ~ x1 + x2, data = d, subset = sub)
  expect_equal(result$N, sum(sub))
  expect_true(result$N < 100L)
})

test_that("NAs dropped by default", {
  d <- make_test_data()
  d$x1[1:5] <- NA
  result <- .parse_formula(y ~ x1 + x2, data = d)
  expect_equal(result$N, 95L)
  expect_false(is.null(result$na.action))
})


# ==========================================================================
# 9. Dimension counts
# ==========================================================================
test_that("OLS dimensions: N, K, K1=0, K2=K, L=K, L1=0", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 + x2 + x3, data = d)
  expect_equal(result$N, 100L)
  expect_equal(result$K, 4L)   # intercept + 3
  expect_equal(result$K1, 0L)
  expect_equal(result$K2, 4L)
  expect_equal(result$L, 4L)
  expect_equal(result$L1, 0L)
})

test_that("IV just-identified dimensions", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 + x2 | endo1 | z1, data = d)
  expect_equal(result$N, 100L)
  expect_equal(result$K, 4L)    # intercept + x1 + x2 + endo1
  expect_equal(result$K1, 1L)   # endo1
  expect_equal(result$K2, 3L)   # intercept + x1 + x2
  expect_equal(result$L, 4L)    # intercept + x1 + x2 + z1
  expect_equal(result$L1, 1L)   # z1
  expect_equal(result$overid_df, 0L)
})

test_that("IV overidentified dimensions", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 | endo1 + endo2 | z1 + z2 + z3, data = d)
  expect_equal(result$K, 4L)    # intercept + x1 + endo1 + endo2
  expect_equal(result$K1, 2L)   # endo1 + endo2
  expect_equal(result$K2, 2L)   # intercept + x1
  expect_equal(result$L, 5L)    # intercept + x1 + z1 + z2 + z3
  expect_equal(result$L1, 3L)   # z1 + z2 + z3
  expect_equal(result$overid_df, 1L)
})


# ==========================================================================
# 10. Integration with Card data
# ==========================================================================
test_that("Card data: just-identified IV spec matches expected dimensions", {
  card_path <- test_path("..", "stata-benchmarks", "fixtures", "card_data.csv")
  skip_if(!file.exists(card_path), "card data not found")
  card <- read.csv(card_path)
  result <- .parse_formula(
    lwage ~ exper + expersq + black + south | educ | nearc4,
    data = card
  )
  expect_equal(result$N, 3010L)
  expect_equal(result$K, 6L)    # intercept + exper + expersq + black + south + educ
  expect_equal(result$K1, 1L)   # educ
  expect_equal(result$K2, 5L)   # intercept + 4 exogenous
  expect_equal(result$L, 6L)    # intercept + exper + expersq + black + south + nearc4
  expect_equal(result$L1, 1L)   # nearc4
  expect_true(result$is_iv)
  expect_false(result$is_overid)
  expect_equal(result$overid_df, 0L)
  expect_equal(result$y_name, "lwage")
  expect_equal(result$endo_names, "educ")
  expect_equal(result$excluded_names, "nearc4")
  expect_true(result$has_intercept)
})

test_that("Card data: OLS spec produces correct structure", {
  card_path <- test_path("..", "stata-benchmarks", "fixtures", "card_data.csv")
  skip_if(!file.exists(card_path), "card data not found")
  card <- read.csv(card_path)
  result <- .parse_formula(
    lwage ~ exper + expersq + black + south + educ,
    data = card
  )
  expect_equal(result$N, 3010L)
  expect_equal(result$K, 6L)
  expect_null(result$Z)
  expect_false(result$is_iv)
})

test_that("Card data: overidentified IV spec", {
  card_path <- test_path("..", "stata-benchmarks", "fixtures", "card_data.csv")
  skip_if(!file.exists(card_path), "card data not found")
  card <- read.csv(card_path)
  result <- .parse_formula(
    lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
    data = card
  )
  expect_equal(result$N, 3010L)
  expect_equal(result$K1, 1L)
  expect_equal(result$L1, 2L)
  expect_true(result$is_overid)
  expect_equal(result$overid_df, 1L)
})

test_that("Card data: response vector matches data column", {
  card_path <- test_path("..", "stata-benchmarks", "fixtures", "card_data.csv")
  skip_if(!file.exists(card_path), "card data not found")
  card <- read.csv(card_path)
  result <- .parse_formula(lwage ~ educ + exper, data = card)
  expect_equal(as.numeric(result$y), card$lwage)
})

test_that("Card data: variable names in return match formula", {
  card_path <- test_path("..", "stata-benchmarks", "fixtures", "card_data.csv")
  skip_if(!file.exists(card_path), "card data not found")
  card <- read.csv(card_path)
  result <- .parse_formula(
    lwage ~ exper + expersq + black + south | educ | nearc4,
    data = card
  )
  expect_equal(result$exog_names, c("exper", "expersq", "black", "south"))
  expect_equal(result$endo_names, "educ")
  expect_equal(result$excluded_names, "nearc4")
})

# ==========================================================================
# Edge cases: response validation
# ==========================================================================
test_that("no-response formula is rejected", {
  d <- make_test_data()
  expect_error(
    .parse_formula(~ x1 + x2, data = d),
    "Formula must have a response variable"
  )
})

test_that("multivariate response (cbind) is rejected", {
  d <- make_test_data()
  expect_error(
    .parse_formula(cbind(y, x1) ~ x2, data = d),
    "Multivariate responses.*not supported"
  )
})

# ==========================================================================
# Edge cases: terms$regressors includes endogenous variables
# ==========================================================================
test_that("IV: terms$regressors includes endogenous term labels", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 | endo1 | z1, data = d)
  reg_labels <- attr(result$terms$regressors, "term.labels")
  expect_true("endo1" %in% reg_labels)
  expect_true("x1" %in% reg_labels)
})

test_that("IV: terms$regressors is usable with model.matrix", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 | endo1 | z1, data = d)
  # Should not error — this failed when "(Intercept)" leaked into term labels
  mm <- model.matrix(result$terms$regressors, d)
  expect_true("(Intercept)" %in% colnames(mm))
  expect_true("x1" %in% colnames(mm))
  expect_true("endo1" %in% colnames(mm))
  expect_false("Intercept" %in% attr(result$terms$regressors, "term.labels"))
})

test_that("OLS: terms$regressors matches exogenous terms", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 + x2, data = d)
  reg_labels <- attr(result$terms$regressors, "term.labels")
  expect_equal(reg_labels, c("x1", "x2"))
})

# ==========================================================================
# Edge cases: general
# ==========================================================================
test_that("factor variable in exogenous regressors expands correctly", {
  d <- make_test_data()
  d$grp <- factor(d$grp)
  result <- .parse_formula(y ~ x1 + grp, data = d)
  # Factor should expand to K-1 dummies
  expect_true(any(grepl("^grp", colnames(result$X))))
  # Original variable name is in exog_names
  expect_true("grp" %in% result$exog_names)
})

test_that("return list has all expected fields", {
  d <- make_test_data()
  result <- .parse_formula(y ~ x1 | endo1 | z1, data = d)
  expected_names <- c(
    "y", "X", "Z", "y_name", "exog_names", "endo_names", "excluded_names",
    "X_names", "Z_names", "N", "K", "K1", "K2", "L", "L1",
    "has_intercept", "model_frame", "terms", "formula", "weights",
    "na.action", "dropped_regressors", "dropped_instruments",
    "reclassified_endogenous", "is_iv", "is_overid", "overid_df"
  )
  expect_true(all(expected_names %in% names(result)))
})
