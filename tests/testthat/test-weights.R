# ============================================================================
# Tests: Analytic Weights (Ticket C3)
# ============================================================================

# --- Helpers ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)

read_vcov_fixture <- function(path) {
  fixture <- read.csv(path)
  stata_names <- fixture$term
  r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)
  vcov_cols <- grep("^vcov_", names(fixture), value = TRUE)
  V <- as.matrix(fixture[, vcov_cols])
  rownames(V) <- r_names
  col_stata <- sub("^vcov_", "", vcov_cols)
  colnames(V) <- ifelse(col_stata == "_cons", "(Intercept)", col_stata)
  V
}

# --- Load Card data ---
card_path <- file.path(fixture_dir, "card_data.csv")
if (file.exists(card_path)) {
  card <- read.csv(card_path)
}


# ============================================================================
# Weighted OLS vs lm(): exact match (iid VCV)
# ============================================================================

test_that("weighted OLS coefficients match lm() with weights", {
  fit_iv2 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp)
  fit_lm <- lm(mpg ~ wt + hp, data = mtcars, weights = disp)
  expect_equal(coef(fit_iv2), coef(fit_lm),
               tolerance = .Machine$double.eps^0.5)
})

test_that("weighted OLS SEs match lm() with weights (iid, small=TRUE)", {
  fit_iv2 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
                    small = TRUE)
  fit_lm <- lm(mpg ~ wt + hp, data = mtcars, weights = disp)
  se_iv2 <- sqrt(diag(fit_iv2$vcov))
  se_lm <- sqrt(diag(vcov(fit_lm)))
  expect_equal(se_iv2, se_lm, tolerance = .Machine$double.eps^0.5)
})

test_that("weighted OLS sigma is scale-invariant in weights", {
  # Stata aweight convention: normalizing weights to sum to N makes sigma
  # independent of weight scale. This differs from lm() where sigma scales
  # with sqrt(c) when weights are multiplied by c.
  fit1 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp, small = TRUE)
  d2 <- mtcars
  d2$disp10 <- mtcars$disp * 10
  fit2 <- ivreg2(mpg ~ wt + hp, data = d2, weights = disp10, small = TRUE)
  expect_equal(fit1$sigma, fit2$sigma, tolerance = .Machine$double.eps^0.5)
})

test_that("weighted OLS sigma matches Stata convention (differs from lm by sqrt(N/sum(w)))", {
  fit_iv2 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
                    small = TRUE)
  fit_lm <- lm(mpg ~ wt + hp, data = mtcars, weights = disp)
  N <- nrow(mtcars)
  correction <- sqrt(N / sum(mtcars$disp))
  expect_equal(fit_iv2$sigma, summary(fit_lm)$sigma * correction,
               tolerance = .Machine$double.eps^0.5)
})

test_that("weighted OLS R-squared matches lm() summary", {
  fit_iv2 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
                    small = TRUE)
  fit_lm <- lm(mpg ~ wt + hp, data = mtcars, weights = disp)
  expect_equal(fit_iv2$r.squared, summary(fit_lm)$r.squared,
               tolerance = .Machine$double.eps^0.5)
  expect_equal(fit_iv2$adj.r.squared, summary(fit_lm)$adj.r.squared,
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Weighted OLS vs lm(): HC1 VCV
# ============================================================================

test_that("weighted OLS HC1 VCV matches sandwich::vcovHC with lm weights", {
  skip_if_not_installed("sandwich")
  fit_iv2 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
                    vcov = "HC1")
  fit_lm <- lm(mpg ~ wt + hp, data = mtcars, weights = disp)
  V_sand <- sandwich::vcovHC(fit_lm, type = "HC1")
  expect_equal(fit_iv2$vcov, V_sand, tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Weighted 2SLS vs Stata: coefficients (iid, small=FALSE and TRUE)
# ============================================================================

test_that("weighted 2SLS coefficients match Stata fixture (iid)", {
  skip_if(!file.exists(card_path), "Card data not found")
  coef_path <- file.path(fixture_dir, "card_just_id_weighted_coef_iid.csv")
  skip_if(!file.exists(coef_path), "Fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight)
  fixture <- read.csv(coef_path)
  fixture$r_name <- ifelse(fixture$term == "_cons", "(Intercept)", fixture$term)

  for (i in seq_len(nrow(fixture))) {
    nm <- fixture$r_name[i]
    expect_equal(
      unname(coef(fit)[nm]), fixture$estimate[i],
      tolerance = stata_tol$coef,
      info = paste("Coef mismatch:", nm)
    )
  }
})

test_that("weighted 2SLS SEs match Stata fixture (iid, small=FALSE)", {
  skip_if(!file.exists(card_path), "Card data not found")
  coef_path <- file.path(fixture_dir, "card_just_id_weighted_coef_iid.csv")
  skip_if(!file.exists(coef_path), "Fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight)
  fixture <- read.csv(coef_path)
  fixture$r_name <- ifelse(fixture$term == "_cons", "(Intercept)", fixture$term)

  for (i in seq_len(nrow(fixture))) {
    nm <- fixture$r_name[i]
    expect_equal(
      sqrt(fit$vcov[nm, nm]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", nm)
    )
  }
})

test_that("weighted 2SLS SEs match Stata fixture (iid, small=TRUE)", {
  skip_if(!file.exists(card_path), "Card data not found")
  coef_path <- file.path(fixture_dir, "card_just_id_weighted_coef_iid_small.csv")
  skip_if(!file.exists(coef_path), "Fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight, small = TRUE)
  fixture <- read.csv(coef_path)
  fixture$r_name <- ifelse(fixture$term == "_cons", "(Intercept)", fixture$term)

  for (i in seq_len(nrow(fixture))) {
    nm <- fixture$r_name[i]
    expect_equal(
      sqrt(fit$vcov[nm, nm]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", nm)
    )
  }
})


# ============================================================================
# Weighted 2SLS vs Stata: robust VCV
# ============================================================================

test_that("weighted 2SLS SEs match Stata fixture (robust/HC0, small=FALSE)", {
  skip_if(!file.exists(card_path), "Card data not found")
  coef_path <- file.path(fixture_dir, "card_just_id_weighted_coef_hc1.csv")
  skip_if(!file.exists(coef_path), "Fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight, vcov = "HC0")
  fixture <- read.csv(coef_path)
  fixture$r_name <- ifelse(fixture$term == "_cons", "(Intercept)", fixture$term)

  for (i in seq_len(nrow(fixture))) {
    nm <- fixture$r_name[i]
    expect_equal(
      sqrt(fit$vcov[nm, nm]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", nm)
    )
  }
})

test_that("weighted 2SLS SEs match Stata fixture (robust/HC1, small=TRUE)", {
  skip_if(!file.exists(card_path), "Card data not found")
  coef_path <- file.path(fixture_dir, "card_just_id_weighted_coef_hc1_small.csv")
  skip_if(!file.exists(coef_path), "Fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight, vcov = "HC1", small = TRUE)
  fixture <- read.csv(coef_path)
  fixture$r_name <- ifelse(fixture$term == "_cons", "(Intercept)", fixture$term)

  for (i in seq_len(nrow(fixture))) {
    nm <- fixture$r_name[i]
    expect_equal(
      sqrt(fit$vcov[nm, nm]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", nm)
    )
  }
})


# ============================================================================
# Weighted 2SLS vs Stata: Cluster VCV
# ============================================================================

test_that("weighted 2SLS SEs match Stata fixture (cluster, small=FALSE)", {
  skip_if(!file.exists(card_path), "Card data not found")
  coef_path <- file.path(fixture_dir, "card_just_id_weighted_coef_cl.csv")
  skip_if(!file.exists(coef_path), "Fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight, clusters = ~smsa66)
  fixture <- read.csv(coef_path)
  fixture$r_name <- ifelse(fixture$term == "_cons", "(Intercept)", fixture$term)

  for (i in seq_len(nrow(fixture))) {
    nm <- fixture$r_name[i]
    expect_equal(
      sqrt(fit$vcov[nm, nm]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", nm)
    )
  }
})

test_that("weighted 2SLS SEs match Stata fixture (cluster, small=TRUE)", {
  skip_if(!file.exists(card_path), "Card data not found")
  coef_path <- file.path(fixture_dir, "card_just_id_weighted_coef_cl_small.csv")
  skip_if(!file.exists(coef_path), "Fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight, clusters = ~smsa66,
                small = TRUE)
  fixture <- read.csv(coef_path)
  fixture$r_name <- ifelse(fixture$term == "_cons", "(Intercept)", fixture$term)

  for (i in seq_len(nrow(fixture))) {
    nm <- fixture$r_name[i]
    expect_equal(
      sqrt(fit$vcov[nm, nm]), fixture$std_error[i],
      tolerance = stata_tol$se,
      info = paste("SE mismatch:", nm)
    )
  }
})


# ============================================================================
# Weighted 2SLS vs Stata: full VCV matrix
# ============================================================================

test_that("weighted 2SLS VCV matches Stata fixture (iid)", {
  skip_if(!file.exists(card_path), "Card data not found")
  vcov_path <- file.path(fixture_dir, "card_just_id_weighted_vcov_iid.csv")
  skip_if(!file.exists(vcov_path), "VCV fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight)
  V_stata <- read_vcov_fixture(vcov_path)
  shared <- intersect(rownames(fit$vcov), rownames(V_stata))
  for (rn in shared) {
    for (cn in shared) {
      expect_equal(
        fit$vcov[rn, cn], V_stata[rn, cn],
        tolerance = stata_tol$vcov,
        info = paste("VCV mismatch:", rn, cn)
      )
    }
  }
})


# ============================================================================
# Weighted 2SLS vs Stata: sigma (RMSE) and RSS
# ============================================================================

test_that("weighted 2SLS sigma matches Stata RMSE fixture", {
  skip_if(!file.exists(card_path), "Card data not found")
  diag_path <- file.path(fixture_dir, "card_just_id_weighted_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "Diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight)
  fixture <- read.csv(diag_path)
  expect_equal(fit$sigma, fixture$rmse, tolerance = stata_tol$se)
})

test_that("weighted 2SLS RSS matches Stata RSS fixture", {
  skip_if(!file.exists(card_path), "Card data not found")
  diag_path <- file.path(fixture_dir, "card_just_id_weighted_diagnostics_iid.csv")
  skip_if(!file.exists(diag_path), "Diagnostics fixture not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight)
  fixture <- read.csv(diag_path)
  expect_equal(fit$rss, fixture$rss, tolerance = stata_tol$stat)
})


# ============================================================================
# Weight validation
# ============================================================================

test_that("negative weights give error", {
  d <- mtcars
  d$w <- c(-1, rep(1, nrow(d) - 1))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w),
    "strictly positive"
  )
})

test_that("zero weights give error", {
  d <- mtcars
  d$w <- c(0, rep(1, nrow(d) - 1))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w),
    "strictly positive"
  )
})

test_that("NA weights give error", {
  d <- mtcars
  d$w <- c(NA, rep(1, nrow(d) - 1))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w, na.action = na.pass),
    "finite and non-missing"
  )
})

test_that("Inf weights give error", {
  d <- mtcars
  d$w <- c(Inf, rep(1, nrow(d) - 1))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w),
    "finite and non-missing"
  )
})

test_that("NULL weights produce unweighted results", {
  fit_null <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit_null$weights)
  expect_equal(coef(fit_null), coef(lm(mpg ~ wt + hp, data = mtcars)),
               tolerance = .Machine$double.eps^0.5)
})


# ============================================================================
# Edge cases
# ============================================================================

test_that("equal weights produce same results as unweighted", {
  d <- mtcars
  d$w <- rep(1, nrow(d))
  fit_w <- ivreg2(mpg ~ wt + hp, data = d, weights = w, small = TRUE)
  fit_u <- ivreg2(mpg ~ wt + hp, data = d, small = TRUE)
  expect_equal(coef(fit_w), coef(fit_u),
               tolerance = .Machine$double.eps^0.5)
  expect_equal(fit_w$vcov, fit_u$vcov,
               tolerance = .Machine$double.eps^0.5)
  expect_equal(fit_w$sigma, fit_u$sigma,
               tolerance = .Machine$double.eps^0.5)
  expect_equal(fit_w$r.squared, fit_u$r.squared,
               tolerance = .Machine$double.eps^0.5)
})

test_that("weighted OLS 1-part formula works", {
  fit <- ivreg2(mpg ~ wt + hp + cyl, data = mtcars, weights = disp)
  expect_s3_class(fit, "ivreg2")
  expect_length(coef(fit), 4L)
  expect_equal(fit$weights, mtcars$disp)
})

test_that("weighted IV with clustering works", {
  skip_if(!file.exists(card_path), "Card data not found")
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight, clusters = ~smsa66)
  expect_s3_class(fit, "ivreg2")
  expect_identical(fit$vcov_type, "CL")
  expect_false(is.null(fit$weights))
})

test_that("weights are stored in the returned object", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp)
  expect_equal(fit$weights, mtcars$disp)
})
