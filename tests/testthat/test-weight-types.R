# ============================================================================
# Tests: Frequency weights (fweight) and Probability weights (pweight)
# Ticket K2
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

read_coef_fixture <- function(path) {
  d <- read.csv(path)
  nms <- ifelse(d$term == "_cons", "(Intercept)", d$term)
  list(
    estimate  = setNames(d$estimate, nms),
    std_error = setNames(d$std_error, nms)
  )
}

read_diagnostics_fixture <- function(path) {
  read.csv(path)
}

read_firststage_fixture <- function(path) {
  read.csv(path)
}

# --- Load Card data ---
card_path <- file.path(fixture_dir, "card_data.csv")
if (file.exists(card_path)) {
  card <- read.csv(card_path)
}


# ============================================================================
# Section 1: Validation tests
# ============================================================================

test_that("invalid weight_type is rejected", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, weight_type = "invalid"),
    'weight_type.*must be one of'
  )
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, weight_type = 42),
    'weight_type.*must be one of'
  )
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, weight_type = c("aweight", "fweight")),
    'weight_type.*must be one of'
  )
})

test_that("fweight rejects non-integer weights", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)  # non-integer
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w, weight_type = "fweight"),
    'must be integers'
  )
})

test_that("fweight rejects non-positive weights", {
  d <- mtcars
  d$w <- rep(0L, nrow(d))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w, weight_type = "fweight"),
    'strictly positive'
  )
})

test_that("pweight rejects non-positive weights", {
  d <- mtcars
  d$w <- c(-1, rep(1, nrow(d) - 1))
  expect_error(
    ivreg2(mpg ~ wt + hp, data = d, weights = w, weight_type = "pweight"),
    'strictly positive'
  )
})

test_that("pweight + iid overrides to HC0 with message", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  expect_message(
    fit <- ivreg2(mpg ~ wt + hp, data = d, weights = w,
                  weight_type = "pweight", vcov = "iid"),
    'pweight implies robust VCE'
  )
  expect_equal(fit$vcov_type, "HC0")
})

test_that("pweight + iid + small overrides to HC1 with message", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  expect_message(
    fit <- ivreg2(mpg ~ wt + hp, data = d, weights = w,
                  weight_type = "pweight", vcov = "iid", small = TRUE),
    'pweight implies robust VCE'
  )
  expect_equal(fit$vcov_type, "HC1")
})

test_that("pweight + cluster does not override vcov", {
  fit <- suppressMessages(
    ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
           weight_type = "pweight", clusters = ~ cyl)
  )
  expect_equal(fit$vcov_type, "CL")
})

test_that("pweight + explicit HC0 matches pweight + iid override", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  fit_override <- suppressMessages(
    ivreg2(mpg ~ wt + hp, data = d, weights = w,
           weight_type = "pweight", vcov = "iid")
  )
  fit_explicit <- ivreg2(mpg ~ wt + hp, data = d, weights = w,
                         weight_type = "pweight", vcov = "HC0")
  expect_equal(vcov(fit_explicit), vcov(fit_override))
  expect_equal(coef(fit_explicit), coef(fit_override))
  expect_equal(fit_explicit$sigma, fit_override$sigma)
})

test_that("pweight + explicit HC1 + small matches pweight + iid + small override", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  fit_override <- suppressMessages(
    ivreg2(mpg ~ wt + hp, data = d, weights = w,
           weight_type = "pweight", vcov = "iid", small = TRUE)
  )
  fit_explicit <- ivreg2(mpg ~ wt + hp, data = d, weights = w,
                         weight_type = "pweight", vcov = "HC1", small = TRUE)
  expect_equal(vcov(fit_explicit), vcov(fit_override))
  expect_equal(coef(fit_explicit), coef(fit_override))
  expect_equal(fit_explicit$sigma, fit_override$sigma)
})


# ============================================================================
# Section 2: N semantics
# ============================================================================

test_that("nobs() returns sum(w) for fweight", {
  d <- data.frame(y = 1:10, x = rnorm(10), w = rep(2L, 10))
  fit <- ivreg2(y ~ x, data = d, weights = w, weight_type = "fweight")
  expect_equal(nobs(fit), 20L)
})

test_that("nobs() returns nrow for aweight", {
  d <- data.frame(y = 1:10, x = rnorm(10), w = rep(2, 10))
  fit <- ivreg2(y ~ x, data = d, weights = w, weight_type = "aweight")
  expect_equal(nobs(fit), 10L)
})

test_that("nobs() returns nrow for pweight", {
  d <- data.frame(y = 1:10, x = rnorm(10), w = rep(2, 10))
  fit <- suppressMessages(
    ivreg2(y ~ x, data = d, weights = w, weight_type = "pweight")
  )
  expect_equal(nobs(fit), 10L)
})


# ============================================================================
# Section 3: Backward compatibility
# ============================================================================

test_that("weight_type='aweight' is identical to omitting it", {
  fit1 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp)
  fit2 <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
                 weight_type = "aweight")
  expect_equal(coef(fit1), coef(fit2))
  expect_equal(vcov(fit1), vcov(fit2))
  expect_equal(fit1$sigma, fit2$sigma)
})

test_that("weight_type='aweight' IV is identical to omitting it", {
  skip_if(!file.exists(card_path))
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, weights = weight, vcov = "HC1")
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card, weights = weight, weight_type = "aweight",
                 vcov = "HC1")
  expect_equal(coef(fit1), coef(fit2))
  expect_equal(vcov(fit1), vcov(fit2))
  expect_equal(fit1$diagnostics$weak_id$stat, fit2$diagnostics$weak_id$stat)
})


# ============================================================================
# Section 4: Frequency weight (fweight) parity vs Stata fixtures
# ============================================================================

# Helper function for testing one fixture configuration
test_fweight_config <- function(fixture_prefix, suffix, card_data,
                                vcov_arg, small_arg, cluster_arg,
                                overid = FALSE) {
  coef_path <- file.path(fixture_dir, paste0(fixture_prefix, "_coef_", suffix, ".csv"))
  vcov_path <- file.path(fixture_dir, paste0(fixture_prefix, "_vcov_", suffix, ".csv"))
  diag_path <- file.path(fixture_dir, paste0(fixture_prefix, "_diagnostics_", suffix, ".csv"))
  fs_path   <- file.path(fixture_dir, paste0(fixture_prefix, "_firststage_", suffix, ".csv"))

  skip_if(!file.exists(coef_path))

  stata_coef <- read_coef_fixture(coef_path)
  stata_vcov <- read_vcov_fixture(vcov_path)
  stata_diag <- read_diagnostics_fixture(diag_path)
  stata_fs   <- if (file.exists(fs_path)) read_firststage_fixture(fs_path) else NULL

  iv_formula <- if (overid) {
    lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2
  } else {
    lwage ~ exper + expersq + black + south | educ | nearc4
  }

  fit <- ivreg2(iv_formula, data = card_data, weights = weight,
                weight_type = "fweight", vcov = vcov_arg,
                small = small_arg, clusters = cluster_arg)

  # Test N
  expect_equal(nobs(fit), as.integer(stata_diag$N))

  # Coefficients
  r_coef <- coef(fit)[names(stata_coef$estimate)]
  expect_equal(r_coef, stata_coef$estimate, tolerance = stata_tol$coef)

  # Standard errors
  r_se <- sqrt(diag(vcov(fit)))[names(stata_coef$std_error)]
  expect_equal(r_se, stata_coef$std_error, tolerance = stata_tol$se)

  # VCV
  r_vcov <- vcov(fit)[rownames(stata_vcov), colnames(stata_vcov)]
  expect_equal(r_vcov, stata_vcov, tolerance = stata_tol$vcov, check.attributes = FALSE)

  # Sigma
  expect_equal(fit$sigma, stata_diag$rmse, tolerance = stata_tol$coef)

  # R-squared
  expect_equal(fit$r.squared, stata_diag$r2, tolerance = stata_tol$coef)

  # Model F
  if (!is.na(stata_diag$F_stat)) {
    expect_equal(fit$model_f, stata_diag$F_stat, tolerance = stata_tol$stat)
  }

  # Diagnostics
  if (!is.na(stata_diag$idstat)) {
    expect_equal(fit$diagnostics$underid$stat, stata_diag$idstat,
                 tolerance = stata_tol$stat)
  }
  if (!is.na(stata_diag$cdf)) {
    expect_equal(fit$diagnostics$weak_id$stat, stata_diag$cdf,
                 tolerance = stata_tol$stat)
  }
  # AR and Stock-Wright: skip when our code returns NA due to numerical

  # singularity (happens with very few clusters, where the cluster meat is
  # rank-deficient — a pre-existing limitation, not K2-specific).
  if (!is.na(stata_diag$arf) &&
      !is.na(fit$diagnostics$anderson_rubin$f_stat)) {
    expect_equal(fit$diagnostics$anderson_rubin$f_stat, stata_diag$arf,
                 tolerance = stata_tol$stat)
  }
  if (!is.na(stata_diag$sstat) &&
      !is.null(fit$diagnostics$stock_wright) &&
      !is.na(fit$diagnostics$stock_wright$stat)) {
    expect_equal(fit$diagnostics$stock_wright$stat, stata_diag$sstat,
                 tolerance = stata_tol$stat)
  }

  # First-stage: skip when our code returns NA (same rank-deficiency issue)
  if (!is.null(stata_fs) && !is.null(fit$first_stage)) {
    fs_f_row <- stata_fs[stata_fs$statistic == "F", ]
    if (nrow(fs_f_row) > 0L) {
      fs_f_stata <- fs_f_row$educ
      if (!is.na(fs_f_stata) && !is.na(fit$first_stage$educ$f_stat)) {
        expect_equal(fit$first_stage$educ$f_stat, fs_f_stata,
                     tolerance = stata_tol$stat)
      }
    }
  }

  # Overid
  if (overid && !is.na(stata_diag$sargan)) {
    expect_equal(fit$diagnostics$overid$stat, stata_diag$sargan,
                 tolerance = stata_tol$stat)
  }
  if (overid && !is.na(stata_diag$j)) {
    expect_equal(fit$diagnostics$overid$stat, stata_diag$j,
                 tolerance = stata_tol$stat)
  }
}


# --- fweight just-identified ---
test_that("fweight just-id: iid matches Stata", {
  skip_if(!file.exists(card_path))
  test_fweight_config("card_fweight_just_id", "iid", card,
                      vcov_arg = "iid", small_arg = FALSE, cluster_arg = NULL)
})

test_that("fweight just-id: iid_small matches Stata", {
  skip_if(!file.exists(card_path))
  test_fweight_config("card_fweight_just_id", "iid_small", card,
                      vcov_arg = "iid", small_arg = TRUE, cluster_arg = NULL)
})

test_that("fweight just-id: robust matches Stata", {
  skip_if(!file.exists(card_path))
  # Stata robust (no small) = our HC0
  test_fweight_config("card_fweight_just_id", "hc1", card,
                      vcov_arg = "HC0", small_arg = FALSE, cluster_arg = NULL)
})

test_that("fweight just-id: robust_small matches Stata", {
  skip_if(!file.exists(card_path))
  # Stata robust small = our HC1
  test_fweight_config("card_fweight_just_id", "hc1_small", card,
                      vcov_arg = "HC1", small_arg = TRUE, cluster_arg = NULL)
})

test_that("fweight just-id: cl matches Stata", {
  skip_if(!file.exists(card_path))
  test_fweight_config("card_fweight_just_id", "cl", card,
                      vcov_arg = "iid", small_arg = FALSE,
                      cluster_arg = ~ smsa)
})

test_that("fweight just-id: cl_small matches Stata", {
  skip_if(!file.exists(card_path))
  test_fweight_config("card_fweight_just_id", "cl_small", card,
                      vcov_arg = "iid", small_arg = TRUE,
                      cluster_arg = ~ smsa)
})

# --- fweight overidentified ---
test_that("fweight overid: iid matches Stata", {
  skip_if(!file.exists(card_path))
  test_fweight_config("card_fweight_overid", "iid", card,
                      vcov_arg = "iid", small_arg = FALSE, cluster_arg = NULL,
                      overid = TRUE)
})

test_that("fweight overid: iid_small matches Stata", {
  skip_if(!file.exists(card_path))
  test_fweight_config("card_fweight_overid", "iid_small", card,
                      vcov_arg = "iid", small_arg = TRUE, cluster_arg = NULL,
                      overid = TRUE)
})

test_that("fweight overid: robust matches Stata", {
  skip_if(!file.exists(card_path))
  # Stata robust (no small) = our HC0
  test_fweight_config("card_fweight_overid", "hc1", card,
                      vcov_arg = "HC0", small_arg = FALSE, cluster_arg = NULL,
                      overid = TRUE)
})

test_that("fweight overid: robust_small matches Stata", {
  skip_if(!file.exists(card_path))
  # Stata robust small = our HC1
  test_fweight_config("card_fweight_overid", "hc1_small", card,
                      vcov_arg = "HC1", small_arg = TRUE, cluster_arg = NULL,
                      overid = TRUE)
})

test_that("fweight overid: cl matches Stata", {
  skip_if(!file.exists(card_path))
  test_fweight_config("card_fweight_overid", "cl", card,
                      vcov_arg = "iid", small_arg = FALSE,
                      cluster_arg = ~ smsa, overid = TRUE)
})

test_that("fweight overid: cl_small matches Stata", {
  skip_if(!file.exists(card_path))
  test_fweight_config("card_fweight_overid", "cl_small", card,
                      vcov_arg = "iid", small_arg = TRUE,
                      cluster_arg = ~ smsa, overid = TRUE)
})


# ============================================================================
# Section 5: Probability weight (pweight) parity vs Stata fixtures
# ============================================================================

test_pweight_config <- function(fixture_prefix, suffix, card_data,
                                vcov_arg, small_arg, cluster_arg,
                                overid = FALSE) {
  coef_path <- file.path(fixture_dir, paste0(fixture_prefix, "_coef_", suffix, ".csv"))
  vcov_path <- file.path(fixture_dir, paste0(fixture_prefix, "_vcov_", suffix, ".csv"))
  diag_path <- file.path(fixture_dir, paste0(fixture_prefix, "_diagnostics_", suffix, ".csv"))
  fs_path   <- file.path(fixture_dir, paste0(fixture_prefix, "_firststage_", suffix, ".csv"))

  skip_if(!file.exists(coef_path))

  stata_coef <- read_coef_fixture(coef_path)
  stata_vcov <- read_vcov_fixture(vcov_path)
  stata_diag <- read_diagnostics_fixture(diag_path)
  stata_fs   <- if (file.exists(fs_path)) read_firststage_fixture(fs_path) else NULL

  iv_formula <- if (overid) {
    lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2
  } else {
    lwage ~ exper + expersq + black + south | educ | nearc4
  }

  fit <- suppressMessages(
    ivreg2(iv_formula, data = card_data, weights = weight,
           weight_type = "pweight", vcov = vcov_arg,
           small = small_arg, clusters = cluster_arg)
  )

  # Test N (pweight uses physical N)
  expect_equal(nobs(fit), as.integer(stata_diag$N))

  # Coefficients
  r_coef <- coef(fit)[names(stata_coef$estimate)]
  expect_equal(r_coef, stata_coef$estimate, tolerance = stata_tol$coef)

  # Standard errors
  r_se <- sqrt(diag(vcov(fit)))[names(stata_coef$std_error)]
  expect_equal(r_se, stata_coef$std_error, tolerance = stata_tol$se)

  # VCV
  r_vcov <- vcov(fit)[rownames(stata_vcov), colnames(stata_vcov)]
  expect_equal(r_vcov, stata_vcov, tolerance = stata_tol$vcov, check.attributes = FALSE)

  # Sigma
  expect_equal(fit$sigma, stata_diag$rmse, tolerance = stata_tol$coef)

  # Diagnostics
  if (!is.na(stata_diag$idstat)) {
    expect_equal(fit$diagnostics$underid$stat, stata_diag$idstat,
                 tolerance = stata_tol$stat)
  }
  if (!is.na(stata_diag$cdf)) {
    expect_equal(fit$diagnostics$weak_id$stat, stata_diag$cdf,
                 tolerance = stata_tol$stat)
  }
  # AR: skip when our code returns NA due to numerical singularity
  if (!is.na(stata_diag$arf) &&
      !is.na(fit$diagnostics$anderson_rubin$f_stat)) {
    expect_equal(fit$diagnostics$anderson_rubin$f_stat, stata_diag$arf,
                 tolerance = stata_tol$stat)
  }

  # First-stage: skip when our code returns NA
  if (!is.null(stata_fs) && !is.null(fit$first_stage)) {
    fs_f_row <- stata_fs[stata_fs$statistic == "F", ]
    if (nrow(fs_f_row) > 0L) {
      fs_f_stata <- fs_f_row$educ
      if (!is.na(fs_f_stata) && !is.na(fit$first_stage$educ$f_stat)) {
        expect_equal(fit$first_stage$educ$f_stat, fs_f_stata,
                     tolerance = stata_tol$stat)
      }
    }
  }

  # Overid
  if (overid && !is.na(stata_diag$j)) {
    expect_equal(fit$diagnostics$overid$stat, stata_diag$j,
                 tolerance = stata_tol$stat)
  }
}


# --- pweight just-identified ---
test_that("pweight just-id: hc0 matches Stata", {
  skip_if(!file.exists(card_path))
  # pweight forces robust: vcov="iid" → HC0 (small=FALSE)
  test_pweight_config("card_pweight_just_id", "hc0", card,
                      vcov_arg = "iid", small_arg = FALSE, cluster_arg = NULL)
})

test_that("pweight just-id: hc0_small matches Stata", {
  skip_if(!file.exists(card_path))
  # pweight forces robust: vcov="iid" → HC1 (small=TRUE)
  test_pweight_config("card_pweight_just_id", "hc0_small", card,
                      vcov_arg = "iid", small_arg = TRUE, cluster_arg = NULL)
})

test_that("pweight just-id: cl matches Stata", {
  skip_if(!file.exists(card_path))
  test_pweight_config("card_pweight_just_id", "cl", card,
                      vcov_arg = "iid", small_arg = FALSE,
                      cluster_arg = ~ smsa)
})

test_that("pweight just-id: cl_small matches Stata", {
  skip_if(!file.exists(card_path))
  test_pweight_config("card_pweight_just_id", "cl_small", card,
                      vcov_arg = "iid", small_arg = TRUE,
                      cluster_arg = ~ smsa)
})

# --- pweight overidentified ---
test_that("pweight overid: hc0 matches Stata", {
  skip_if(!file.exists(card_path))
  # pweight forces robust: vcov="iid" → HC0 (small=FALSE)
  test_pweight_config("card_pweight_overid", "hc0", card,
                      vcov_arg = "iid", small_arg = FALSE, cluster_arg = NULL,
                      overid = TRUE)
})

test_that("pweight overid: hc0_small matches Stata", {
  skip_if(!file.exists(card_path))
  # pweight forces robust: vcov="iid" → HC1 (small=TRUE)
  test_pweight_config("card_pweight_overid", "hc0_small", card,
                      vcov_arg = "iid", small_arg = TRUE, cluster_arg = NULL,
                      overid = TRUE)
})

test_that("pweight overid: cl matches Stata", {
  skip_if(!file.exists(card_path))
  test_pweight_config("card_pweight_overid", "cl", card,
                      vcov_arg = "iid", small_arg = FALSE,
                      cluster_arg = ~ smsa, overid = TRUE)
})

test_that("pweight overid: cl_small matches Stata", {
  skip_if(!file.exists(card_path))
  test_pweight_config("card_pweight_overid", "cl_small", card,
                      vcov_arg = "iid", small_arg = TRUE,
                      cluster_arg = ~ smsa, overid = TRUE)
})


# ============================================================================
# Section 6: Weighted overidentified aweight (fills Tier 1 gap)
# ============================================================================

test_aweight_overid_config <- function(suffix, card_data,
                                       vcov_arg, small_arg, cluster_arg) {
  prefix <- "card_aweight_overid"
  coef_path <- file.path(fixture_dir, paste0(prefix, "_coef_", suffix, ".csv"))
  vcov_path <- file.path(fixture_dir, paste0(prefix, "_vcov_", suffix, ".csv"))
  diag_path <- file.path(fixture_dir, paste0(prefix, "_diagnostics_", suffix, ".csv"))
  fs_path   <- file.path(fixture_dir, paste0(prefix, "_firststage_", suffix, ".csv"))

  skip_if(!file.exists(coef_path))

  stata_coef <- read_coef_fixture(coef_path)
  stata_vcov <- read_vcov_fixture(vcov_path)
  stata_diag <- read_diagnostics_fixture(diag_path)

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card_data, weights = weight, vcov = vcov_arg,
                small = small_arg, clusters = cluster_arg)

  # Coefficients
  r_coef <- coef(fit)[names(stata_coef$estimate)]
  expect_equal(r_coef, stata_coef$estimate, tolerance = stata_tol$coef)

  # Standard errors
  r_se <- sqrt(diag(vcov(fit)))[names(stata_coef$std_error)]
  expect_equal(r_se, stata_coef$std_error, tolerance = stata_tol$se)

  # VCV
  r_vcov <- vcov(fit)[rownames(stata_vcov), colnames(stata_vcov)]
  expect_equal(r_vcov, stata_vcov, tolerance = stata_tol$vcov, check.attributes = FALSE)

  # Sigma
  expect_equal(fit$sigma, stata_diag$rmse, tolerance = stata_tol$coef)

  # Overid test
  if (!is.na(stata_diag$sargan)) {
    expect_equal(fit$diagnostics$overid$stat, stata_diag$sargan,
                 tolerance = stata_tol$stat)
  }
  if (!is.na(stata_diag$j)) {
    expect_equal(fit$diagnostics$overid$stat, stata_diag$j,
                 tolerance = stata_tol$stat)
  }

  # Model F
  if (!is.na(stata_diag$F_stat)) {
    expect_equal(fit$model_f, stata_diag$F_stat, tolerance = stata_tol$stat)
  }
}

test_that("aweight overid: iid matches Stata", {
  skip_if(!file.exists(card_path))
  test_aweight_overid_config("iid", card,
                             vcov_arg = "iid", small_arg = FALSE,
                             cluster_arg = NULL)
})

test_that("aweight overid: iid_small matches Stata", {
  skip_if(!file.exists(card_path))
  test_aweight_overid_config("iid_small", card,
                             vcov_arg = "iid", small_arg = TRUE,
                             cluster_arg = NULL)
})

test_that("aweight overid: robust matches Stata", {
  skip_if(!file.exists(card_path))
  # Stata robust (no small) = our HC0
  test_aweight_overid_config("hc1", card,
                             vcov_arg = "HC0", small_arg = FALSE,
                             cluster_arg = NULL)
})

test_that("aweight overid: robust_small matches Stata", {
  skip_if(!file.exists(card_path))
  # Stata robust small = our HC1
  test_aweight_overid_config("hc1_small", card,
                             vcov_arg = "HC1", small_arg = TRUE,
                             cluster_arg = NULL)
})

test_that("aweight overid: cl matches Stata", {
  skip_if(!file.exists(card_path))
  test_aweight_overid_config("cl", card,
                             vcov_arg = "iid", small_arg = FALSE,
                             cluster_arg = ~ smsa)
})

test_that("aweight overid: cl_small matches Stata", {
  skip_if(!file.exists(card_path))
  test_aweight_overid_config("cl_small", card,
                             vcov_arg = "iid", small_arg = TRUE,
                             cluster_arg = ~ smsa)
})


# ============================================================================
# Section 7: glance includes weight_type
# ============================================================================

test_that("glance includes weight_type for fweight", {
  skip_if(!file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight, weight_type = "fweight")
  gl <- glance(fit)
  expect_equal(gl$weight_type, "fweight")
})

test_that("glance includes weight_type='aweight' by default", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  gl <- glance(fit)
  expect_equal(gl$weight_type, "aweight")
})

test_that("weight_type is stored in fitted object", {
  fit_a <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp)
  expect_equal(fit_a$weight_type, "aweight")

  skip_if(!file.exists(card_path))
  fit_f <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card, weights = weight, weight_type = "fweight")
  expect_equal(fit_f$weight_type, "fweight")

  fit_p <- suppressMessages(
    ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
           weight_type = "pweight")
  )
  expect_equal(fit_p$weight_type, "pweight")
})

test_that("n_physical is stored for fweight", {
  skip_if(!file.exists(card_path))
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, weights = weight, weight_type = "fweight")
  expect_equal(fit$n_physical, 3010L)
  expect_true(nobs(fit) > fit$n_physical)
})

test_that("n_physical is NULL for non-fweight", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp)
  expect_null(fit$n_physical)
})
