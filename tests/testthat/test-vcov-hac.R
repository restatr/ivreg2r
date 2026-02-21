# ============================================================================
# Tests: HAC and AC VCV estimation (Ticket M2)
# ============================================================================

# --- Helpers and data ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)

hac_data_path <- file.path(fixture_dir, "ts_hac_data.csv")
gap_data_path <- file.path(fixture_dir, "ts_gap_hac_data.csv")

if (file.exists(hac_data_path)) {
  ts_data <- read.csv(hac_data_path)
}
if (file.exists(gap_data_path)) {
  ts_gap_data <- read.csv(gap_data_path)
}

# Helper: read HAC VCV fixture (no term column — just v1, v2, v3)
read_hac_vcov <- function(path) {
  as.matrix(read.csv(path))
}

# Helper: compare VCV matrices element-wise (both unnamed numeric matrices)
expect_vcov_match <- function(V_r, V_stata, tol = stata_tol$vcov, label = "") {
  V_r <- unname(as.matrix(V_r))
  V_stata <- unname(as.matrix(V_stata))
  P <- nrow(V_stata)
  for (i in seq_len(P)) {
    for (j in seq_len(P)) {
      expect_equal(
        V_r[i, j], V_stata[i, j],
        tolerance = tol,
        info = paste0(label, " VCV[", i, ",", j, "]")
      )
    }
  }
}

# Helper: load fixture and compare a HAC/AC fit
check_hac_fixture <- function(fit, prefix, suffix, check_vcov = TRUE,
                              check_diag = TRUE) {
  coef_path <- file.path(fixture_dir, paste0(prefix, "_coef_", suffix, ".csv"))
  vcov_path <- file.path(fixture_dir, paste0(prefix, "_vcov_", suffix, ".csv"))
  diag_path <- file.path(fixture_dir,
                         paste0(prefix, "_diagnostics_", suffix, ".csv"))
  skip_if(!file.exists(coef_path), "Fixture not found")

  stata_coef <- read.csv(coef_path, stringsAsFactors = FALSE)
  stata_coef$term[stata_coef$term == "_cons"] <- "(Intercept)"

  # R variable names and ordering
  r_names <- names(coef(fit))

  # Match coefficients by name
  m <- match(r_names, stata_coef$term)
  expect_equal(unname(coef(fit)), stata_coef$estimate[m],
               tolerance = stata_tol$coef,
               info = paste(suffix, "coefficients"))
  expect_equal(unname(sqrt(diag(fit$vcov))), stata_coef$std_error[m],
               tolerance = stata_tol$se,
               info = paste(suffix, "std errors"))

  # Match VCV (reorder R to match Stata's variable ordering)
  if (check_vcov && file.exists(vcov_path)) {
    V_stata <- read_hac_vcov(vcov_path)
    # Stata ordering (from coef file): map to R positions
    stata_to_r <- match(stata_coef$term, r_names)
    V_r_reordered <- unname(fit$vcov[stata_to_r, stata_to_r])
    expect_vcov_match(V_r_reordered, V_stata, label = suffix)
  }

  # Match diagnostics
  if (check_diag && file.exists(diag_path)) {
    stata_diag <- read.csv(diag_path)
    diag <- fit$diagnostics

    # N
    expect_identical(fit$nobs, as.integer(stata_diag$N))

    # Overidentification test
    if (!is.na(stata_diag$overid_stat) && stata_diag$overid_stat != 0) {
      expect_equal(diag$overid$stat, stata_diag$overid_stat,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "overid stat"))
      expect_equal(diag$overid$p, stata_diag$overid_p,
                   tolerance = stata_tol$pval,
                   info = paste(suffix, "overid p"))
      expect_identical(diag$overid$df, as.integer(stata_diag$overid_df),
                       info = paste(suffix, "overid df"))
    }

    # Underidentification test (KP rk LM)
    if (!is.na(stata_diag$underid_stat)) {
      expect_equal(diag$underid$stat, stata_diag$underid_stat,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "underid stat"))
      expect_equal(diag$underid$p, stata_diag$underid_p,
                   tolerance = stata_tol$pval,
                   info = paste(suffix, "underid p"))
    }

    # Weak identification (CD F)
    if (!is.na(stata_diag$weak_id_cd_f)) {
      expect_equal(diag$weak_id$stat, stata_diag$weak_id_cd_f,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "CD F"))
    }

    # Weak identification (KP Wald F)
    if (!is.na(stata_diag$weak_id_kp_f)) {
      expect_equal(diag$weak_id_robust$stat, stata_diag$weak_id_kp_f,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "KP Wald F"))
    }
  }
}


# ============================================================================
# Unit tests: .build_time_index()
# ============================================================================

test_that(".build_time_index works for simple time series", {
  tvar <- c(1, 2, 3, 4, 5)
  ti <- ivreg2r:::.build_time_index(tvar)
  expect_equal(ti$T_span, 5L)
  expect_equal(ti$tdelta, 1)
  expect_equal(ti$n_gaps, 0L)
  expect_equal(ti$tvar_sorted, 1:5)
  expect_null(ti$panel_info)
})

test_that(".build_time_index detects gaps", {
  tvar <- c(1, 2, 4, 5)  # gap at t=3
  ti <- ivreg2r:::.build_time_index(tvar)
  expect_equal(ti$T_span, 5L)
  expect_equal(ti$tdelta, 1)
  expect_equal(ti$n_gaps, 1L)
})

test_that(".build_time_index sorts unsorted data", {
  tvar <- c(5, 3, 1, 4, 2)
  ti <- ivreg2r:::.build_time_index(tvar)
  expect_equal(ti$tvar_sorted, 1:5)
  expect_equal(ti$sort_order, c(3, 5, 2, 4, 1))
})

test_that(".build_time_index handles non-unit spacing", {
  tvar <- c(2, 4, 6, 8, 10)
  ti <- ivreg2r:::.build_time_index(tvar)
  expect_equal(ti$T_span, 9L)
  expect_equal(ti$tdelta, 2)
  expect_equal(ti$n_gaps, 0L)
})

test_that(".build_time_index works with panel data", {
  tvar <- c(1, 2, 3, 1, 2, 3)
  ivar <- c("a", "a", "a", "b", "b", "b")
  ti <- ivreg2r:::.build_time_index(tvar, ivar)
  expect_equal(ti$T_span, 3L)
  expect_equal(ti$tdelta, 1)
  expect_equal(ti$n_gaps, 0L)
  expect_length(ti$panel_info, 2L)
})

test_that(".build_time_index errors on duplicate time values", {
  expect_error(
    ivreg2r:::.build_time_index(c(1, 2, 2, 3)),
    "duplicate"
  )
})

test_that(".build_time_index errors on duplicate panel-time pairs", {
  expect_error(
    ivreg2r:::.build_time_index(c(1, 2, 2, 1, 2, 3),
                                c("a", "a", "a", "b", "b", "b")),
    "Duplicate"
  )
})

test_that(".build_time_index allows same time across different panels", {
  tvar <- c(1, 2, 1, 2)
  ivar <- c("a", "a", "b", "b")
  ti <- ivreg2r:::.build_time_index(tvar, ivar)
  expect_equal(ti$T_span, 2L)
})


# ============================================================================
# Unit tests: .lag_pairs()
# ============================================================================

test_that(".lag_pairs returns correct pairs for simple series", {
  tvar <- 1:5
  ti <- ivreg2r:::.build_time_index(tvar)
  pairs <- ivreg2r:::.lag_pairs(ti, 1)
  # At tau=1: obs 2 pairs with 1, 3 pairs with 2, etc.
  expect_equal(nrow(pairs), 4L)
  expect_equal(pairs[, "i_now"], 2:5)
  expect_equal(pairs[, "i_lag"], 1:4)
})

test_that(".lag_pairs handles gaps correctly", {
  tvar <- c(1, 2, 4, 5)  # gap at t=3
  ti <- ivreg2r:::.build_time_index(tvar)
  pairs <- ivreg2r:::.lag_pairs(ti, 1)
  # Tau=1: (2,1), (5,4) = indices (2,1), (4,3) in sorted data
  expect_equal(nrow(pairs), 2L)
})

test_that(".lag_pairs returns zero rows for too-large lag", {
  tvar <- 1:5
  ti <- ivreg2r:::.build_time_index(tvar)
  pairs <- ivreg2r:::.lag_pairs(ti, 10)
  expect_equal(nrow(pairs), 0L)
})


# ============================================================================
# Identity check: HAC Bartlett bw=1 should equal HC0
# ============================================================================

test_that("HAC Bartlett bw=1 equals HC0 (zero off-diagonal lags)", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit_hac <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                    vcov = "HAC", kernel = "bartlett", bw = 1, tvar = "t")
  fit_hc0 <- ivreg2(y ~ w | x | z1 + z2, data = ts_data, vcov = "HC0")
  expect_equal(fit_hac$vcov, fit_hc0$vcov, tolerance = 1e-12)
})


# ============================================================================
# Input validation
# ============================================================================

test_that("HAC/AC requires tvar", {
  expect_error(
    ivreg2(y ~ w | x | z1 + z2, data = ts_data,
           vcov = "HAC", kernel = "bartlett", bw = 3),
    "tvar"
  )
})

test_that("kernel without bw is an error", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  expect_error(
    ivreg2(y ~ w | x | z1 + z2, data = ts_data,
           kernel = "bartlett", tvar = "t"),
    "bw"
  )
})

test_that("fweight + kernel is an error", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  ts_data$fw <- rep(1L, nrow(ts_data))
  expect_error(
    ivreg2(y ~ w | x | z1 + z2, data = ts_data,
           vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t",
           weights = fw, weight_type = "fweight"),
    "fweight"
  )
})

test_that("bandwidth span check: bw too large is an error", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  expect_error(
    ivreg2(y ~ w | x | z1 + z2, data = ts_data,
           vcov = "HAC", kernel = "bartlett", bw = 500, tvar = "t"),
    "bandwidth"
  )
})

test_that("kernel + clusters is an error (HAC-CL not implemented)", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  ts_data$cl <- rep(1:10, length.out = nrow(ts_data))
  expect_error(
    ivreg2(y ~ w | x | z1 + z2, data = ts_data,
           vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t",
           clusters = ~cl),
    "HAC-CL"
  )
})

test_that("duplicate time values are an error", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  dup_data <- ts_data
  dup_data$t[2] <- dup_data$t[1]  # introduce duplicate
  expect_error(
    ivreg2(y ~ w | x | z1 + z2, data = dup_data,
           vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t"),
    "duplicate"
  )
})

test_that("HAC fit y slot is in original row order", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t")
  # y should equal fitted + resid (both in original order)
  expect_equal(unname(fit$y), unname(fit$fitted.values + fit$residuals),
               tolerance = 1e-14)
  # y should match the original data column
  expect_equal(unname(fit$y), ts_data$y)
})

test_that("HAC fit x slot is in original row order", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t",
                x = TRUE)
  # X and Z rows should align with model frame / original data
  expect_equal(unname(fit$x$X[, "w"]), ts_data$w)
  expect_equal(unname(fit$x$Z[, "z1"]), ts_data$z1)
  expect_equal(unname(fit$x$Z[, "z2"]), ts_data$z2)
})


# ============================================================================
# VCE inference: kernel implies VCE type
# ============================================================================

test_that("kernel + iid infers AC", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                kernel = "bartlett", bw = 3, tvar = "t")
  expect_equal(fit$vcov_type, "AC")
})

test_that("kernel + robust infers HAC", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HC1", kernel = "bartlett", bw = 3, tvar = "t")
  expect_equal(fit$vcov_type, "HAC")
})


# ============================================================================
# Object metadata
# ============================================================================

test_that("HAC fit stores kernel, bw, tvar metadata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t")
  expect_equal(fit$kernel, "Bartlett")
  expect_equal(fit$bw, 3)
  expect_equal(fit$tvar, "t")
  expect_null(fit$ivar)
})

test_that("glance includes kernel and bw", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t")
  gl <- glance(fit)
  expect_equal(gl$kernel, "Bartlett")
  expect_equal(gl$bw, 3)
})

test_that("glance kernel/bw are NA for non-HAC fit", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, small = TRUE)
  gl <- glance(fit)
  expect_true(is.na(gl$kernel))
  expect_true(is.na(gl$bw))
})


# ============================================================================
# HAC Stata fixtures: overidentified model (y w (x = z1 z2))
# ============================================================================

test_that("HAC Bartlett bw=3 matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t")
  check_hac_fixture(fit, "ts_hac", "bartlett_bw3")
})

test_that("HAC Bartlett bw=5 matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 5, tvar = "t")
  check_hac_fixture(fit, "ts_hac", "bartlett_bw5")
})

test_that("HAC Parzen bw=4 matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "parzen", bw = 4, tvar = "t")
  check_hac_fixture(fit, "ts_hac", "parzen_bw4")
})

test_that("HAC Truncated bw=3 matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "truncated", bw = 3, tvar = "t")
  check_hac_fixture(fit, "ts_hac", "truncated_bw3")
})

test_that("HAC Tukey-Hanning bw=3 matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "tukey-hanning", bw = 3, tvar = "t")
  check_hac_fixture(fit, "ts_hac", "thanning_bw3")
})

test_that("HAC Tukey-Hamming bw=3 matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "tukey-hamming", bw = 3, tvar = "t")
  check_hac_fixture(fit, "ts_hac", "thamming_bw3")
})

test_that("HAC Daniell bw=4 matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "daniell", bw = 4, tvar = "t")
  check_hac_fixture(fit, "ts_hac", "daniell_bw4")
})

test_that("HAC QS bw=3 just-identified matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1, data = ts_data,
                vcov = "HAC", kernel = "quadratic spectral", bw = 3, tvar = "t")
  check_hac_fixture(fit, "ts_hac", "qs_bw3_justid")
})

test_that("HAC QS bw=3 overidentified matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "quadratic spectral", bw = 3, tvar = "t")
  check_hac_fixture(fit, "ts_hac", "qs_bw3")
})


# ============================================================================
# AC Stata fixtures: overidentified model
# ============================================================================

test_that("AC Bartlett bw=3 matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "AC", kernel = "bartlett", bw = 3, tvar = "t")
  check_hac_fixture(fit, "ts_ac", "bartlett_bw3")
})

test_that("AC Bartlett bw=5 just-identified matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1, data = ts_data,
                vcov = "AC", kernel = "bartlett", bw = 5, tvar = "t")
  check_hac_fixture(fit, "ts_ac", "bartlett_bw5_justid")
})

test_that("AC Parzen bw=4 matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "AC", kernel = "parzen", bw = 4, tvar = "t")
  check_hac_fixture(fit, "ts_ac", "parzen_bw4")
})

test_that("AC uses Sargan (not Hansen J) for overid test", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "AC", kernel = "bartlett", bw = 3, tvar = "t")
  expect_equal(fit$diagnostics$overid$test_name, "Sargan")
})

test_that("HAC uses Hansen J for overid test", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t")
  expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
})


# ============================================================================
# Gappy time-series HAC fixture
# ============================================================================

test_that("HAC Bartlett bw=3 with gappy data matches Stata", {
  skip_if(!file.exists(gap_data_path), "Gap HAC data not found")
  fit <- suppressWarnings(
    ivreg2(y ~ w | x | z1 + z2, data = ts_gap_data,
           vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t")
  )
  check_hac_fixture(fit, "ts_gap_hac", "bartlett_bw3")
})

test_that("Gappy data produces a warning about gaps", {
  skip_if(!file.exists(gap_data_path), "Gap HAC data not found")
  expect_warning(
    ivreg2(y ~ w | x | z1 + z2, data = ts_gap_data,
           vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t"),
    "gap"
  )
})


# ============================================================================
# Symmetry and PSD checks
# ============================================================================

test_that("HAC VCV matrix is symmetric", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t")
  expect_equal(fit$vcov, t(fit$vcov), tolerance = 1e-15)
})

test_that("AC VCV matrix is symmetric", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "AC", kernel = "bartlett", bw = 3, tvar = "t")
  expect_equal(fit$vcov, t(fit$vcov), tolerance = 1e-15)
})

test_that("HAC VCV diagonal is positive", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t")
  expect_true(all(diag(fit$vcov) > 0))
})


# ============================================================================
# Auto-bandwidth: unit tests for .auto_bandwidth()
# ============================================================================

test_that(".auto_bandwidth returns positive numeric scalar", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t")
  ti <- ivreg2r:::.build_time_index(ts_data$t)
  bw <- ivreg2r:::.auto_bandwidth(
    resid = fit$residuals[ti$sort_order],
    Z = model.matrix(~ w + z1 + z2, data = ts_data)[ti$sort_order, ],
    time_index = ti, kernel = "Bartlett", has_intercept = TRUE,
    N = nrow(ts_data)
  )
  expect_true(is.numeric(bw))
  expect_length(bw, 1L)
  expect_true(bw >= 1)
})

test_that(".auto_bandwidth Bartlett/Parzen return integer bw", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t")
  ti <- ivreg2r:::.build_time_index(ts_data$t)
  Z <- model.matrix(~ w + z1 + z2, data = ts_data)[ti$sort_order, ]
  resid_sorted <- fit$residuals[ti$sort_order]

  bw_bart <- ivreg2r:::.auto_bandwidth(resid_sorted, Z, ti,
                                        "Bartlett", TRUE, nrow(ts_data))
  bw_parz <- ivreg2r:::.auto_bandwidth(resid_sorted, Z, ti,
                                        "Parzen", TRUE, nrow(ts_data))
  expect_equal(bw_bart, trunc(bw_bart))
  expect_equal(bw_parz, trunc(bw_parz))
})

test_that(".auto_bandwidth QS may return non-integer bw", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = 3, tvar = "t")
  ti <- ivreg2r:::.build_time_index(ts_data$t)
  Z <- model.matrix(~ w + z1 + z2, data = ts_data)[ti$sort_order, ]
  resid_sorted <- fit$residuals[ti$sort_order]

  bw_qs <- ivreg2r:::.auto_bandwidth(resid_sorted, Z, ti,
                                      "Quadratic Spectral", TRUE,
                                      nrow(ts_data))
  # QS can return fractional; just check it's positive
  expect_true(bw_qs > 0)
})

test_that(".auto_bandwidth warns and returns 1 for mstar=0", {
  # mstar = trunc(20*(tobs/100)^(2/9)) = 0 requires tobs < 0.00014.
  # Construct a contrived time_index with very small T_span/tdelta ratio.
  ti <- list(
    sort_order = 1L, tvar_sorted = 1, ivar_sorted = NULL,
    T_span = 1L, tdelta = 10000, n_gaps = 0L, panel_info = NULL
  )
  # tobs = T_span / tdelta = 1/10000 = 0.0001 → mstar = 0
  expect_warning(
    bw <- ivreg2r:::.auto_bandwidth(1, matrix(1, 1, 1), ti,
                                     "Bartlett", FALSE, 1L),
    "Time span too short"
  )
  expect_equal(bw, 1)
})

test_that(".auto_bandwidth warns and returns 1 for zero-variance scores", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  ti <- ivreg2r:::.build_time_index(1:20)
  # All-zero residuals → zero scores → shat_0 ≈ 0
  resid0 <- rep(0, 20)
  Z <- matrix(rnorm(40), 20, 2)
  expect_warning(
    bw <- ivreg2r:::.auto_bandwidth(resid0, Z, ti,
                                     "Bartlett", FALSE, 20L),
    "No autocorrelation"
  )
  expect_equal(bw, 1)
})

test_that("bw = 'auto' with panel data (ivar) gives error", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  # Create a minimal panel dataset
  panel_data <- rbind(
    data.frame(id = 1, t = 1:50, y = rnorm(50), w = rnorm(50),
               x = rnorm(50), z1 = rnorm(50), z2 = rnorm(50)),
    data.frame(id = 2, t = 1:50, y = rnorm(50), w = rnorm(50),
               x = rnorm(50), z1 = rnorm(50), z2 = rnorm(50))
  )
  expect_error(
    ivreg2(y ~ w | x | z1 + z2, data = panel_data,
           kernel = "bartlett", bw = "auto", tvar = "t", ivar = "id"),
    "panel data"
  )
})

test_that("bw = 'auto' with unsupported kernel gives error", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  expect_error(
    ivreg2(y ~ w | x | z1 + z2, data = ts_data,
           kernel = "truncated", bw = "auto", tvar = "t"),
    "not compatible"
  )
})

test_that("bw = 'auto' caps to max_bw with warning for short series", {
  set.seed(999)
  n <- 10  # Very short series: max_bw = (10-1)/1 = 9
  short_data <- data.frame(
    t = seq_len(n),
    y = rnorm(n),
    w = rnorm(n),
    x = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n)
  )
  expect_warning(
    ivreg2(y ~ w | x | z1 + z2, data = short_data,
           vcov = "HAC", kernel = "bartlett", bw = "auto", tvar = "t"),
    "exceeds time-span limit"
  )
  # Re-run with suppressWarnings to inspect the result
  fit <- suppressWarnings(
    ivreg2(y ~ w | x | z1 + z2, data = short_data,
           vcov = "HAC", kernel = "bartlett", bw = "auto", tvar = "t")
  )
  expect_true(fit$bw <= 9)
})


# ============================================================================
# Auto-bandwidth: Stata fixture tests
# ============================================================================

# Helper: read auto-bw fixture bandwidth value
read_autobiw <- function(prefix, suffix) {
  path <- file.path(fixture_dir, paste0(prefix, "_bw_", suffix, ".csv"))
  if (!file.exists(path)) return(NA_real_)
  as.numeric(read.csv(path)$bw)
}

test_that("HAC auto Bartlett bandwidth matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  stata_bw <- read_autobiw("ts_hac", "auto_bartlett")
  skip_if(is.na(stata_bw), "Auto-bw fixture not found")

  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = "auto", tvar = "t")
  expect_equal(fit$bw, stata_bw, info = "auto Bartlett bw")
  check_hac_fixture(fit, "ts_hac", "auto_bartlett")
})

test_that("HAC auto Parzen bandwidth matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  stata_bw <- read_autobiw("ts_hac", "auto_parzen")
  skip_if(is.na(stata_bw), "Auto-bw fixture not found")

  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "parzen", bw = "auto", tvar = "t")
  expect_equal(fit$bw, stata_bw, info = "auto Parzen bw")
  check_hac_fixture(fit, "ts_hac", "auto_parzen")
})

test_that("HAC auto QS bandwidth matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  stata_bw <- read_autobiw("ts_hac", "auto_qs")
  skip_if(is.na(stata_bw), "Auto-bw fixture not found")

  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "qs", bw = "auto", tvar = "t")
  # QS bandwidth is fractional; compare within stat tolerance
  # (Stata truncates via strofreal → ~7 sig digits, so ~1e-5 bw noise)
  expect_equal(fit$bw, stata_bw, tolerance = stata_tol$stat,
               info = "auto QS bw")
  # Use stat tolerance for VCV/SEs (fractional bw noise cascades through
  # spectral window → ~1e-6 VCV differences, at the vcov tolerance boundary)
  coef_path <- file.path(fixture_dir, "ts_hac_coef_auto_qs.csv")
  stata_coef <- read.csv(coef_path, stringsAsFactors = FALSE)
  stata_coef$term[stata_coef$term == "_cons"] <- "(Intercept)"
  r_names <- names(coef(fit))
  m <- match(r_names, stata_coef$term)
  expect_equal(unname(coef(fit)), stata_coef$estimate[m],
               tolerance = stata_tol$coef, info = "auto QS coefficients")
  expect_equal(unname(sqrt(diag(fit$vcov))), stata_coef$std_error[m],
               tolerance = stata_tol$stat, info = "auto QS std errors")
})

test_that("AC auto Bartlett bandwidth matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  stata_bw <- read_autobiw("ts_ac", "auto_bartlett")
  skip_if(is.na(stata_bw), "Auto-bw fixture not found")

  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                kernel = "bartlett", bw = "auto", tvar = "t")
  expect_equal(fit$bw, stata_bw, info = "auto AC Bartlett bw")
  check_hac_fixture(fit, "ts_ac", "auto_bartlett")
})

test_that("HAC auto Bartlett just-identified bandwidth matches Stata", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  stata_bw <- read_autobiw("ts_hac", "auto_bartlett_justid")
  skip_if(is.na(stata_bw), "Auto-bw fixture not found")

  fit <- ivreg2(y ~ w | x | z1, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = "auto", tvar = "t")
  expect_equal(fit$bw, stata_bw, info = "auto Bartlett justid bw")
  check_hac_fixture(fit, "ts_hac", "auto_bartlett_justid")
})

test_that("HAC auto Bartlett gappy data bandwidth matches Stata", {
  skip_if(!file.exists(gap_data_path), "Gap HAC data not found")
  stata_bw <- read_autobiw("ts_gap_hac", "auto_bartlett")
  skip_if(is.na(stata_bw), "Auto-bw fixture not found")

  fit <- suppressWarnings(
    ivreg2(y ~ w | x | z1 + z2, data = ts_gap_data,
           vcov = "HAC", kernel = "bartlett", bw = "auto", tvar = "t")
  )
  expect_equal(fit$bw, stata_bw, info = "auto Bartlett gappy bw")
  check_hac_fixture(fit, "ts_gap_hac", "auto_bartlett")
})

test_that("bw = 'auto' stores resolved numeric in fit$bw", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = "auto", tvar = "t")
  expect_true(is.numeric(fit$bw))
  expect_true(fit$bw >= 1)
})

test_that("glance reports resolved bw for auto-bandwidth", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = "auto", tvar = "t")
  gl <- glance(fit)
  expect_true(is.numeric(gl$bw))
  expect_true(gl$bw >= 1)
})
