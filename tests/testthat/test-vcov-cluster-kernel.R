# ============================================================================
# Tests: Cluster+Kernel VCE — Kiefer, Driscoll-Kraay, Thompson (Ticket M4)
# ============================================================================

# --- Helpers and data ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)

wp_data_path <- file.path(fixture_dir, "wp_ck_data.csv")
if (file.exists(wp_data_path)) {
  wp <- read.csv(wp_data_path)
}

# Helper: read VCV fixture (no term column — just v1, v2, ...)
read_ck_vcov <- function(path) {
  as.matrix(read.csv(path))
}

# Helper: compare VCV matrices element-wise
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

# Helper: load fixture and compare a cluster+kernel fit
check_ck_fixture <- function(fit, prefix, suffix, check_vcov = TRUE,
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

  # Match VCV
  if (check_vcov && file.exists(vcov_path)) {
    V_stata <- read_ck_vcov(vcov_path)
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

    # Overid (when available)
    if (!is.na(stata_diag$overid_stat) && !is.null(diag$overid) &&
        diag$overid$df > 0L) {
      expect_equal(diag$overid$stat, stata_diag$overid_stat,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "overid stat"))
      expect_equal(diag$overid$p, stata_diag$overid_p,
                   tolerance = stata_tol$pval,
                   info = paste(suffix, "overid p"))
    }

    # Underid (when available)
    if (!is.na(stata_diag$underid_stat) && !is.null(diag$underid)) {
      expect_equal(diag$underid$stat, stata_diag$underid_stat,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "underid stat"))
      expect_equal(diag$underid$p, stata_diag$underid_p,
                   tolerance = stata_tol$pval,
                   info = paste(suffix, "underid p"))
    }

    # Weak ID: Cragg-Donald F
    if (!is.na(stata_diag$weak_id_cd_f) && !is.null(diag$weak_id)) {
      expect_equal(diag$weak_id$stat, stata_diag$weak_id_cd_f,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "CD F"))
    }

    # Weak ID: KP F
    if (!is.na(stata_diag$weak_id_kp_f) && !is.null(diag$weak_id_robust)) {
      expect_equal(diag$weak_id_robust$stat, stata_diag$weak_id_kp_f,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "KP F"))
    }

    # Model F
    if (!is.na(stata_diag$model_f) && !is.null(fit$model_f)) {
      expect_equal(fit$model_f, stata_diag$model_f,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "model F"))
    }

    # AR test (when available)
    if (!is.na(stata_diag$ar_f) && !is.null(diag$anderson_rubin)) {
      expect_equal(diag$anderson_rubin$f_stat, stata_diag$ar_f,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "AR F"))
    }

    # Endogeneity test (when available)
    if (!is.na(stata_diag$endog_stat) && !is.null(diag$endogeneity) &&
        !is.na(diag$endogeneity$stat)) {
      expect_equal(diag$endogeneity$stat, stata_diag$endog_stat,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "endogeneity stat"))
    }
  }
}


# ============================================================================
# Guard tests
# ============================================================================

test_that("kiefer requires panel data (tvar + ivar)", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, kiefer = TRUE),
    "kiefer requires panel data"
  )
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, kiefer = TRUE, tvar = "wt"),
    "kiefer requires panel data"
  )
})

test_that("kiefer is incompatible with robust VCE or clustering", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | hours | educ, data = wp,
           kiefer = TRUE, vcov = "HC1", tvar = "year", ivar = "nr"),
    "kiefer is incompatible with robust VCE or clustering"
  )
  expect_error(
    ivreg2(lwage ~ exper + expersq | hours | educ, data = wp,
           kiefer = TRUE, clusters = ~nr, tvar = "year", ivar = "nr"),
    "kiefer is incompatible with robust VCE or clustering"
  )
})

test_that("kiefer is incompatible with non-Truncated kernel", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | hours | educ, data = wp,
           kiefer = TRUE, kernel = "Bartlett", tvar = "year", ivar = "nr"),
    "kiefer is incompatible with kernel='Bartlett'"
  )
})

test_that("kiefer is incompatible with explicit bw", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | hours | educ, data = wp,
           kiefer = TRUE, bw = 3, tvar = "year", ivar = "nr"),
    "kiefer is incompatible with explicit `bw`"
  )
})

test_that("dkraay must be positive numeric scalar", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | hours | educ, data = wp,
           dkraay = -1, tvar = "year", ivar = "nr"),
    "positive numeric scalar"
  )
  expect_error(
    ivreg2(lwage ~ exper + expersq | hours | educ, data = wp,
           dkraay = "foo", tvar = "year", ivar = "nr"),
    "positive numeric scalar"
  )
})

test_that("dkraay requires panel data (tvar + ivar)", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, dkraay = 3),
    "dkraay requires panel data"
  )
})

test_that("cannot specify both dkraay and bw", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | hours | educ, data = wp,
           dkraay = 3, bw = 5, tvar = "year", ivar = "nr"),
    "cannot specify both `dkraay` and `bw`"
  )
})

test_that("cluster+kernel requires cluster vars to match ivar and tvar", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  wp$other <- rep(1:10, length.out = nrow(wp))
  expect_error(
    ivreg2(lwage ~ exper + expersq | hours | educ, data = wp,
           clusters = ~nr + other, kernel = "bartlett", bw = 3,
           tvar = "year", ivar = "nr"),
    "cluster variables to match ivar and tvar"
  )
})

test_that("one-way cluster+kernel requires clustering on tvar", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  expect_error(
    ivreg2(lwage ~ exper + expersq | hours | educ, data = wp,
           clusters = ~nr, kernel = "bartlett", bw = 3,
           tvar = "year", ivar = "nr"),
    "clustering on the time variable"
  )
})


# ============================================================================
# Kiefer equivalence: kiefer=TRUE == vcov="AC", kernel="Truncated", bw=T
# ============================================================================

test_that("kiefer=TRUE is equivalent to AC+Truncated+fullBW", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")

  fit_kiefer <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, kiefer = TRUE, tvar = "year", ivar = "nr"
  )

  # Kiefer internally uses bw = T_span = 8, but the max_bw check for
  # non-kiefer AC is (T_span-1)/tdelta = 7. The VCVs are identical
  # because lag = T_span has zero observation pairs. Use bw = 7 for
  # the manual fit to stay within the allowed range.
  # Note: diagnostics differ (kiefer uses IID path, manual AC uses KP
  # path), so we only compare VCV here.
  years <- sort(unique(wp$year))
  T_span <- max(years) - min(years) + 1L
  max_bw <- (T_span - 1L) / 1L  # tdelta = 1

  fit_manual <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, vcov = "AC", kernel = "truncated", bw = max_bw,
    tvar = "year", ivar = "nr"
  )

  # Same VCV to machine precision
  expect_equal(fit_kiefer$vcov, fit_manual$vcov,
               tolerance = 1e-12)
  expect_equal(coef(fit_kiefer), coef(fit_manual))
})


# ============================================================================
# DK equivalence: dkraay=3 == clusters=~year, kernel="bar", bw=3
# ============================================================================

test_that("dkraay=3 is equivalent to explicit cluster+kernel on tvar", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")

  fit_dk <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, dkraay = 3, tvar = "year", ivar = "nr"
  )

  fit_manual <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, clusters = ~year, kernel = "bartlett", bw = 3,
    tvar = "year", ivar = "nr"
  )

  # Same VCV to machine precision
  expect_equal(fit_dk$vcov, fit_manual$vcov,
               tolerance = 1e-12)
  expect_equal(coef(fit_dk), coef(fit_manual))
})


# ============================================================================
# Kiefer Stata parity
# ============================================================================

test_that("Kiefer overid matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, kiefer = TRUE, tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_kiefer", "overid")
})

test_that("Kiefer overid small matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, kiefer = TRUE, small = TRUE, tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_kiefer", "overid_small")
})

test_that("Kiefer justid matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ,
    data = wp, kiefer = TRUE, tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_kiefer", "justid")
})


# ============================================================================
# Driscoll-Kraay Stata parity
# ============================================================================

test_that("DK Bartlett bw=3 overid matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, dkraay = 3, tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_dk", "bar3_overid")
})

test_that("DK Bartlett bw=3 overid small matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, dkraay = 3, small = TRUE, tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_dk", "bar3_overid_small")
})

test_that("DK Parzen bw=4 overid matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, dkraay = 4, kernel = "parzen", tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_dk", "par4_overid")
})

test_that("DK Truncated bw=2 overid matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, dkraay = 2, kernel = "truncated", tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_dk", "tru2_overid")
})

test_that("DK Bartlett bw=3 justid matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ,
    data = wp, dkraay = 3, tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_dk", "bar3_justid")
})


# ============================================================================
# Thompson (two-way cluster + kernel) Stata parity
# ============================================================================

test_that("Thompson Bartlett bw=3 overid matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, clusters = ~nr + year, kernel = "bartlett", bw = 3,
    tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_ck", "bar3_overid")
})

test_that("Thompson Bartlett bw=3 overid small matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, clusters = ~nr + year, kernel = "bartlett", bw = 3,
    small = TRUE, tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_ck", "bar3_overid_small")
})

test_that("Thompson Parzen bw=4 overid matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, clusters = ~nr + year, kernel = "parzen", bw = 4,
    tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_ck", "par4_overid")
})

test_that("Thompson Bartlett bw=3 justid matches Stata", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ,
    data = wp, clusters = ~nr + year, kernel = "bartlett", bw = 3,
    tvar = "year", ivar = "nr"
  )
  check_ck_fixture(fit, "wp_ck", "bar3_justid")
})


# ============================================================================
# Structural checks
# ============================================================================

test_that("kiefer sets vcov_type to AC", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, kiefer = TRUE, tvar = "year", ivar = "nr"
  )
  expect_equal(fit$vcov_type, "AC")
  expect_true(fit$kiefer)
})

test_that("dkraay sets vcov_type to CL and stores dkraay", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, dkraay = 3, tvar = "year", ivar = "nr"
  )
  expect_equal(fit$vcov_type, "CL")
  expect_equal(fit$dkraay, 3)
  expect_equal(fit$kernel, "Bartlett")
  expect_equal(fit$bw, 3)
})

test_that("Thompson two-way sets vcov_type to CL", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, clusters = ~nr + year, kernel = "bartlett", bw = 3,
    tvar = "year", ivar = "nr"
  )
  expect_equal(fit$vcov_type, "CL")
  expect_equal(fit$kernel, "Bartlett")
  expect_equal(fit$bw, 3)
})

test_that("glance includes kiefer and dkraay columns", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, dkraay = 3, tvar = "year", ivar = "nr"
  )
  gl <- glance(fit)
  expect_true("kiefer" %in% names(gl))
  expect_true("dkraay" %in% names(gl))
  expect_equal(gl$dkraay, 3)
})

test_that("Thompson normalizes cluster order (ivar first)", {
  skip_if(!file.exists(wp_data_path), "wagepan data not found")

  # Specify in reverse order: year + nr instead of nr + year
  fit <- ivreg2(
    lwage ~ exper + expersq + married + union | hours | educ + black,
    data = wp, clusters = ~year + nr, kernel = "bartlett", bw = 3,
    tvar = "year", ivar = "nr"
  )

  # Should be normalized so cluster_var[1] = ivar
  expect_equal(fit$cluster_var[1L], "nr")
  expect_equal(fit$cluster_var[2L], "year")
})
