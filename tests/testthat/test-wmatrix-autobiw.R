# ============================================================================
# Tests: wmatrix + bw="auto" — auto-bandwidth uses W-step residuals
# ============================================================================
# Verifies that when wmatrix is supplied with bw="auto", the bandwidth is
# selected using W-step residuals (matching Stata ivreg2.ado:970-980), not
# 2SLS residuals.

# --- Helpers and data ---
hac_data_path <- fixture_path("ts_hac_data.csv")
if (file.exists(hac_data_path)) {
  ts_data <- read.csv(hac_data_path)
}

read_autobiw <- function(prefix, suffix) {
  path <- fixture_path(paste0(prefix, "_bw_", suffix, ".csv"))
  if (!file.exists(path)) return(NA_real_)
  as.numeric(read.csv(path)$bw)
}

check_fixture <- function(fit, prefix, suffix) {
  coef_path <- fixture_path(paste0(prefix, "_coef_", suffix, ".csv"))
  vcov_path <- fixture_path(paste0(prefix, "_vcov_", suffix, ".csv"))
  diag_path <- fixture_path(paste0(prefix, "_diagnostics_", suffix, ".csv"))
  skip_if(!file.exists(coef_path), "Fixture not found")

  stata_coef <- read.csv(coef_path, stringsAsFactors = FALSE)
  stata_coef$term[stata_coef$term == "_cons"] <- "(Intercept)"

  r_names <- names(coef(fit))
  m <- match(r_names, stata_coef$term)

  expect_equal(unname(coef(fit)), stata_coef$estimate[m],
               tolerance = stata_tol$coef,
               info = paste(suffix, "coefficients"))
  expect_equal(unname(sqrt(diag(fit$vcov))), stata_coef$std_error[m],
               tolerance = stata_tol$se,
               info = paste(suffix, "std errors"))

  # VCV
  if (file.exists(vcov_path)) {
    V_stata <- unname(as.matrix(read.csv(vcov_path)))
    stata_to_r <- match(stata_coef$term, r_names)
    V_r_reordered <- unname(fit$vcov[stata_to_r, stata_to_r])
    for (i in seq_len(nrow(V_stata))) {
      for (j in seq_len(ncol(V_stata))) {
        expect_equal(V_r_reordered[i, j], V_stata[i, j],
                     tolerance = stata_tol$vcov,
                     info = paste(suffix, "VCV[", i, ",", j, "]"))
      }
    }
  }

  # Diagnostics
  if (file.exists(diag_path)) {
    stata_diag <- read.csv(diag_path)
    diag <- fit$diagnostics

    expect_identical(fit$nobs, as.integer(stata_diag$N))

    if (!is.na(stata_diag$overid_stat) && stata_diag$overid_stat != 0) {
      expect_equal(diag$overid$stat, stata_diag$overid_stat,
                   tolerance = stata_tol$stat,
                   info = paste(suffix, "overid stat"))
      expect_equal(diag$overid$p, stata_diag$overid_p,
                   tolerance = stata_tol$pval,
                   info = paste(suffix, "overid p"))
    }
  }
}

# ============================================================================
# wmatrix alone (gmmw) + auto-bw
# ============================================================================

test_that("wmatrix + auto Bartlett: bandwidth uses W-step residuals", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  stata_bw <- read_autobiw("ts_wmat_autobiw", "ident_bartlett")
  skip_if(is.na(stata_bw), "Auto-bw fixture not found")

  W <- diag(4)
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = "auto", tvar = "t",
                wmatrix = W)
  expect_equal(fit$bw, stata_bw, info = "wmatrix auto Bartlett bw")
  check_fixture(fit, "ts_wmat_autobiw", "ident_bartlett")
})

test_that("wmatrix + auto Parzen: bandwidth uses W-step residuals", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  stata_bw <- read_autobiw("ts_wmat_autobiw", "ident_parzen")
  skip_if(is.na(stata_bw), "Auto-bw fixture not found")

  W <- diag(4)
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "parzen", bw = "auto", tvar = "t",
                wmatrix = W)
  expect_equal(fit$bw, stata_bw, info = "wmatrix auto Parzen bw")
  check_fixture(fit, "ts_wmat_autobiw", "ident_parzen")
})

# ============================================================================
# wmatrix + gmm2s + auto-bw
# ============================================================================

test_that("wmatrix + gmm2s + auto Bartlett: bandwidth uses W-step residuals", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  stata_bw <- read_autobiw("ts_wmat_autobiw", "gmm2s_bartlett")
  skip_if(is.na(stata_bw), "Auto-bw fixture not found")

  W <- diag(4)
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "bartlett", bw = "auto", tvar = "t",
                wmatrix = W, method = "gmm2s")
  expect_equal(fit$bw, stata_bw, info = "wmatrix+gmm2s auto Bartlett bw")
  check_fixture(fit, "ts_wmat_autobiw", "gmm2s_bartlett")
})

test_that("wmatrix + gmm2s + auto Parzen: bandwidth uses W-step residuals", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  stata_bw <- read_autobiw("ts_wmat_autobiw", "gmm2s_parzen")
  skip_if(is.na(stata_bw), "Auto-bw fixture not found")

  W <- diag(4)
  fit <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                vcov = "HAC", kernel = "parzen", bw = "auto", tvar = "t",
                wmatrix = W, method = "gmm2s")
  expect_equal(fit$bw, stata_bw, info = "wmatrix+gmm2s auto Parzen bw")
  check_fixture(fit, "ts_wmat_autobiw", "gmm2s_parzen")
})

# ============================================================================
# Verify W-step residuals differ from 2SLS residuals
# ============================================================================

test_that("wmatrix auto-bw differs from non-wmatrix auto-bw (Bartlett)", {
  skip_if(!file.exists(hac_data_path), "HAC data not found")
  stata_bw_wmat <- read_autobiw("ts_wmat_autobiw", "ident_bartlett")
  stata_bw_nowmat <- read_autobiw("ts_hac", "auto_bartlett")
  skip_if(is.na(stata_bw_wmat) || is.na(stata_bw_nowmat), "Fixtures not found")

  # Assert fixtures differ — this is the precondition we depend on
  expect_false(stata_bw_wmat == stata_bw_nowmat,
               info = "Fixtures should differ: wmatrix vs non-wmatrix auto-bw")

  # Verify R reproduces the difference
  W <- diag(4)
  fit_wmat <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                     vcov = "HAC", kernel = "bartlett", bw = "auto",
                     tvar = "t", wmatrix = W)
  fit_nowmat <- ivreg2(y ~ w | x | z1 + z2, data = ts_data,
                       vcov = "HAC", kernel = "bartlett", bw = "auto",
                       tvar = "t", method = "gmm2s")
  expect_false(fit_wmat$bw == fit_nowmat$bw,
               info = "wmatrix vs non-wmatrix auto-bw should differ")
})
