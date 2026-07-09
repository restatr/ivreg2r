# ============================================================================
# Tests: Stock-Watson (2008) Panel-Robust VCE (Ticket Q3)
# ============================================================================

# --- Data loading ---
data(wagepan)

# Helper: check SW fixture against fitted model
check_sw_fixture <- function(fit, prefix) {
  coef_path <- fixture_path(paste0(prefix, "_coef.csv"))
  vcov_path <- fixture_path(paste0(prefix, "_vcov.csv"))
  diag_path <- fixture_path(paste0(prefix, "_diagnostics.csv"))
  skip_if(!file.exists(coef_path), "SW fixture not found")

  stata_coef <- read.csv(coef_path, stringsAsFactors = FALSE)
  stata_coef$term[stata_coef$term == "_cons"] <- "(Intercept)"

  r_names <- names(coef(fit))

  # Coefficients
  for (i in seq_len(nrow(stata_coef))) {
    nm <- stata_coef$term[i]
    if (nm %in% r_names) {
      expect_equal(unname(coef(fit)[nm]), stata_coef$estimate[i],
                   tolerance = stata_tol$coef,
                   info = paste("coef mismatch:", nm))
    }
  }

  # Standard errors
  r_se <- sqrt(diag(vcov(fit)))
  for (i in seq_len(nrow(stata_coef))) {
    nm <- stata_coef$term[i]
    if (nm %in% r_names) {
      expect_equal(unname(r_se[nm]), stata_coef$std_error[i],
                   tolerance = stata_tol$se,
                   info = paste("SE mismatch:", nm))
    }
  }

  # VCV matrix (match by name since R and Stata have different ordering)
  if (file.exists(vcov_path)) {
    stata_V <- as.matrix(read.csv(vcov_path))
    # Map Stata term names to R names
    stata_names <- stata_coef$term  # already translated _cons -> (Intercept)
    # Build a permutation index from Stata order to R order
    perm <- match(r_names, stata_names)
    if (!anyNA(perm)) {
      stata_V_reordered <- stata_V[perm, perm, drop = FALSE]
      r_V <- vcov(fit)
      for (i in seq_along(r_names)) {
        for (j in seq_along(r_names)) {
          expect_equal(unname(r_V[i, j]), unname(stata_V_reordered[i, j]),
                       tolerance = stata_tol$vcov,
                       info = paste("VCV mismatch:", r_names[i], r_names[j]))
        }
      }
    }
  }

  # Diagnostics
  if (file.exists(diag_path)) {
    stata_diag <- read.csv(diag_path, check.names = FALSE)
    d <- fit$diagnostics

    # Hansen J / overid
    if (!is.null(d$overid) && !is.na(stata_diag$overid_stat) &&
        d$overid$df > 0L) {
      expect_equal(d$overid$stat, stata_diag$overid_stat,
                   tolerance = stata_tol$stat,
                   info = "overid stat")
      expect_equal(d$overid$p, stata_diag$overid_p,
                   tolerance = stata_tol$pval,
                   info = "overid p")
    }

    # KP rk (uses standard HC, not SW — matches Stata behavior)
    if (!is.null(d$underid) && !is.na(stata_diag$underid_stat)) {
      expect_equal(d$underid$stat, stata_diag$underid_stat,
                   tolerance = stata_tol$stat,
                   info = "underid stat")
    }

    # KP Wald F
    if (!is.null(d$weak_id_robust) && !is.na(stata_diag$weak_id_kp_f)) {
      expect_equal(d$weak_id_robust$stat, stata_diag$weak_id_kp_f,
                   tolerance = stata_tol$stat,
                   info = "weak_id KP F")
    }

    # Anderson-Rubin
    if (!is.null(d$anderson_rubin) && !is.na(stata_diag$ar_f)) {
      expect_equal(d$anderson_rubin$f_stat, stata_diag$ar_f,
                   tolerance = stata_tol$stat,
                   info = "AR F stat")
      expect_equal(d$anderson_rubin$chi2_stat, stata_diag$ar_chi2,
                   tolerance = stata_tol$stat,
                   info = "AR chi2 stat")
    }

    # Stock-Wright S
    if (!is.null(d$stock_wright) && !is.na(stata_diag$sw_stat)) {
      expect_equal(d$stock_wright$stat, stata_diag$sw_stat,
                   tolerance = stata_tol$stat,
                   info = "Stock-Wright S stat")
    }

    # Endogeneity C-stat
    if (!is.null(d$endogeneity) && !is.na(stata_diag$endog_stat)) {
      expect_equal(d$endogeneity$stat, stata_diag$endog_stat,
                   tolerance = stata_tol$stat,
                   info = "endogeneity stat")
    }

    # Model F
    if (!is.null(fit$model_f) && !is.na(stata_diag$model_f)) {
      expect_equal(fit$model_f, stata_diag$model_f,
                   tolerance = stata_tol$stat,
                   info = "model F")
    }
  }
}

# ============================================================================
# Input validation tests
# ============================================================================

test_that("sw must be logical", {
  expect_error(
    ivreg2(lwage ~ exper | hours | educ, data = wagepan, sw = "yes", ivar = "nr"),
    "must be TRUE or FALSE"
  )
  expect_error(
    ivreg2(lwage ~ exper | hours | educ, data = wagepan, sw = NA, ivar = "nr"),
    "must be TRUE or FALSE"
  )
})

test_that("sw requires ivar", {
  expect_error(
    ivreg2(lwage ~ exper | hours | educ, data = wagepan, sw = TRUE),
    "requires panel data"
  )
})

test_that("sw blocks clustering", {
  expect_error(
    ivreg2(lwage ~ exper | hours | educ, data = wagepan,
           sw = TRUE, ivar = "nr", clusters = ~nr),
    "not supported with clustering"
  )
})

test_that("sw blocks kernel", {
  expect_error(
    ivreg2(lwage ~ exper | hours | educ, data = wagepan,
           sw = TRUE, ivar = "nr", kernel = "bartlett", bw = 3, tvar = "year"),
    "not supported with kernel"
  )
})

test_that("sw blocks fweight", {
  wagepan$fw <- 1L
  expect_error(
    ivreg2(lwage ~ exper | hours | educ, data = wagepan, weights = fw,
           sw = TRUE, ivar = "nr", weight_type = "fweight"),
    "fweights not supported"
  )
})

test_that("sw blocks iweight", {
  expect_error(
    ivreg2(lwage ~ exper | hours | educ, data = wagepan, weights = hours,
           sw = TRUE, ivar = "nr", weight_type = "iweight"),
    "iweights not supported"
  )
})

test_that("sw blocks kiefer", {
  expect_error(
    ivreg2(lwage ~ exper | hours | educ, data = wagepan,
           sw = TRUE, ivar = "nr", kiefer = TRUE, tvar = "year"),
    "not supported with kernel|not supported with Kiefer"
  )
})

test_that("sw blocks dkraay", {
  expect_error(
    ivreg2(lwage ~ exper | hours | educ, data = wagepan,
           sw = TRUE, ivar = "nr", dkraay = 3, tvar = "year"),
    "not supported with kernel|not supported with Driscoll-Kraay"
  )
})

test_that("sw forces robust VCE", {
  fit <- ivreg2(lwage ~ exper | hours | educ, data = wagepan,
                sw = TRUE, ivar = "nr")
  expect_equal(fit$vcov_type, "robust")
})

# ============================================================================
# Fixture-based parity tests
# ============================================================================

test_that("SW just-identified matches Stata", {
  fit <- ivreg2(lwage ~ exper | hours | educ, data = wagepan,
                sw = TRUE, ivar = "nr")
  check_sw_fixture(fit, "wp_sw_justid")
})

test_that("SW just-identified + small matches Stata", {
  fit <- ivreg2(lwage ~ exper | hours | educ, data = wagepan,
                sw = TRUE, ivar = "nr", small = TRUE)
  check_sw_fixture(fit, "wp_sw_justid_small")
})

test_that("SW overidentified matches Stata", {
  fit <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                sw = TRUE, ivar = "nr")
  check_sw_fixture(fit, "wp_sw_overid")
})

test_that("SW overidentified + small matches Stata", {
  fit <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                sw = TRUE, ivar = "nr", small = TRUE)
  check_sw_fixture(fit, "wp_sw_overid_small")
})

test_that("SW + aweight matches Stata", {
  wagepan$aw <- wagepan$hours + 10
  fit <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                weights = aw, sw = TRUE, ivar = "nr")
  check_sw_fixture(fit, "wp_sw_aweight")
})

test_that("SW + LIML matches Stata", {
  fit <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                sw = TRUE, ivar = "nr", method = "liml")
  check_sw_fixture(fit, "wp_sw_liml")
})

test_that("SW + dofminus matches Stata", {
  fit <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                sw = TRUE, ivar = "nr", dofminus = 1L)
  check_sw_fixture(fit, "wp_sw_dofminus")
})

test_that("SW + endogeneity test matches Stata", {
  fit <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                sw = TRUE, ivar = "nr", endog = "hours")
  check_sw_fixture(fit, "wp_sw_endog")
})

test_that("SW + center matches Stata", {
  fit <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                sw = TRUE, ivar = "nr", center = TRUE)
  check_sw_fixture(fit, "wp_sw_center")
})

test_that("SW + gmm2s matches Stata", {
  fit <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                sw = TRUE, ivar = "nr", method = "gmm2s")
  check_sw_fixture(fit, "wp_sw_gmm2s")
})

test_that("SW + cue matches Stata", {
  fit <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                sw = TRUE, ivar = "nr", method = "cue")
  check_sw_fixture(fit, "wp_sw_cue")
})

# ============================================================================
# SW + center + reduced_form = "system" (cross-equation blocks)
# ============================================================================

test_that("SW + center + system RF produces symmetric system VCV", {
  fit <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                sw = TRUE, ivar = "nr", center = TRUE,
                reduced_form = "system")
  rf <- fit$reduced_form
  expect_equal(rf$mode, "system")
  # System VCV must be symmetric
  expect_true(isSymmetric(rf$vcov, tol = sqrt(.Machine$double.eps)))
  # Dimensions: (1 + K1) equations * L instruments
  n_eq <- 1L + length(fit$endogenous)
  L <- fit$rankzz
  expect_equal(nrow(rf$vcov), n_eq * L)
})

test_that("SW + center + system RF matches non-center up to centering effect", {
  fit_c <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                  sw = TRUE, ivar = "nr", center = TRUE,
                  reduced_form = "system")
  fit_nc <- ivreg2(lwage ~ exper | hours | educ + married, data = wagepan,
                   sw = TRUE, ivar = "nr", center = FALSE,
                   reduced_form = "system")
  # Both should produce valid system VCVs with same dimensions
  expect_equal(dim(fit_c$reduced_form$vcov), dim(fit_nc$reduced_form$vcov))
  # Coefficients identical (centering only affects VCV, not point estimates)
  expect_equal(coef(fit_c), coef(fit_nc))
})

# ============================================================================
# Metadata tests
# ============================================================================

test_that("sw is stored on the return object", {
  fit <- ivreg2(lwage ~ exper | hours | educ, data = wagepan,
                sw = TRUE, ivar = "nr")
  expect_true(fit$sw)
})

test_that("sw is stored on the fitted object", {
  fit <- ivreg2(lwage ~ exper | hours | educ, data = wagepan,
                sw = TRUE, ivar = "nr")
  expect_true(fit$sw)
})

test_that("summary shows SW VCV description", {
  fit <- ivreg2(lwage ~ exper | hours | educ, data = wagepan,
                sw = TRUE, ivar = "nr")
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Stock-Watson", out)))
})
