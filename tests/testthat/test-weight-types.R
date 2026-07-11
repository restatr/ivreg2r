# ============================================================================
# Tests: Frequency weights (fweight) and Probability weights (pweight)
# Ticket K2
# ============================================================================

# read_coef_fixture and read_diagnostics come from helper-fixtures.R
# card_wt (card data + deterministic region/fwt columns) comes from
# helper-fixtures.R; card-based tests use it directly (bundled data, full
# precision — see the M-11 re-base note in helper-fixtures.R).


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
  # The pweight-forces-robust warning fires (vcov defaults to "iid") before
  # the weight-positivity check errors out.
  expect_warning(
    expect_error(
      ivreg2(mpg ~ wt + hp, data = d, weights = w, weight_type = "pweight"),
      'strictly positive'
    ),
    'pweight implies robust VCE'
  )
})

test_that("pweight + iid overrides to HC0 with a warning", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  expect_warning(
    fit <- ivreg2(mpg ~ wt + hp, data = d, weights = w,
                  weight_type = "pweight", vcov = "iid"),
    'pweight implies robust VCE'
  )
  expect_equal(fit$vcov_type, "robust")
})

test_that("pweight + iid + small overrides to robust with a warning", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  expect_warning(
    fit <- ivreg2(mpg ~ wt + hp, data = d, weights = w,
                  weight_type = "pweight", vcov = "iid", small = TRUE),
    'pweight implies robust VCE'
  )
  expect_equal(fit$vcov_type, "robust")
})

test_that("pweight + explicit AC promotes to HAC with a warning (Stata parity)", {
  # Stata forces robust unconditionally for [pw=] (ivreg2.ado:353-357), so
  # pw + kernel yields HAC even when the user asked for AC. Verified against
  # Stata: educ SE matches to all printed digits (2026-06-10 session).
  d <- mtcars
  d$t <- seq_len(nrow(d))
  expect_warning(
    fit <- ivreg2(mpg ~ wt + hp, data = d, weights = disp,
                  weight_type = "pweight", vcov = "AC",
                  kernel = "bartlett", bw = 3, tvar = "t"),
    'pweight implies robust VCE'
  )
  expect_equal(fit$vcov_type, "HAC")
})

test_that("pweight + kiefer yields heteroskedasticity-robust Kiefer with a warning", {
  # Stata: [pw=] + kiefer reports "robust to heteroskedasticity and
  # within-cluster autocorrelation (Kiefer)"; our equivalent is the
  # kiefer-implied AC promoted to HAC. Verified against Stata numerically.
  d <- mtcars
  d$id <- rep(1:8, each = 4)
  d$t <- rep(1:4, times = 8)
  expect_warning(
    fit <- ivreg2(mpg ~ wt + hp, data = d, weights = disp,
                  weight_type = "pweight", kiefer = TRUE,
                  tvar = "t", ivar = "id"),
    'pweight implies robust VCE'
  )
  expect_equal(fit$vcov_type, "HAC")
  expect_true(fit$kiefer)
})

test_that("pweight + kiefer IV fit uses robust (KP) identification tests", {
  # Stata's ranktest receives the pweight-forced robust flag (but not the
  # truncated kernel) under kiefer, so the id diagnostics are KP rk
  # statistics with no kernel. Plain kiefer stays on the iid (Anderson/CD)
  # path. Both branches verified against Stata e(idstat)/e(widstat) to all
  # printed digits, 2026-06-10 session. The synthetic panel index mirrors
  # the Stata cross-check exactly; only the dispatch is asserted here.
  d <- card_wt
  d$id2 <- ceiling(seq_len(nrow(d)) / 10)
  d$tt <- (seq_len(nrow(d)) - 1) %% 10 + 1

  expect_warning(
    fit_pw <- ivreg2(lwage ~ exper + expersq | educ | nearc4, data = d,
                     weights = weight, weight_type = "pweight",
                     kiefer = TRUE, tvar = "tt", ivar = "id2"),
    'pweight implies robust VCE'
  )
  expect_match(fit_pw$diagnostics$underid$test_name, "Kleibergen-Paap")
  expect_true(is.finite(fit_pw$diagnostics$underid$stat))
  expect_true(is.finite(fit_pw$diagnostics$weak_id_robust$stat))

  # The dispatch contract is "robust, no kernel": the kiefer fit's id tests
  # must coincide exactly with those of a plain robust pweight fit of the
  # same model (Stata's ranktest receives only the robust flag under
  # kiefer). This also rules out the NA error-fallback path, which shares
  # the KP labels. vcov is explicit "robust" here, so pweight has nothing
  # to override and no warning fires.
  fit_rob <- ivreg2(lwage ~ exper + expersq | educ | nearc4, data = d,
                    weights = weight, weight_type = "pweight", vcov = "robust")
  expect_equal(fit_pw$diagnostics$underid$stat,
               fit_rob$diagnostics$underid$stat)
  expect_equal(fit_pw$diagnostics$weak_id_robust$stat,
               fit_rob$diagnostics$weak_id_robust$stat)

  fit_plain <- ivreg2(lwage ~ exper + expersq | educ | nearc4, data = d,
                      kiefer = TRUE, tvar = "tt", ivar = "id2")
  expect_match(fit_plain$diagnostics$underid$test_name, "Anderson")
  expect_true(is.finite(fit_plain$diagnostics$underid$stat))
  expect_null(fit_plain$diagnostics$weak_id_robust)
})

test_that("pweight + cluster does not override vcov", {
  # clusters is set, so the pweight-forces-robust branch (which requires
  # is.null(clusters)) never runs and no warning fires.
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
                weight_type = "pweight", clusters = ~ cyl)
  expect_equal(fit$vcov_type, "CL")
})

test_that("pweight + explicit HC0 matches pweight + iid override", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  expect_warning(
    fit_override <- ivreg2(mpg ~ wt + hp, data = d, weights = w,
                           weight_type = "pweight", vcov = "iid"),
    'pweight implies robust VCE'
  )
  fit_explicit <- ivreg2(mpg ~ wt + hp, data = d, weights = w,
                         weight_type = "pweight", vcov = "robust")
  expect_equal(vcov(fit_explicit), vcov(fit_override))
  expect_equal(coef(fit_explicit), coef(fit_override))
  expect_equal(fit_explicit$sigma, fit_override$sigma)
})

test_that("pweight + explicit HC1 + small matches pweight + iid + small override", {
  d <- mtcars
  d$w <- runif(nrow(d), 1, 5)
  expect_warning(
    fit_override <- ivreg2(mpg ~ wt + hp, data = d, weights = w,
                           weight_type = "pweight", vcov = "iid", small = TRUE),
    'pweight implies robust VCE'
  )
  fit_explicit <- ivreg2(mpg ~ wt + hp, data = d, weights = w,
                         weight_type = "pweight", vcov = "robust", small = TRUE)
  expect_equal(vcov(fit_explicit), vcov(fit_override))
  expect_equal(coef(fit_explicit), coef(fit_override))
  expect_equal(fit_explicit$sigma, fit_override$sigma)
})

test_that("pweight + kernel + iid overrides to HAC (not AC), with a warning", {
  phil <- na.omit(ivreg2r::phillips)
  expect_warning(
    fit <- ivreg2(cinf ~ unem, data = phil, weights = unem,
                  weight_type = "pweight", vcov = "iid",
                  kernel = "bartlett", bw = 3, tvar = "year"),
    'pweight implies robust VCE'
  )

  expect_equal(fit$vcov_type, "HAC")
})

test_that("pweight + kernel + iid + small overrides to HAC with correction, with a warning", {
  phil <- na.omit(ivreg2r::phillips)
  expect_warning(
    fit <- ivreg2(cinf ~ unem, data = phil, weights = unem,
                  weight_type = "pweight", vcov = "iid",
                  kernel = "bartlett", bw = 3, tvar = "year", small = TRUE),
    'pweight implies robust VCE'
  )
  expect_equal(fit$vcov_type, "HAC")
  expect_true(fit$small)
})

test_that("pweight + kernel + iid + small=FALSE overrides to HAC without correction, with a warning", {
  phil <- na.omit(ivreg2r::phillips)
  expect_warning(
    fit <- ivreg2(cinf ~ unem, data = phil, weights = unem,
                  weight_type = "pweight", vcov = "iid",
                  kernel = "bartlett", bw = 3, tvar = "year", small = FALSE),
    'pweight implies robust VCE'
  )
  expect_equal(fit$vcov_type, "HAC")
  expect_false(fit$small)
})

test_that("pweight + explicit robust + kernel routes to HAC", {
  phil <- na.omit(ivreg2r::phillips)
  fit <- ivreg2(cinf ~ unem, data = phil, weights = unem,
                weight_type = "pweight", vcov = "robust",
                kernel = "bartlett", bw = 3, tvar = "year")
  expect_equal(fit$vcov_type, "HAC")
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
  expect_warning(
    fit <- ivreg2(y ~ x, data = d, weights = w, weight_type = "pweight"),
    'pweight implies robust VCE'
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
  fit1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card_wt, weights = weight, vcov = "robust")
  fit2 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                 data = card_wt, weights = weight, weight_type = "aweight",
                 vcov = "robust")
  expect_equal(coef(fit1), coef(fit2))
  expect_equal(vcov(fit1), vcov(fit2))
  expect_equal(fit1$diagnostics$weak_id$stat, fit2$diagnostics$weak_id$stat)
})


# ============================================================================
# Section 4: Frequency weight (fweight) parity vs Stata fixtures (card_fw_*)
# ============================================================================
# card_fw: lwage ~ exper+expersq+black+south | educ | nearc4+nearc2,
# fweight fwt = mod(age,5)+1, cluster(region) M=9 (family M-11, re-based).
# `first endog(educ)` populates estat/arf/sstat on every cell.

test_fweight_config <- function(suffix, vcov_arg, small_arg, cluster_arg) {
  coef_path <- fixture_path(paste0("card_fw_coef_", suffix, ".csv"))
  vcov_path <- fixture_path(paste0("card_fw_vcov_", suffix, ".csv"))
  diag_path <- fixture_path(paste0("card_fw_diagnostics_", suffix, ".csv"))

  skip_if(!file.exists(coef_path))

  stata_coef <- read_coef_fixture(coef_path)
  stata_vcov <- read_vcov_fixture(vcov_path)
  stata_diag <- read_diagnostics(diag_path)

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card_wt, weights = fwt, weight_type = "fweight",
                vcov = vcov_arg, small = small_arg, clusters = cluster_arg)

  # N (fweight uses expanded N = sum(fwt))
  expect_equal(nobs(fit), as.integer(stata_diag$N))

  expect_stata_parity_core(fit, stata_coef, stata_vcov, stata_diag)

  # R-squared
  expect_equal(fit$r.squared, stata_diag$r2, tolerance = stata_tol$coef)

  # Underidentification / weak-identification
  expect_equal(fit$diagnostics$underid$stat, stata_diag$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$weak_id$stat, stata_diag$cdf,
               tolerance = stata_tol$stat)

  # Anderson-Rubin and Stock-Wright (weak-instrument-robust)
  expect_equal(fit$diagnostics$anderson_rubin$f_stat, stata_diag$arf,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$stock_wright$stat, stata_diag$sstat,
               tolerance = stata_tol$stat)

  # Endogeneity (C-statistic), from `endog(educ)`
  expect_equal(fit$diagnostics$endogeneity$stat, stata_diag$estat,
               tolerance = stata_tol$stat)

  # Overidentification: Sargan under iid, Hansen J under robust/cluster
  if (!is.na(stata_diag$sargan)) {
    expect_equal(fit$diagnostics$overid$stat, stata_diag$sargan,
                 tolerance = stata_tol$stat)
  }
  if (!is.na(stata_diag$j)) {
    expect_equal(fit$diagnostics$overid$stat, stata_diag$j,
                 tolerance = stata_tol$stat)
  }
}

for (cell in list(
  list(suffix = "iid",       vcov_arg = "iid",    small_arg = FALSE, cluster_arg = NULL),
  list(suffix = "iid_small", vcov_arg = "iid",    small_arg = TRUE,  cluster_arg = NULL),
  list(suffix = "hc1",       vcov_arg = "robust", small_arg = FALSE, cluster_arg = NULL),
  list(suffix = "hc1_small", vcov_arg = "robust", small_arg = TRUE,  cluster_arg = NULL),
  list(suffix = "cl",        vcov_arg = "iid",    small_arg = FALSE, cluster_arg = ~ region),
  list(suffix = "cl_small",  vcov_arg = "iid",    small_arg = TRUE,  cluster_arg = ~ region)
)) {
  test_that(paste("fweight overid matches Stata:", cell$suffix), {
    test_fweight_config(cell$suffix, cell$vcov_arg, cell$small_arg, cell$cluster_arg)
  })
}

test_that("fweight identity: weights = fwt equals unweighted fit on duplicated rows", {
  # Mirrors the .do generator's self-check: [fw=fwt] on the micro data must
  # equal the unweighted fit on the row-duplicated (expand fwt) data.
  card_expanded <- card_wt[rep(seq_len(nrow(card_wt)), card_wt$fwt), ]

  fit_fw <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                   data = card_wt, weights = fwt, weight_type = "fweight")
  fit_expanded <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                         data = card_expanded)

  expect_equal(coef(fit_fw), coef(fit_expanded))
  expect_equal(vcov(fit_fw), vcov(fit_expanded))
  expect_equal(nobs(fit_fw), nobs(fit_expanded))

  # This restores fweight first-stage coverage via the equivalence chain
  # (unweighted first stage is pinned against Stata by the M-10/M-25
  # fixtures). Compare all first-stage statistic components except the
  # per-observation residuals/fitted.values, which differ in length (3010
  # rows for fit_fw vs 9510 for fit_expanded).
  fs_fw <- fit_fw$first_stage$educ
  fs_expanded <- fit_expanded$first_stage$educ
  fs_scalar_names <- setdiff(names(fs_fw), c("residuals", "fitted.values"))
  expect_equal(fs_fw[fs_scalar_names], fs_expanded[fs_scalar_names])
})


# ============================================================================
# Section 5: Probability weight (pweight) parity vs Stata fixtures (card_pw_*)
# ============================================================================
# card_pw: lwage ~ exper+expersq+black+south | educ | nearc4+nearc2,
# pweight = genuine NLS sampling weight (pweight forces robust; no iid
# cells), cluster(region) M=9 (family M-11, re-based).
# The cl cells below pass vcov = "robust" together with clusters: pweight +
# vcov = "iid" emits the override message, and the CL dispatch is selected
# by the clusters argument regardless of vcov (this differs cosmetically
# from card_fw's cl cells, which use vcov = "iid" per the older convention).

test_pw_style_config <- function(suffix, weight_type_arg, vcov_arg, small_arg,
                                 cluster_arg) {
  coef_path <- fixture_path(paste0("card_pw_coef_", suffix, ".csv"))
  vcov_path <- fixture_path(paste0("card_pw_vcov_", suffix, ".csv"))
  diag_path <- fixture_path(paste0("card_pw_diagnostics_", suffix, ".csv"))

  skip_if(!file.exists(coef_path))

  stata_coef <- read_coef_fixture(coef_path)
  stata_vcov <- read_vcov_fixture(vcov_path)
  stata_diag <- read_diagnostics(diag_path)

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
                data = card_wt, weights = weight, weight_type = weight_type_arg,
                vcov = vcov_arg, small = small_arg, clusters = cluster_arg)

  # N (pweight/aweight use physical N)
  expect_equal(nobs(fit), as.integer(stata_diag$N))

  expect_stata_parity_core(fit, stata_coef, stata_vcov, stata_diag)

  # Underidentification / weak-identification
  expect_equal(fit$diagnostics$underid$stat, stata_diag$idstat,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$weak_id$stat, stata_diag$cdf,
               tolerance = stata_tol$stat)

  # Anderson-Rubin and Stock-Wright (weak-instrument-robust)
  expect_equal(fit$diagnostics$anderson_rubin$f_stat, stata_diag$arf,
               tolerance = stata_tol$stat)
  expect_equal(fit$diagnostics$stock_wright$stat, stata_diag$sstat,
               tolerance = stata_tol$stat)

  # Endogeneity (C-statistic), from `endog(educ)`
  expect_equal(fit$diagnostics$endogeneity$stat, stata_diag$estat,
               tolerance = stata_tol$stat)

  # Overidentification: Hansen J (pweight/aweight+robust always route here)
  expect_equal(fit$diagnostics$overid$stat, stata_diag$j,
               tolerance = stata_tol$stat)

  invisible(fit)
}

for (cell in list(
  list(suffix = "hc0",       vcov_arg = "robust", small_arg = FALSE, cluster_arg = NULL),
  list(suffix = "hc0_small", vcov_arg = "robust", small_arg = TRUE,  cluster_arg = NULL),
  list(suffix = "cl",        vcov_arg = "robust", small_arg = FALSE, cluster_arg = ~ region),
  list(suffix = "cl_small",  vcov_arg = "robust", small_arg = TRUE,  cluster_arg = ~ region)
)) {
  test_that(paste("pweight overid matches Stata:", cell$suffix), {
    test_pw_style_config(cell$suffix, "pweight", cell$vcov_arg,
                         cell$small_arg, cell$cluster_arg)
  })
}


# ============================================================================
# aweight+robust == pweight (Stata byte-identity, R direct equality)
# ============================================================================
# Byte-identity verified in the retired card_aweight_overid vs
# card_pweight_overid fixtures (2026-07-05); Stata treats [pw=w] as [aw=w]
# with robust. No separate aweight fixtures are generated (D6) — the
# aweight fit is checked against the SAME card_pw_* fixture as the pweight
# fit, and then the two R fits are asserted directly equal (the M-22
# lesson: comparing both variants only against the fixture would loosen
# the invariant to Stata tolerance rather than pinning it exactly).

for (cell in list(
  list(suffix = "hc0",       vcov_arg = "robust", small_arg = FALSE, cluster_arg = NULL),
  list(suffix = "hc0_small", vcov_arg = "robust", small_arg = TRUE,  cluster_arg = NULL),
  list(suffix = "cl",        vcov_arg = "robust", small_arg = FALSE, cluster_arg = ~ region),
  list(suffix = "cl_small",  vcov_arg = "robust", small_arg = TRUE,  cluster_arg = ~ region)
)) {
  test_that(paste("aweight+robust == pweight, both match Stata:", cell$suffix), {
    fit_pw <- test_pw_style_config(cell$suffix, "pweight", cell$vcov_arg,
                                   cell$small_arg, cell$cluster_arg)
    fit_aw <- test_pw_style_config(cell$suffix, "aweight", cell$vcov_arg,
                                   cell$small_arg, cell$cluster_arg)
    expect_equal(coef(fit_aw), coef(fit_pw))
    expect_equal(vcov(fit_aw), vcov(fit_pw))
    expect_equal(fit_aw$sigma, fit_pw$sigma)

    # The retired Stata fixtures' first-stage CSVs were byte-identical
    # between [aw=w] robust and [pw=w] (verified 2026-07-05), so this pins
    # pweight first-stage behavior through the aweight coverage that M-25
    # owns. Both fits see identical 3010-row data, so whole-object equality
    # holds.
    expect_equal(fit_aw$first_stage, fit_pw$first_stage)
  })
}


# ============================================================================
# Section 7: weight_type stored on the fitted object
# ============================================================================

test_that("weight_type is stored on the fitted object for fweight", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card_wt, weights = fwt, weight_type = "fweight")
  expect_equal(fit$weight_type, "fweight")
})

test_that("weight_type='aweight' is stored on the fitted object by default", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$weight_type, "aweight")
})

test_that("weight_type is stored in fitted object", {
  fit_a <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp)
  expect_equal(fit_a$weight_type, "aweight")

  fit_f <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card_wt, weights = fwt, weight_type = "fweight")
  expect_equal(fit_f$weight_type, "fweight")

  expect_warning(
    fit_p <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp,
                    weight_type = "pweight"),
    'pweight implies robust VCE'
  )
  expect_equal(fit_p$weight_type, "pweight")
})

test_that("n_physical is stored for fweight", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card_wt, weights = fwt, weight_type = "fweight")
  expect_equal(fit$n_physical, 3010L)
  expect_true(nobs(fit) > fit$n_physical)
})

test_that("n_physical is NULL for non-fweight", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, weights = disp)
  expect_null(fit$n_physical)
})
