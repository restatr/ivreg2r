# ============================================================================
# Tests: Continuously Updated GMM Estimator (CUE) — Ticket N3; M-17 fixture re-base
#
# The M-17 fixture re-base moved the CUE Stata-parity cells off the Card + synthetic seed-12345 fixtures onto the canonical bases (planning/22-spec-matrix.md): the stockwatson BSS 2007 pp. 476/480 worked CUE arc (Baum, Schaffer & Stillman 2007, Stata Journal 7(4):465-506, Section 4.3) carries the robust/HAC/small/center/dofminus/weighted/endog/orthog option-variations, and abdata H88 minus gmm2s carries the CUE x cluster and CUE x multi-endogenous intersections.
#
# Reuse rulings: klein H74 CUE-IID coefficient parity and the CUE = LIML identity live in test-ts-operators.R's H74 test — this file adds only the vcov and Hansen-J assertions that test does not make. mroz H68 b0 Stata parity lives in test-helpfile-examples.R (compare_hf asserts coef/SE/VCV/diagnostics) — this file keeps only fixture-free b0 behavioral tests.
#
# Delete rulings: all card CUE cells were retired as the known-dirty optimizer-basin class (the cluster(age) cells are also on the spec-matrix delete table); the card just-identified CUE cell was retired to the fixture-free just-id CUE = 2SLS identity below (2SLS Stata parity is owned by M-10); the non-deterministic [aw=age] weighted cell was replaced by the deterministic swwt cell; the synthetic seed-12345 ts DGP and its ts_gmm_data.csv export were retired in favor of the BSS p. 480 HAC cell.
# ============================================================================

data(stockwatson, package = "ivreg2r")
data(abdata, package = "ivreg2r")
data(klein, package = "ivreg2r")

# sw_formula (BSS 2007 p. 480), ab_formula (H88), klein_formula (H72-H76), and
# the stockwatson_swwt weighted data object are provided by helper-fixtures.R.

# Shared fits for the fixture-free behavioral sections below (CUE optimization is expensive; M-22 hoisting precedent). The fixture-comparison cells keep their own fits behind skip guards.
sw_cue_fit_iid <- ivreg2(sw_formula, data = stockwatson, method = "cue")
sw_cue_fit_robust <- ivreg2(sw_formula, data = stockwatson, method = "cue", vcov = "robust")

# CUE cells use the STANDARD tolerances since the end-of-grind re-audit
# (2026-07-06): the 5e-6 exception this file carried was calibrated on the
# retired card cells, and on the re-based cells the audit-measured worst CUE
# gap is 7.6e-7 coef / 2.0e-7 se / 4.1e-7 vcov against the 1e-6 standard.
# A future breach is investigate-first per the CLAUDE.md escalation rule --
# do not reintroduce a wider tolerance without a root cause. (The klein
# N=21 vcov cell below keeps its own documented 1e-4, the H74 ruling.)

# ============================================================================
# 1. Full-cell Stata parity (coef + vcov + diagnostics)
# ============================================================================
#
# Eight full cells across the two canonical bases. A table drives the repeated coef/vcov/diagnostics shell; the heterogeneous fit arguments (kernel/bw/tvar, weights, clusters, ivar) live in each cell's fit_args and are spliced in with do.call, keeping every cell readable.
cue_full_cells <- list(
  list(name = "stockwatson HAC bartlett bw=5 (BSS p.480)",
       prefix = "sw_cue", suffix = "hac_bw5",
       fit_args = list(sw_formula, data = stockwatson, method = "cue",
                       vcov = "robust", kernel = "bartlett", bw = 5,
                       tvar = "date")),
  list(name = "stockwatson robust",
       prefix = "sw_cue", suffix = "robust",
       fit_args = list(sw_formula, data = stockwatson, method = "cue",
                       vcov = "robust")),
  list(name = "stockwatson robust small",
       prefix = "sw_cue", suffix = "robust_small",
       fit_args = list(sw_formula, data = stockwatson, method = "cue",
                       vcov = "robust", small = TRUE)),
  list(name = "stockwatson robust center",
       prefix = "sw_cue_center", suffix = "robust",
       fit_args = list(sw_formula, data = stockwatson, method = "cue",
                       vcov = "robust", center = TRUE)),
  list(name = "stockwatson robust dofminus(2)",
       prefix = "sw_cue_dof", suffix = "robust",
       fit_args = list(sw_formula, data = stockwatson, method = "cue",
                       vcov = "robust", dofminus = 2L)),
  list(name = "stockwatson [aw=swwt] robust",
       prefix = "sw_cue_aw", suffix = "robust",
       fit_args = list(sw_formula, data = stockwatson_swwt, method = "cue",
                       vcov = "robust", weights = quote(swwt))),
  list(name = "abdata cluster(id)",
       prefix = "ab_cue", suffix = "cl",
       fit_args = list(ab_formula, data = abdata, tvar = "year", ivar = "id",
                       method = "cue", clusters = ~id)),
  list(name = "abdata cluster(id) small",
       prefix = "ab_cue", suffix = "cl_small",
       fit_args = list(ab_formula, data = abdata, tvar = "year", ivar = "id",
                       method = "cue", clusters = ~id, small = TRUE))
)

for (cell in cue_full_cells) {
  test_that(paste0("CUE parity: ", cell$name), {
    coef_file <- paste0(cell$prefix, "_coef_", cell$suffix, ".csv")
    skip_if(!file.exists(fixture_path(coef_file)), "CUE fixture not found")

    fit <- do.call(ivreg2, cell$fit_args)
    expect_coef_fixture(fit, coef_file)
    expect_vcov_fixture(fit, paste0(cell$prefix, "_vcov_", cell$suffix, ".csv"))
    expect_diagnostics_fixture(
      fit, paste0(cell$prefix, "_diagnostics_", cell$suffix, ".csv")
    )
    # For CUE the overid statistic is always the Hansen J (covers the abdata
    # cluster cell; asserted unconditionally since it holds for every cell).
    expect_equal(fit$diagnostics$overid$test_name, "Hansen J")
  })
}

# ============================================================================
# 2. klein CUE-IID (H74 reuse complement)
# ============================================================================

test_that("klein CUE-IID VCV and Hansen J match Stata (H74 complement)", {
  # The coef assertion and the CUE = LIML identity are owned by test-ts-operators.R's H74 test; this block adds only what that test does not assert: the full VCV and the Hansen J stat/p. The vcov tolerance is 1e-4 for the same reason H74 uses it on the coef fixture — klein N=21 CUE optimizer noise.
  skip_if(!file.exists(fixture_path("tsop_klein_vcov_cue.csv")),
          "klein CUE fixture not found")

  fit <- ivreg2(klein_formula, data = klein, tvar = "yr", method = "cue")
  expect_vcov_fixture(fit, "tsop_klein_vcov_cue.csv", tol = 1e-4)
  dx <- read_diagnostics(fixture_path("tsop_klein_diagnostics_cue.csv"))
  expect_equal(fit$diagnostics$overid$stat, dx$j,
               tolerance = stata_tol$stat, info = "Hansen J stat")
  expect_equal(fit$diagnostics$overid$p, dx$jp,
               tolerance = stata_tol$pval, info = "Hansen J p")
})

# ============================================================================
# 3. endog() and orthog() C-statistics under CUE (small-invariance harness)
# ============================================================================
#
# These two families' Stata statistics (e(estat), e(cstat)) are invariant to `small`. The generator (generate-cue-fixtures.do) re-ran both cells with `small` added and certified at generation time that e(estat)/e(estatp) and e(cstat)/e(cstatp) reproduce to reldif < 1e-12 with exact df equality, so a single fixture backs both fits; the shared test_stata_fixture_cells() harness fits small = FALSE and small = TRUE and additionally pins the two fits' diagnostics equal at machine precision.

cue_endog_cells <- list(
  list(name = "stockwatson robust CUE endog(UR)",
       fixture = "sw_cue_endog_diagnostics_robust.csv",
       fit_args = list(sw_formula, data = stockwatson, method = "cue",
                       vcov = "robust", endog = "UR"))
)

test_stata_fixture_cells(cue_endog_cells, compare_endog_fixture,
                         slot = "endogeneity",
                         label_prefix = "CUE endogeneity matches Stata")

compare_cue_orthog <- function(fit, fixture) {
  compare_orthog_fixture(fit, fixture)
  # Sanity that the C-statistic rides on the CUE full-model J.
  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat, info = "overid J (CUE)")
}

cue_orthog_cells <- list(
  list(name = "stockwatson robust CUE orthog(TBON_1)",
       fixture = "sw_cue_orthog_robust.csv",
       fit_args = list(sw_formula, data = stockwatson, method = "cue",
                       vcov = "robust", orthog = "TBON_1"))
)

test_stata_fixture_cells(cue_orthog_cells, compare_cue_orthog,
                         slot = "orthog",
                         label_prefix = "CUE orthogonality matches Stata")

# ============================================================================
# 4. Just-identified CUE = 2SLS (fixture-free identity)
# ============================================================================
#
# The retired card_justid_cue fixture is not replaced: for a just-identified model CUE coincides with 2SLS (the moment conditions hold exactly), and 2SLS Stata parity is owned by M-10, so a fixture here would be redundant.

test_that("just-identified CUE equals 2SLS with J = 0", {
  fit_2sls <- ivreg2(dinf ~ 1 | UR | ggdp_2, data = stockwatson,
                     vcov = "robust")
  fit_cue <- ivreg2(dinf ~ 1 | UR | ggdp_2, data = stockwatson,
                    method = "cue", vcov = "robust")

  expect_equal(coef(fit_cue), coef(fit_2sls), tolerance = 1e-6)
  expect_equal(fit_cue$diagnostics$overid$stat, 0)
  expect_true(is.na(fit_cue$diagnostics$overid$p))
  expect_identical(fit_cue$diagnostics$overid$df, 0L)
})

# ============================================================================
# 5. endog()/orthog() leave CUE estimation invariant (fixture-free)
# ============================================================================

test_that("endog() leaves CUE estimation invariant", {
  # endog() is diagnostics-only and does not enter the CUE objective, so the point estimates and VCE are untouched; the optimizer is deterministic from identical inputs, so the comparison is exact (this is why the endog cell exports diagnostics only; M-16/M-18 pattern).
  fit_endog <- ivreg2(sw_formula, data = stockwatson, method = "cue",
                      vcov = "robust", endog = "UR")

  expect_identical(coef(fit_endog), coef(sw_cue_fit_robust))
  expect_identical(vcov(fit_endog), vcov(sw_cue_fit_robust))
})

test_that("orthog() leaves CUE estimation invariant", {
  # orthog() is likewise diagnostics-only: the reported C-statistic re-optimizes a restricted model, but the fitted full model does not change.
  fit_orthog <- ivreg2(sw_formula, data = stockwatson, method = "cue",
                       vcov = "robust", orthog = "TBON_1")

  expect_identical(coef(fit_orthog), coef(sw_cue_fit_robust))
  expect_identical(vcov(fit_orthog), vcov(sw_cue_fit_robust))
})

# ============================================================================
# 6. Stock-Yogo uses LIML tables for CUE
# ============================================================================

test_that("CUE uses LIML Stock-Yogo tables", {
  fit <- sw_cue_fit_iid

  sy <- fit$diagnostics$weak_id_sy
  expect_false(is.null(sy))
  # CUE should use LIML size tables
  expect_true("LIML size" %in% sy$type)
  # Should NOT have IV size tables (those are for 2SLS/GMM2S)
  expect_false("IV size" %in% sy$type)
})

# ============================================================================
# 7. CUE overid test name is always Hansen J
# ============================================================================

test_that("CUE overid is always Hansen J", {
  # IID and robust via the shared fits; the cluster case is asserted in the
  # section-1 fixture loop below (for CUE the overid test_name is always
  # "Hansen J", so the loop asserts it unconditionally for every cell).
  expect_equal(sw_cue_fit_iid$diagnostics$overid$test_name, "Hansen J")
  expect_equal(sw_cue_fit_robust$diagnostics$overid$test_name, "Hansen J")
})

# ============================================================================
# 8. CUE diagnostics populated
# ============================================================================

test_that("CUE populates identification diagnostics", {
  fit <- sw_cue_fit_robust

  # Overid
  expect_false(is.null(fit$diagnostics$overid))

  # Underid
  expect_false(is.null(fit$diagnostics$underid))

  # Weak-id
  expect_false(is.null(fit$diagnostics$weak_id))
  expect_false(is.null(fit$diagnostics$weak_id_robust))

  # First-stage
  expect_false(is.null(fit$first_stage))

  # AR
  expect_false(is.null(fit$diagnostics$anderson_rubin))

  # Stock-Wright
  expect_false(is.null(fit$diagnostics$stock_wright))

  # Endogeneity
  expect_false(is.null(fit$diagnostics$endogeneity))
})

# ============================================================================
# 9. b0 mode — evaluate at a fixed parameter vector (fixture-free)
# ============================================================================
#
# The old "b0 mode J(b0) matches Stata fixture" test and its card_overid_b0_vector.csv read are DELETED: b0 Stata parity is owned by test-helpfile-examples.R's H68 test (mroz, compare_hf). This section keeps only the fixture-free behavioral contracts, re-pointed to the stockwatson model.

test_that("b0 mode produces coefficients equal to b0", {
  # Get 2SLS estimates to use as b0
  fit_2sls <- ivreg2(sw_formula, data = stockwatson, vcov = "robust")
  b0_vec <- coef(fit_2sls)

  fit_b0 <- ivreg2(sw_formula, data = stockwatson, vcov = "robust",
                   b0 = b0_vec)

  # Coefficients should be exactly b0
  expect_equal(coef(fit_b0), b0_vec, tolerance = 1e-14)
})

test_that("b0 mode suppresses identification diagnostics", {
  fit_2sls <- ivreg2(sw_formula, data = stockwatson, vcov = "robust")
  b0_vec <- coef(fit_2sls)

  fit_b0 <- ivreg2(sw_formula, data = stockwatson, vcov = "robust",
                   b0 = b0_vec)

  # Overid should still be populated (J-stat)
  expect_false(is.null(fit_b0$diagnostics$overid))

  # All identification diagnostics should be NULL
  expect_null(fit_b0$diagnostics$underid)
  expect_null(fit_b0$diagnostics$weak_id)
  expect_null(fit_b0$diagnostics$weak_id_robust)
  expect_null(fit_b0$diagnostics$weak_id_sy)
  expect_null(fit_b0$diagnostics$anderson_rubin)
  expect_null(fit_b0$diagnostics$stock_wright)
  expect_null(fit_b0$diagnostics$endogeneity)

  # First stage should be NULL
  expect_null(fit_b0$first_stage)
})

test_that("b0 implies method = cue", {
  fit_2sls <- ivreg2(sw_formula, data = stockwatson, vcov = "robust")
  b0_vec <- coef(fit_2sls)

  fit_b0 <- ivreg2(sw_formula, data = stockwatson, vcov = "robust",
                   b0 = b0_vec)

  expect_equal(fit_b0$method, "cue")
})

test_that("b0 convergence is 0", {
  fit_2sls <- ivreg2(sw_formula, data = stockwatson, vcov = "robust")
  b0_vec <- coef(fit_2sls)

  fit_b0 <- ivreg2(sw_formula, data = stockwatson, vcov = "robust",
                   b0 = b0_vec)

  expect_equal(fit_b0$cue_convergence, 0L)
})

test_that("b0 unnamed vector works in column order", {
  fit_2sls <- ivreg2(sw_formula, data = stockwatson, vcov = "robust")
  b0_named <- coef(fit_2sls)
  b0_unnamed <- unname(b0_named)

  fit_named <- ivreg2(sw_formula, data = stockwatson, vcov = "robust",
                      b0 = b0_named)
  fit_unnamed <- ivreg2(sw_formula, data = stockwatson, vcov = "robust",
                        b0 = b0_unnamed)

  expect_equal(coef(fit_named), coef(fit_unnamed))
})

# ============================================================================
# 10. Input validation
# ============================================================================
#
# The stockwatson model has K = 2 coefficients (UR, intercept) and L = 5 instrument columns, so the b0 length checks use length-2 vectors (the bad-length test passes length 3) and the wmatrix/smatrix mutual-exclusion tests use diag(5).

test_that("CUE mutual exclusion with fuller", {
  expect_error(
    ivreg2(sw_formula, data = stockwatson, method = "cue", fuller = 1),
    "Cannot specify `fuller` with `method = \"cue\"`"
  )
})

test_that("CUE mutual exclusion with kclass", {
  expect_error(
    ivreg2(sw_formula, data = stockwatson, method = "cue", kclass = 0.5),
    "Cannot specify `kclass` with `method = \"cue\"`"
  )
})

test_that("CUE mutual exclusion with wmatrix", {
  W <- diag(5)
  expect_error(
    ivreg2(sw_formula, data = stockwatson, method = "cue", wmatrix = W),
    "Cannot specify `wmatrix` with `method = \"cue\"`"
  )
})

test_that("CUE mutual exclusion with smatrix", {
  S <- diag(5)
  expect_error(
    ivreg2(sw_formula, data = stockwatson, method = "cue", smatrix = S),
    "Cannot specify `smatrix` with `method = \"cue\"`"
  )
})

test_that("CUE requires IV model", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, method = "cue"),
    "requires an IV model"
  )
})

test_that("b0 mutual exclusion with gmm2s", {
  expect_error(
    ivreg2(sw_formula, data = stockwatson, method = "gmm2s", b0 = rep(0, 2)),
    "Cannot specify `b0` with `method = \"gmm2s\"`"
  )
})

test_that("b0 mutual exclusion with liml", {
  expect_error(
    ivreg2(sw_formula, data = stockwatson, method = "liml", b0 = rep(0, 2)),
    "Cannot specify `b0` with `method = \"liml\"`"
  )
})

test_that("b0 mutual exclusion with fuller", {
  expect_error(
    ivreg2(sw_formula, data = stockwatson, fuller = 1, b0 = rep(0, 2)),
    "Cannot specify `b0` with `fuller`"
  )
})

test_that("b0 mutual exclusion with wmatrix", {
  W <- diag(5)
  expect_error(
    ivreg2(sw_formula, data = stockwatson, b0 = rep(0, 2), wmatrix = W),
    "Cannot specify `b0` with `wmatrix`"
  )
})

test_that("b0 dimension check", {
  expect_error(
    ivreg2(sw_formula, data = stockwatson, b0 = c(1, 2, 3)),
    "`b0` length .* must equal"
  )
})

test_that("b0 must be numeric", {
  expect_error(
    ivreg2(sw_formula, data = stockwatson, b0 = "not numeric"),
    "`b0` must be a numeric vector"
  )
})

test_that("b0 must be finite", {
  expect_error(
    ivreg2(sw_formula, data = stockwatson, b0 = c(1, NA)),
    "`b0` must contain only finite values"
  )
})

test_that("b0 named vector with wrong names errors", {
  b0_bad <- c(a = 1, b = 2)
  expect_error(
    ivreg2(sw_formula, data = stockwatson, b0 = b0_bad),
    "not matching model columns"
  )
})

test_that("b0 + uppercase method errors correctly (case-insensitive)", {
  # Regression test: uppercase method must still trigger incompatibility errors
  expect_error(
    ivreg2(sw_formula, data = stockwatson, method = "LIML", b0 = rep(0, 2)),
    "Cannot specify `b0` with `method = \"liml\"`"
  )
  expect_error(
    ivreg2(sw_formula, data = stockwatson, method = "GMM2S", b0 = rep(0, 2)),
    "Cannot specify `b0` with `method = \"gmm2s\"`"
  )
})

test_that("b0 + uppercase '2SLS' auto-promotes to CUE", {
  fit_2sls <- ivreg2(sw_formula, data = stockwatson, vcov = "robust")
  b0_vec <- coef(fit_2sls)
  # method = "2SLS" (uppercase) should still auto-promote to CUE
  fit_b0 <- ivreg2(sw_formula, data = stockwatson, vcov = "robust",
                   method = "2SLS", b0 = b0_vec)
  expect_equal(fit_b0$method, "cue")
})

# ============================================================================
# 11. Convergence
# ============================================================================

test_that("CUE convergence code = 0 for well-specified model", {
  fit <- sw_cue_fit_iid

  expect_equal(fit$cue_convergence, 0L)
  expect_false(is.null(fit$cue_message))
})

test_that("CUE convergence fields are NULL for non-CUE", {
  fit <- ivreg2(sw_formula, data = stockwatson)

  expect_null(fit$cue_convergence)
  expect_null(fit$cue_message)
})

# ============================================================================
# 12. Broom methods
# ============================================================================

test_that("tidy works for CUE", {
  fit <- sw_cue_fit_iid
  td <- tidy(fit)

  expect_s3_class(td, "tbl_df")
  expect_equal(nrow(td), 2L)
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value") %in%
                  names(td)))
})

test_that("glance works for CUE and includes cue_convergence", {
  fit <- sw_cue_fit_iid
  gl <- glance(fit)

  expect_s3_class(gl, "tbl_df")
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$method, "cue")
  expect_true("cue_convergence" %in% names(gl))
  expect_equal(gl$cue_convergence, 0L)
})

test_that("glance cue_convergence is NA for non-CUE", {
  fit <- ivreg2(sw_formula, data = stockwatson)
  gl <- glance(fit)

  expect_true(is.na(gl$cue_convergence))
})

test_that("augment works for CUE", {
  fit <- sw_cue_fit_iid
  aug <- augment(fit)

  expect_s3_class(aug, "tbl_df")
  expect_true(all(c(".fitted", ".resid") %in% names(aug)))
})

# ============================================================================
# 13. Method field
# ============================================================================

test_that("CUE method field is 'cue'", {
  fit <- sw_cue_fit_iid

  expect_equal(fit$method, "cue")
})

# The old "CUE IID equals LIML" section is DELETED: the identity is owned by test-ts-operators.R's H74 test on klein (the help file's own CUE = LIML demonstration).

# ============================================================================
# 14. Summary/print methods
# ============================================================================

test_that("summary prints for CUE without error", {
  fit <- sw_cue_fit_iid

  out <- capture.output(summary(fit))
  expect_true(any(grepl("CUE Estimation", out)))
  expect_true(any(grepl("Estimates efficient for arbitrary homoskedasticity", out)))
})

test_that("summary prints for CUE with b0", {
  fit_2sls <- ivreg2(sw_formula, data = stockwatson, vcov = "robust")
  fit_b0 <- ivreg2(sw_formula, data = stockwatson, vcov = "robust",
                   b0 = coef(fit_2sls))

  out <- capture.output(summary(fit_b0))
  expect_true(any(grepl("User-supplied Parameter Vector", out)))
})

test_that("print method works for CUE", {
  fit <- sw_cue_fit_iid

  out <- capture.output(print(fit))
  expect_true(any(grepl("CUE Estimation", out)))
})
