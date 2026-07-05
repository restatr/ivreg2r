# ============================================================================
# Tests: Center Option (Ticket N4; M-18 fixture re-base)
#
# The M-18 fixture re-base moved the center-option Stata-parity cells off the Card fixtures onto the M-16 canonical bases (planning/22-spec-matrix.md): griliches H06 (robust, robust+small, gmm2s+robust, dofminus, endog, orthog option-variations), abdata H88 (cluster, cluster+small, gmm2s+cluster), and phillips (kernel-without-robust HAC) -- all D5a option-variations, since Stata documents center with no worked example.
#
# Both CUE x center cells were DELETED, not ported: CUE x center coverage moves to M-17's canonical bases (the retired CUE+cluster+center cell was the known-dirty ubuntu CI basin-pathology cell, and CUE robust+center was the same card-CUE known-dirty class). The just-identified and (vacuous) card orthog/endog center cells were likewise retired -- center enters only through the shared meat helpers, so the identification-specific card cells added no distinct code path; the griliches endog/orthog cells below pin the C-statistics under center for the first time.
# ============================================================================

data(griliches, package = "ivreg2r")
data(abdata, package = "ivreg2r")
data(phillips, package = "ivreg2r")

# gril_formula (H06, help.txt:1154), ab_formula (H88), and phil_formula are provided by helper-fixtures.R.


# ============================================================================
# 1. Input validation (fixture-free; mtcars is fine here)
# ============================================================================

test_that("center must be TRUE or FALSE", {
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, center = 1),
    "`center` must be TRUE or FALSE"
  )
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, center = "yes"),
    "`center` must be TRUE or FALSE"
  )
  expect_error(
    ivreg2(mpg ~ wt + hp, data = mtcars, center = NA),
    "`center` must be TRUE or FALSE"
  )
})


# ============================================================================
# 2. Warning for no-op configurations (fixture-free)
# ============================================================================

test_that("center = TRUE + IID gives warning", {
  expect_warning(
    ivreg2(mpg ~ wt + hp, data = mtcars, center = TRUE),
    "center.*has no effect.*iid"
  )
})

test_that("center = TRUE + robust gives no warning", {
  expect_no_warning(
    ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "robust", center = TRUE)
  )
})

test_that("center = FALSE gives no warning (iid)", {
  expect_no_warning(
    ivreg2(mpg ~ wt + hp, data = mtcars, center = FALSE)
  )
})


# ============================================================================
# 3. Identity: center = FALSE identical to omitting center; center does not
#    move the point estimates (re-pointed to griliches H06)
# ============================================================================

test_that("center = FALSE gives identical results to default", {
  fit_default <- ivreg2(gril_formula, data = griliches, vcov = "robust")
  fit_false <- ivreg2(gril_formula, data = griliches, vcov = "robust",
                      center = FALSE)

  expect_identical(coef(fit_default), coef(fit_false))
  expect_identical(vcov(fit_default), vcov(fit_false))
  expect_identical(fit_default$sigma, fit_false$sigma)
  expect_identical(fit_default$diagnostics$overid$stat,
                   fit_false$diagnostics$overid$stat)
})

test_that("center does not change 2SLS coefficients, sigma, or rss", {
  fit_no <- ivreg2(gril_formula, data = griliches, vcov = "robust",
                   center = FALSE)
  fit_yes <- ivreg2(gril_formula, data = griliches, vcov = "robust",
                    center = TRUE)

  expect_identical(coef(fit_no), coef(fit_yes))
  expect_identical(fit_no$sigma, fit_yes$sigma)
  expect_identical(fit_no$rss, fit_yes$rss)
})


# ============================================================================
# 4. Property tests: where center DOES bite (griliches H06)
# ============================================================================

test_that("GMM2S: center changes the coefficients", {
  # center reaches the second-step GMM weighting matrix, so it moves the point estimates (unlike plain 2SLS where it touches only the VCE). Assert against a small positive floor rather than a bare inequality so a genuine (tiny but nonzero) shift is required.
  fit_no <- ivreg2(gril_formula, data = griliches, method = "gmm2s",
                   vcov = "robust", center = FALSE)
  fit_yes <- ivreg2(gril_formula, data = griliches, method = "gmm2s",
                    vcov = "robust", center = TRUE)

  max_abs_diff <- max(abs(coef(fit_yes) - coef(fit_no)))
  expect_gt(max_abs_diff, 1e-8)
})

test_that("endog() leaves estimation invariant under center", {
  # endog() is diagnostics-only and Stata does not forward center to the recursive endog call, so the point estimates and VCE are untouched (this is why the endog cell exports diagnostics only; M-16 pattern).
  fit_endog <- ivreg2(gril_formula, data = griliches, vcov = "robust",
                      center = TRUE, endog = "iq")
  fit_plain <- ivreg2(gril_formula, data = griliches, vcov = "robust",
                      center = TRUE)

  expect_identical(coef(fit_endog), coef(fit_plain))
  expect_identical(vcov(fit_endog), vcov(fit_plain))
})

test_that("orthog() leaves estimation invariant under center", {
  # orthog() is likewise diagnostics-only: the reported C-statistic changes but the fitted model does not.
  fit_orthog <- ivreg2(gril_formula, data = griliches, vcov = "robust",
                       center = TRUE, orthog = c("age", "mrt"))
  fit_plain <- ivreg2(gril_formula, data = griliches, vcov = "robust",
                      center = TRUE)

  expect_identical(coef(fit_orthog), coef(fit_plain))
  expect_identical(vcov(fit_orthog), vcov(fit_plain))
})


# ============================================================================
# 5. Storage and exposure (fixture-free)
# ============================================================================

test_that("center is stored in fit object", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "robust", center = TRUE)
  expect_true(fit$center)

  fit2 <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_false(fit2$center)
})

test_that("glance includes center column", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars, vcov = "robust", center = TRUE)
  gl <- glance(fit)
  expect_true("center" %in% names(gl))
  expect_true(gl$center)
})


# ============================================================================
# 6. Full-cell Stata parity (coef + vcov + diagnostics)
# ============================================================================
#
# Eight full cells across the three canonical bases. A table drives the repeated coef/vcov/diagnostics shell; the heterogeneous fit arguments (method, tvar/ivar/clusters, kernel/bw) live in each cell's fit_args and are spliced in with do.call, keeping every cell readable.
center_full_cells <- list(
  list(name = "griliches robust",
       prefix = "gril_center", suffix = "robust",
       fit_args = list(gril_formula, data = griliches,
                       vcov = "robust", center = TRUE)),
  list(name = "griliches robust small",
       prefix = "gril_center", suffix = "robust_small",
       fit_args = list(gril_formula, data = griliches,
                       vcov = "robust", small = TRUE, center = TRUE)),
  list(name = "griliches gmm2s robust",
       prefix = "gril_gmm2s_center", suffix = "robust",
       fit_args = list(gril_formula, data = griliches, method = "gmm2s",
                       vcov = "robust", center = TRUE)),
  list(name = "griliches robust dofminus(2)",
       prefix = "gril_center_dof", suffix = "robust",
       fit_args = list(gril_formula, data = griliches,
                       vcov = "robust", center = TRUE, dofminus = 2L)),
  list(name = "abdata cluster(id)",
       prefix = "ab_center", suffix = "cl",
       fit_args = list(ab_formula, data = abdata, tvar = "year", ivar = "id",
                       clusters = ~id, center = TRUE)),
  list(name = "abdata cluster(id) small",
       prefix = "ab_center", suffix = "cl_small",
       fit_args = list(ab_formula, data = abdata, tvar = "year", ivar = "id",
                       clusters = ~id, small = TRUE, center = TRUE)),
  list(name = "abdata gmm2s cluster(id)",
       prefix = "ab_gmm2s_center", suffix = "cl",
       fit_args = list(ab_formula, data = abdata, tvar = "year", ivar = "id",
                       method = "gmm2s", clusters = ~id, center = TRUE)),
  list(name = "phillips HAC bartlett bw=3",
       prefix = "phil_center", suffix = "hac_bw3",
       fit_args = list(phil_formula, data = phillips, tvar = "year",
                       vcov = "robust", kernel = "bartlett", bw = 3,
                       center = TRUE))
)

for (cell in center_full_cells) {
  test_that(paste0("center parity: ", cell$name), {
    coef_file <- paste0(cell$prefix, "_coef_", cell$suffix, ".csv")
    skip_if(!file.exists(fixture_path(coef_file)), "center fixture not found")

    fit <- do.call(ivreg2, cell$fit_args)
    expect_coef_fixture(fit, coef_file)
    expect_vcov_fixture(fit, paste0(cell$prefix, "_vcov_", cell$suffix, ".csv"))
    expect_diagnostics_fixture(fit, paste0(cell$prefix, "_diagnostics_", cell$suffix, ".csv"))
  })
}


# ============================================================================
# 7. endog() and orthog() C-statistics under center (small-invariance harness)
# ============================================================================
#
# These two families' Stata statistics (e(estat), e(cstat)) are invariant to `small`. The generator re-ran both cells with `small` added and asserted e(estat)/e(estatp) and e(cstat)/e(cstatp) reproduce to reldif < 1e-12 with exact df equality, so a single fixture backs both fits; the shared test_stata_fixture_cells() harness fits small = FALSE and small = TRUE and additionally pins the two fits' diagnostics equal at machine precision.

center_endog_cells <- list(
  list(name = "griliches robust center endog(iq)",
       fixture = "gril_center_endog_diagnostics_robust.csv",
       fit_args = list(gril_formula, data = griliches, vcov = "robust",
                       center = TRUE, endog = "iq"))
)

test_stata_fixture_cells(center_endog_cells, compare_endog_fixture,
                         slot = "endogeneity",
                         label_prefix = "center endogeneity matches Stata")

compare_center_orthog <- function(fit, fixture) {
  compare_orthog_fixture(fit, fixture)
  # Sanity that the C-statistic rides on the centered full-model J.
  expect_equal(fit$diagnostics$overid$stat, fixture$j,
               tolerance = stata_tol$stat, info = "overid J (centered)")
}

center_orthog_cells <- list(
  list(name = "griliches robust center orthog(age mrt)",
       fixture = "gril_center_orthog_robust.csv",
       fit_args = list(gril_formula, data = griliches, vcov = "robust",
                       center = TRUE, orthog = c("age", "mrt")))
)

test_stata_fixture_cells(center_orthog_cells, compare_center_orthog,
                         slot = "orthog",
                         label_prefix = "center orthogonality matches Stata")
