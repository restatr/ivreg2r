# Stata parity tolerances — single source of truth for fixture comparisons.
# See CLAUDE.md "Testing Tolerances" for the policy table and escalation rules.
#
# === Rationale ===
#
# Tolerances are set at the level where exceeding them indicates a likely bug
# (wrong formula, wrong scaling, transposed matrix), NOT calibrated to observed
# discrepancies. The evidence base:
#
# Empirical (Card dataset, just-identified 2SLS, iid VCV):
#   All three R packages (ivreg, estimatr, ivreg2r) agree with each other to
#   near machine epsilon (~1e-14 for coefs). All three differ from Stata by
#   the same amount: ~5e-8 for coefficients, ~4.5e-9 for SEs/VCV. This is a
#   systematic R-vs-Stata gap (different LAPACK builds / accumulation order),
#   not an implementation bug. Condition number kappa(X) ~ 1.3e3.
#
# Ecosystem practices:
#   - estimatr uses 1e-5 for SEs vs Stata, 3e-5 for F-stats, 1e-3 for
#     clustered diagnostics
#   - fixest uses 1e-7 for SEs vs Stata's reghdfe
#   - R's all.equal() default is sqrt(.Machine$double.eps) ~ 1.5e-8
#   - sqrt(eps) would FAIL for our R-vs-Stata coefficient comparison (~5e-8)
#
# Numerical analysis:
#   For condition number ~1e3 and machine epsilon ~2.2e-16, the theoretical
#   bound on cross-implementation disagreement is ~2e-13. The observed ~5e-8
#   gap is larger, likely from differences in LAPACK builds or data loading.
#   Anything above 1e-6 relative would be hard to explain as floating-point
#   noise and warrants investigation.

# Selective warning handler for M=2 cluster diagnostic warnings.
# Binary cluster variables (e.g., smsa with 2 unique values) make diagnostic
# matrices rank-deficient. The warnings are correct and expected, but
# suppressWarnings() would mask unrelated warnings. This handler muffles
# only rank-deficient/singular diagnostic warnings and lets others propagate.
muffle_rank_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("rank-deficient|singular|not computed", conditionMessage(w)))
        invokeRestart("muffleWarning")
    }
  )
}

# Worst observed values below are from the standing tolerance audit (2,022
# comparisons across 36 configs; remeasure via
# `Rscript pkg/tests/audit-tolerances.R` from the repo root).
# Keep in sync: the tolerance table below is quoted in CLAUDE.md ("Testing
# Tolerances") and in the public validation article
# (vignettes/articles/validation.Rmd) — update all three together.
stata_tol <- list(
  coef = 1e-6,   # coefficients        (worst observed: 9.0e-7, a CUE cell)
  se   = 1e-6,   # standard errors     (worst observed: 1.7e-8, a LIML cell)
  vcov = 1e-6,   # VCV matrix elements (worst observed: 5.8e-8, a LIML cell)
  stat = 1e-4,   # test statistics (F, chi-sq, LM, J) — worst observed: 5.0e-7
  pval = 1e-4    # p-values (absolute) — worst observed: 2.5e-7
)

# --- Documented tolerance overrides -----------------------------------------
# Every tolerance literal looser than standard_tol_ceiling must have a row
# here and a root-cause comment at its call site. test-tolerance-policy.R
# enforces this mechanically: each row must match the exact number of
# occurrences of that (file, value) pair, so this inventory can neither
# under- nor over-claim. Counts are per literal, not per line (one call can
# carry two tolerance arguments). To add an override: find the root cause,
# comment the call site, then add the row. Escalation rule: a breach of the
# standard tolerances is investigate-first; never widen a tolerance without a
# root cause. Known limitation, accepted: a loose tolerance passed through a
# local variable escapes the scan; don't do that.
# NOTE: this file is also source()d by pkg/tests/audit-tolerances.R — keep it
# free of top-level testthat calls.
standard_tol_ceiling <- 1e-6  # at or below this is standard-or-tighter for every class

tolerance_overrides <- data.frame(
  file = c("test-ts-operators.R",
           "test-cue.R",
           "test-user-matrices.R",
           "test-gmm2s.R"),
  value = c(1e-5, 1e-5, 1e-5, 1e-4),
  count = c(3L, 1L, 1L, 1L),
  reason = c(
    "klein CUE coef/se + CUE=LIML identity: cross-platform optimizer-endpoint noise on the N=21 model; ~10x the measured cross-platform worst",
    "klein CUE vcov: same optimizer-endpoint noise, quadratic in the coefficient gap",
    "wmatrix-vs-standard just-identified GMM fit-vs-fit identity: different first-step W perturbs step-2 residuals; the VCV squares the gap",
    "inverted assertion (expect_false(all.equal)): looseness strengthens the claim that GMM2S and 2SLS coefficients differ"
  ),
  stringsAsFactors = FALSE
)
