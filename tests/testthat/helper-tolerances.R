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

# Worst observed values below are from the end-of-grind tolerance re-audit
# (2026-07-06, 2022 comparisons across 36 configs on the re-based fixtures).
stata_tol <- list(
  coef = 1e-6,   # coefficients        (worst observed: 7.6e-7, a CUE cell)
  se   = 1e-6,   # standard errors     (worst observed: 2.0e-7, a CUE cell)
  vcov = 1e-6,   # VCV matrix elements (worst observed: 4.1e-7, a CUE cell)
  stat = 1e-4,   # test statistics (F, chi-sq, LM, J) — worst observed: 8.0e-7
  pval = 1e-4    # p-values (absolute) — worst observed: 2.9e-6
)

# The per-file exception tolerances (cue 5e-6, center 2e-5, dofminus 1e-5)
# were retired at the end-of-grind re-audit (2026-07-06): every binding case
# lived on fixture cells deleted in the Tranche 3 re-base (the card CUE and
# CUE+center cells, the M=2 card cluster dofminus cells, the card-era
# LIML+cluster cells). Two documented per-site overrides remain, each with
# its rationale at the call site:
#
#   test-cue.R / test-ts-operators.R: klein CUE coef/vcov = 1e-4
#     (H74 ruling: klein N=21 CUE optimizer noise)
#   test-user-matrices.R: wmatrix-vs-standard GMM fit comparison, vcov = 1e-5
#     (R-internal fit-vs-fit identity, not Stata parity)
#
# Run `Rscript pkg/tests/audit-tolerances.R` to measure actual discrepancies.
# Escalation rule: a breach of the standard tolerances is investigate-first;
# never widen a tolerance without a root cause.
