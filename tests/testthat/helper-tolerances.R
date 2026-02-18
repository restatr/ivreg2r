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

stata_tol <- list(
  coef = 1e-6,   # coefficients       (observed R-vs-Stata: ~5e-8)
  se   = 1e-6,   # standard errors    (observed R-vs-Stata: ~4.5e-9)
  vcov = 1e-6,   # VCV matrix elements (observed R-vs-Stata: ~9e-9)
  stat = 1e-4,   # test statistics (F, chi-sq, LM, J) — not yet measured
  pval = 1e-4    # p-values (absolute) — not yet measured
)
