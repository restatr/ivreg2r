# Structural rank bound of the moment-covariance meat

A one-way cluster-robust meat is a sum of M rank-one outer products (one
per cluster), so its rank is at most M regardless of floating-point
behavior. For UNWEIGHTED models, centering subtracts the plain mean
score, which makes the M cluster sums add to exactly zero and drops the
bound to M - 1. For weighted models the M - 1 refinement does NOT hold:
[`.cl_scores()`](https://restatr.com/ivreg2r/reference/dot-cl_scores.md)
centers by subtracting a weighted mean, so the cluster sums need not add
to zero and the exact-arithmetic rank stays M (verified empirically per
weight type at the 2026-07-06 review: unweighted centered sums are zero
to 1e-16 and the meat has rank M - 1; aweight/fweight/ pweight centered
sums are O(1) and the meat has full rank M). When the bound is below L
the meat is singular by construction, and consumers that must invert it
can refuse deterministically.

## Usage

``` r
.cluster_rank_bound(cluster_vec, kernel, center = FALSE, weights = NULL)
```

## Arguments

- cluster_vec:

  Length-N cluster vector (one-way), list of 2 vectors (two-way), or
  NULL.

- kernel:

  Kernel name, or NULL.

- center:

  Logical: are the moment conditions centered?

- weights:

  Normalized weight vector, or NULL. The M - 1 centered refinement
  applies only when NULL.

## Value

Numeric rank bound (M, or M - 1 for unweighted centered), or `Inf` when
no structural bound applies.

## Details

Stata anchor: Stata detects this singularity NUMERICALLY (rankS from
invsym's rank determination) but then applies an explicit suppression
POLICY — with a robust/cluster/kernel VCE and rankS \< iv1_ct it sets j,
cstat, and estat to missing (ivreg2.ado:1793-1808), and the efficient
2-step GMM path exits with r(506) "matrix not positive definite"
(s_egmm, ivreg2.ado:5436-5440). This helper reproduces those outcomes
through the structural bound instead of numerics because the numeric
detectors proved BLAS-dependent on the rank-deficit-1 case: at the
2026-07-06 CI cycle the Griliches cluster(year) + partial cells
(H28/H29, M = 7 clusters vs L = 8 moments) were detected on macOS
(lambda_min/lambda_max ~ 2e-18) but missed on ubuntu/Windows, whose BLAS
leaves a larger rounding residue in the deficit eigenvalue. Stata runs
on a single Mata implementation and does not face this cross-platform
problem; the structural gate is our platform-stable mechanism for the
same outcomes. The Anderson-Rubin statistic is deliberately NOT gated:
Stata computes it even when rankS \< L (verified live 2026-07-06 — H28
with `first` reports arf = 98.6498, matching ours, while j is missing).

Known scope limit: two-way cluster and cluster-kernel (DK/Thompson)
meats can also be structurally rank-deficient, but no comparably simple
bound holds (the CGM meat is a sum AND difference of cluster meats), so
those shapes return Inf and the numeric checks remain their detectors.
Extend here if a fixture ever exercises that regime.
