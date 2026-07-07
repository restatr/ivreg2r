# Compute J statistic given a pre-computed Omega matrix

Performs closed-form 2-step efficient GMM re-estimation using a provided
Omega (instrument-space score covariance) and returns the J statistic.
This is the inner loop shared by
[`.hansen_j_test()`](https://restatr.com/ivreg2r/reference/dot-hansen_j_test.md)
and
[`.compute_endogeneity_test()`](https://restatr.com/ivreg2r/reference/dot-compute_endogeneity_test.md).

## Usage

``` r
.compute_j_with_omega(Z, X, y, Omega, weights, N, rank_bound = Inf)
```

## Arguments

- Z:

  N x L instrument matrix.

- X:

  N x K regressor matrix.

- y:

  N x 1 response vector.

- Omega:

  L x L symmetric score covariance matrix.

- weights:

  Normalized weights (sum to N), or NULL.

- N:

  Number of observations.

- rank_bound:

  Structural rank bound of Omega from
  [`.cluster_rank_bound()`](https://restatr.com/ivreg2r/reference/dot-cluster_rank_bound.md),
  or `Inf` (default) when none applies. When below L, Omega is singular
  by construction and J is `NA` without consulting the
  platform-dependent numeric checks.

## Value

Scalar J statistic, or `NA_real_` if Omega is rank-deficient or the GMM
Hessian is singular.
