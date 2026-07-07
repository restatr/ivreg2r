# Compute L x L moment covariance matrix Omega (non-iid VCE types)

Computes the instrument-space score covariance for the Hansen J test and
the stored `fit$S`. This is a *new* computation in Z-space, not the K x
K meat from vcov-robust.R. Dispatched per VCE type: cluster,
cluster+kernel (DK/Thompson), AC, HAC, Stock-Watson, and plain HC. The
iid case lives in
[`.compute_moment_cov()`](https://restatr.com/ivreg2r/reference/dot-compute_moment_cov.md),
which shares this parameter list.

## Usage

``` r
.compute_omega(
  Z,
  residuals,
  weights,
  cluster_vec,
  N,
  dofminus = 0L,
  weight_type = "aweight",
  kernel = NULL,
  bw = NULL,
  time_index = NULL,
  center = FALSE,
  psd = NULL,
  vcov_type = "HAC",
  ZwZ = NULL,
  sw = FALSE,
  ivar_vec = NULL
)
```

## Arguments

- Z:

  N x L instrument matrix.

- residuals:

  N x 1 residual vector.

- weights:

  Normalized weights (sum to N), or NULL.

- cluster_vec:

  Cluster membership vector, or NULL.

- N:

  Number of observations.

- dofminus:

  Integer: large-sample DoF adjustment (default 0). HC path divides by
  `N - dofminus` (Stata livreg2.do line 326); cluster path divides by
  `N` (line 545, no dofminus adjustment).

- kernel:

  Canonical kernel name, or NULL for non-HAC.

- bw:

  Numeric bandwidth, or NULL.

- time_index:

  List from
  [`.build_time_index()`](https://restatr.com/ivreg2r/reference/dot-build_time_index.md),
  or NULL.

## Value

L x L symmetric matrix Omega.
