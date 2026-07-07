# Compute cluster-robust VCV (one-way or two-way)

Sandwich estimator with scores aggregated within clusters via
[`.cluster_meat()`](https://restatr.com/ivreg2r/reference/dot-cluster_meat.md).
When `small = TRUE`, applies the Stata `ivreg2` correction
`(N-1)/(N-K-sdofminus) * M/(M-1)`.

## Usage

``` r
.compute_cl_vcov(
  bread,
  X_hat,
  resid,
  cluster_vec,
  N,
  K,
  M,
  small,
  dofminus = 0L,
  sdofminus = 0L,
  weights = NULL,
  center = FALSE,
  weight_type = "aweight"
)
```

## Arguments

- bread:

  K x K bread matrix: \\(X'P_Z X)^{-1}\\ for IV, \\(X'X)^{-1}\\ for OLS.

- X_hat:

  N x K matrix: projected regressors for IV, original X for OLS.

- resid:

  Length-N residual vector.

- cluster_vec:

  Length-N vector (one-way) or list of 2 vectors (two-way) of cluster
  identifiers.

- N:

  Integer: number of observations.

- K:

  Integer: number of regressors.

- M:

  Integer: number of clusters (min(M1, M2) for two-way).

- small:

  Logical: if `TRUE`, apply the finite-sample correction
  `(N-1)/(N-K-sdofminus) * M/(M-1)`.

- dofminus:

  Integer: large-sample DoF adjustment (default 0). Note: dofminus does
  NOT appear in cluster VCV scaling (Stata convention).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

## Value

K x K variance-covariance matrix.
