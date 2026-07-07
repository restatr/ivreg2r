# Compute HAC variance-covariance matrix

Sandwich estimator: `V = bread * meat * bread * N/(N - dofminus)`. This
matches the HC scaling pattern used in
[`.compute_hc_vcov()`](https://restatr.com/ivreg2r/reference/dot-compute_hc_vcov.md).
Note: Stata's m_omega divides by (N-dofminus) in Z-space, but our bread
is `(X'P_Z X)^{-1}` (unscaled), so we use the HC-style `N/(N-dofminus)`
multiplier instead. When `small = TRUE`, applies the additional
correction `(N-dofminus)/(N-K-dofminus-sdofminus)`.

## Usage

``` r
.compute_hac_vcov(
  bread,
  X_hat,
  resid,
  time_index,
  kernel,
  bw,
  N,
  K,
  dofminus = 0L,
  sdofminus = 0L,
  small = FALSE,
  weights = NULL,
  weight_type = "aweight",
  center = FALSE
)
```

## Arguments

- bread:

  K x K bread matrix.

- X_hat:

  N x K projected regressor matrix (sorted by time).

- resid:

  N-vector of residuals (sorted by time).

- time_index:

  List from
  [`.build_time_index()`](https://restatr.com/ivreg2r/reference/dot-build_time_index.md).

- kernel:

  Canonical kernel name.

- bw:

  Numeric bandwidth.

- N:

  Integer: number of observations.

- K:

  Integer: number of regressors.

- dofminus:

  Integer: large-sample DoF adjustment.

- sdofminus:

  Integer: small-sample DoF adjustment.

- small:

  Logical: if `TRUE`, apply the finite-sample correction
  `(N-dofminus)/(N-K-dofminus-sdofminus)`.

- weights:

  Normalized weights (sorted) or NULL.

- weight_type:

  Character: weight type.

## Value

K x K variance-covariance matrix.
