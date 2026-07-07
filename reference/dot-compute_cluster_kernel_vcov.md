# Compute cluster+kernel VCV (Driscoll-Kraay or Thompson)

Sandwich estimator with three internal paths:

- DK (one-way cluster+kernel on tvar): meat = cluster_kernel_meat

- Thompson (two-way): meat = cluster_meat(ivar) +
  cluster_kernel_meat(tvar) - hc_meat Small-sample correction:
  `(N-1)/(N-K-sdofminus) * M/(M-1)` when `small=TRUE`.

## Usage

``` r
.compute_cluster_kernel_vcov(
  bread,
  X_hat,
  resid,
  cluster_vec,
  time_index,
  kernel,
  bw,
  N,
  K,
  M,
  small,
  dofminus = 0L,
  sdofminus = 0L,
  weights = NULL,
  weight_type = "aweight",
  is_twoway = FALSE,
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

- cluster_vec:

  Cluster vector (DK: tvar vector; Thompson: list of 2).

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

- M:

  Integer: effective cluster count (min(M1, M2) for Thompson).

- small:

  Logical: apply small-sample correction.

- dofminus:

  Integer: large-sample DoF adjustment.

- sdofminus:

  Integer: small-sample DoF adjustment.

- weights:

  Normalized weights (sorted) or NULL.

- weight_type:

  Character: weight type.

- is_twoway:

  Logical: TRUE for Thompson, FALSE for DK.

## Value

K x K variance-covariance matrix.
