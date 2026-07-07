# Compute AC variance-covariance matrix

Sandwich estimator: `V = N * bread * (shat/N) * bread` using the AC
(autocorrelation consistent) meat. The
[`.ac_meat()`](https://restatr.com/ivreg2r/reference/dot-ac_meat.md)
returns `shat/N` (per-observation omega), so we multiply by N to get the
full VCV. This is analogous to the HAC/HC pattern where
`V = bread * meat * bread * N/(N-dofminus)`.

## Usage

``` r
.compute_ac_vcov(
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
  weight_type = "aweight"
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

- weights:

  Normalized weights (sorted) or NULL.

- weight_type:

  Character: weight type.

## Value

K x K variance-covariance matrix.
