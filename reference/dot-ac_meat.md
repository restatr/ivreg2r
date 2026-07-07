# Compute AC meat matrix (autocorrelation consistent, homoskedastic)

Returns `shat / N` (the per-observation AC omega). The diagonal is
`sigma^2 * Z'WZ / N` where `sigma^2 = e'We / (N - dofminus)`. The VCV
wrapper multiplies by N: `V = N * bread * (shat/N) * bread`. Diagnostics
use `shat/N` directly as the omega (matching the IID convention
`sigma^2 * ZwZ / N`).

## Usage

``` r
.ac_meat(
  basis,
  resid,
  time_index,
  kernel,
  bw,
  N,
  dofminus,
  weights = NULL,
  weight_type = "aweight",
  ZwZ
)
```

## Arguments

- basis:

  N x P matrix (X_hat for VCV, Z for diagnostics), sorted.

- resid:

  N-vector of residuals, sorted.

- time_index:

  List from
  [`.build_time_index()`](https://restatr.com/ivreg2r/reference/dot-build_time_index.md).

- kernel:

  Canonical kernel name.

- bw:

  Numeric bandwidth.

- N:

  Integer: number of observations.

- dofminus:

  Integer: large-sample DoF adjustment.

- weights:

  Normalized weights (sorted) or NULL.

- weight_type:

  Character: weight type.

- ZwZ:

  P x P precomputed instrument cross-product `Z'WZ`.

## Value

P x P matrix `shat / N` (normalized).

## Details

The AC meat uses the Kronecker structure: `sigma_tau * Z'WZ_tau` at each
lag, where `sigma_tau` is a scalar autocovariance of the residuals and
`Z'WZ_tau` is the cross-product of instruments at lag tau.
