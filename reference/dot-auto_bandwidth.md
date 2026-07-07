# Automatic bandwidth selection via Newey-West (1994)

Implements the Newey-West (1994, REStud 61(4):631-653) plug-in bandwidth
selector, matching Stata's `s_abw()` (ivreg2.ado:4796-4911). Supported
for Bartlett, Parzen, and Quadratic Spectral kernels only.

## Usage

``` r
.auto_bandwidth(resid, Z, time_index, kernel, has_intercept, N)
```

## Arguments

- resid:

  N-vector of residuals (sorted by time).

- Z:

  N x L instrument matrix (sorted by time), including exogenous
  regressors and (if present) the intercept column.

- time_index:

  List from
  [`.build_time_index()`](https://restatr.com/ivreg2r/reference/dot-build_time_index.md).

- kernel:

  Canonical kernel name (already validated as auto-compatible).

- has_intercept:

  Logical: whether the last column of Z is an intercept. If TRUE, the
  intercept column is excluded from the score process (Stata zeros it
  out via the `h` vector).

- N:

  Integer: number of observations (used as autocovariance denominator,
  matching Stata's `nobs`).

## Value

Numeric scalar: selected bandwidth (\>= 1). Integer for Bartlett/Parzen,
possibly fractional for Quadratic Spectral.
