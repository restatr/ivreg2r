# Compute Stock-Watson (2008) panel-robust VCV

Sandwich estimator using
[`.sw_meat()`](https://restatr.com/ivreg2r/reference/dot-sw_meat.md).
The SW meat has its own internal normalization (`/(N - N_panels)`), so
the sandwich multiplies by `N` (not `N/(N-dofminus)` as for HC).
Small-sample correction is the standard robust correction
`(N-dofminus)/(N-K-dofminus-sdofminus)`.

## Usage

``` r
.compute_sw_vcov(
  bread,
  X_hat,
  resid,
  ivar_vec,
  N,
  K,
  small = FALSE,
  dofminus = 0L,
  sdofminus = 0L,
  weights = NULL,
  weight_type = "aweight",
  center = FALSE
)
```

## Arguments

- bread:

  K x K bread matrix.

- X_hat:

  N x K projected regressors (or X for OLS).

- resid:

  N-vector of residuals.

- ivar_vec:

  N-vector of panel identifiers.

- N, K:

  Integer dimensions.

- small:

  Logical: apply finite-sample correction.

- dofminus, sdofminus:

  Integer DoF adjustments.

- weights:

  Normalized weights or NULL.

- weight_type:

  Character.

- center:

  Logical.

## Value

K x K variance-covariance matrix.
