# Extract design matrices from an ivreg2 model

Returns the regressor matrix (X), instrument matrix (Z), or projected
regressors (X_hat = P_Z X). Matrices are retrieved from the stored `x`
component if available (when `ivreg2(..., x = TRUE)` was used),
otherwise reconstructed from the model frame.

## Usage

``` r
# S3 method for class 'ivreg2'
model.matrix(
  object,
  component = c("regressors", "projected", "instruments"),
  ...
)
```

## Arguments

- object:

  An object of class `"ivreg2"`.

- component:

  Character: which matrix to return. `"regressors"` (default) returns
  the regressor matrix X; `"instruments"` returns the full instrument
  matrix Z (`NULL` for OLS); `"projected"` returns the projected
  regressors X_hat (equals X for OLS).

- ...:

  Additional arguments (ignored).

## Value

A numeric matrix, or `NULL` if `component = "instruments"` for an OLS
model.

## Details

For models estimated with `partial`, the stored matrices (when
`x = TRUE`) are the post-partialling matrices. Reconstruction from the
model frame returns pre-partialling matrices.
