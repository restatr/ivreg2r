# Extract design matrices from an ivreg2 model

Returns the regressor matrix (X), instrument matrix (Z), or projected
regressors (X_hat = P_Z X), as used in estimation. For models estimated
with `partial`, these are the post-partialling matrices that
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html), and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) correspond to. The
returned matrices do not depend on whether the model was fitted with
`x = TRUE`: when the matrices were not stored, they are reconstructed
from the model frame and the same partialling projection is re-applied.

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

## See also

[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md)

Other ivreg2 methods:
[`coef.ivreg2()`](https://restatr.com/ivreg2r/reference/coef.ivreg2.md),
[`confint.ivreg2()`](https://restatr.com/ivreg2r/reference/confint.ivreg2.md),
[`diagnostics()`](https://restatr.com/ivreg2r/reference/diagnostics.md),
[`first_stage()`](https://restatr.com/ivreg2r/reference/first_stage.md),
[`fitted.ivreg2()`](https://restatr.com/ivreg2r/reference/fitted.ivreg2.md),
[`formula.ivreg2()`](https://restatr.com/ivreg2r/reference/formula.ivreg2.md),
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md),
[`nobs.ivreg2()`](https://restatr.com/ivreg2r/reference/nobs.ivreg2.md),
[`predict.ivreg2()`](https://restatr.com/ivreg2r/reference/predict.ivreg2.md),
[`print.ivreg2()`](https://restatr.com/ivreg2r/reference/print.ivreg2.md),
[`print.summary.ivreg2()`](https://restatr.com/ivreg2r/reference/print.summary.ivreg2.md),
[`residuals.ivreg2()`](https://restatr.com/ivreg2r/reference/residuals.ivreg2.md),
[`summary.ivreg2()`](https://restatr.com/ivreg2r/reference/summary.ivreg2.md),
[`terms.ivreg2()`](https://restatr.com/ivreg2r/reference/terms.ivreg2.md),
[`update.ivreg2()`](https://restatr.com/ivreg2r/reference/update.ivreg2.md),
[`vcov.ivreg2()`](https://restatr.com/ivreg2r/reference/vcov.ivreg2.md)
