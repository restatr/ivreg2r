# Update and re-fit an ivreg2 model

Updates the formula and/or other arguments of an `ivreg2` call and
(optionally) re-fits the model.

## Usage

``` r
# S3 method for class 'ivreg2'
update(object, formula., ..., evaluate = TRUE)
```

## Arguments

- object:

  An object of class `"ivreg2"`.

- formula.:

  A formula to update the model formula (see
  [update.formula](https://rdrr.io/r/stats/update.formula.html)).
  Multi-part formula updates are supported.

- ...:

  Additional arguments to update in the call (e.g., `vcov = "robust"`,
  `data = new_data`).

- evaluate:

  Logical: if `TRUE` (default), evaluate the updated call; if `FALSE`,
  return the unevaluated call.

## Value

If `evaluate = TRUE`, a new `ivreg2` object. If `evaluate = FALSE`, the
unevaluated call.

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
[`model.matrix.ivreg2()`](https://restatr.com/ivreg2r/reference/model.matrix.ivreg2.md),
[`nobs.ivreg2()`](https://restatr.com/ivreg2r/reference/nobs.ivreg2.md),
[`predict.ivreg2()`](https://restatr.com/ivreg2r/reference/predict.ivreg2.md),
[`print.ivreg2()`](https://restatr.com/ivreg2r/reference/print.ivreg2.md),
[`print.summary.ivreg2()`](https://restatr.com/ivreg2r/reference/print.summary.ivreg2.md),
[`residuals.ivreg2()`](https://restatr.com/ivreg2r/reference/residuals.ivreg2.md),
[`summary.ivreg2()`](https://restatr.com/ivreg2r/reference/summary.ivreg2.md),
[`terms.ivreg2()`](https://restatr.com/ivreg2r/reference/terms.ivreg2.md),
[`vcov.ivreg2()`](https://restatr.com/ivreg2r/reference/vcov.ivreg2.md)
