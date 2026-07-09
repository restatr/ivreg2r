# Confidence intervals for ivreg2 coefficients

Computes confidence intervals using the t distribution when
`small = TRUE` was used at estimation time, and the standard normal
otherwise.

## Usage

``` r
# S3 method for class 'ivreg2'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `"ivreg2"`.

- parm:

  A specification of which parameters to give intervals for, either a
  numeric vector of positions or a character vector of names. If
  missing, all parameters are included.

- level:

  The confidence level (default 0.95).

- ...:

  Additional arguments (ignored).

## Value

A matrix with columns for the lower and upper confidence limits.

## See also

[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md);
[ivreg2r-conventions](https://restatr.com/ivreg2r/reference/ivreg2r-conventions.md)
for when the t vs normal distribution is used.

Other ivreg2 methods:
[`coef.ivreg2()`](https://restatr.com/ivreg2r/reference/coef.ivreg2.md),
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
[`update.ivreg2()`](https://restatr.com/ivreg2r/reference/update.ivreg2.md),
[`vcov.ivreg2()`](https://restatr.com/ivreg2r/reference/vcov.ivreg2.md)
