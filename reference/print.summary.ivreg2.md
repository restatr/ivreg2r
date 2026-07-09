# Print a summary.ivreg2 object

Formats and displays the full estimation output, including coefficient
table, fit statistics, and IV diagnostic tests.

## Usage

``` r
# S3 method for class 'summary.ivreg2'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars", TRUE),
  ...
)
```

## Arguments

- x:

  An object of class `"summary.ivreg2"`.

- digits:

  Minimum number of significant digits.

- signif.stars:

  Logical: print significance stars? Default `TRUE`.

- ...:

  Additional arguments passed to
  [`printCoefmat()`](https://rdrr.io/r/stats/printCoefmat.html).

## Value

`x`, invisibly.

## See also

[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md)

Other ivreg2 methods:
[`coef.ivreg2()`](https://restatr.com/ivreg2r/reference/coef.ivreg2.md),
[`confint.ivreg2()`](https://restatr.com/ivreg2r/reference/confint.ivreg2.md),
[`first_stage()`](https://restatr.com/ivreg2r/reference/first_stage.md),
[`fitted.ivreg2()`](https://restatr.com/ivreg2r/reference/fitted.ivreg2.md),
[`formula.ivreg2()`](https://restatr.com/ivreg2r/reference/formula.ivreg2.md),
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md),
[`model.matrix.ivreg2()`](https://restatr.com/ivreg2r/reference/model.matrix.ivreg2.md),
[`nobs.ivreg2()`](https://restatr.com/ivreg2r/reference/nobs.ivreg2.md),
[`predict.ivreg2()`](https://restatr.com/ivreg2r/reference/predict.ivreg2.md),
[`print.ivreg2()`](https://restatr.com/ivreg2r/reference/print.ivreg2.md),
[`residuals.ivreg2()`](https://restatr.com/ivreg2r/reference/residuals.ivreg2.md),
[`summary.ivreg2()`](https://restatr.com/ivreg2r/reference/summary.ivreg2.md),
[`terms.ivreg2()`](https://restatr.com/ivreg2r/reference/terms.ivreg2.md),
[`update.ivreg2()`](https://restatr.com/ivreg2r/reference/update.ivreg2.md),
[`vcov.ivreg2()`](https://restatr.com/ivreg2r/reference/vcov.ivreg2.md)
