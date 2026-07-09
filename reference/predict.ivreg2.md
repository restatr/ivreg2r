# Predict from an ivreg2 model

Predict from an ivreg2 model

## Usage

``` r
# S3 method for class 'ivreg2'
predict(object, newdata, se.fit = FALSE, na.action = stats::na.pass, ...)
```

## Arguments

- object:

  An object of class `"ivreg2"`.

- newdata:

  An optional data frame for prediction. If omitted, fitted values from
  the original data are returned. If the model uses time-series
  operators (see
  [`ts-operators`](https://restatr.com/ivreg2r/reference/ts-operators.md))
  *among the regressors*, `newdata` must contain the time variable (and
  panel variable, if any) plus the history rows needed to compute the
  lags; rows whose lags are missing within `newdata` get `NA`
  predictions (matching Stata). Operators confined to the instrument
  part impose no such requirement — prediction scores only the
  regressors.

- se.fit:

  Logical: if `TRUE`, return prediction standard errors alongside fitted
  values. Standard errors are computed as `sqrt(diag(X V X'))` where
  `V = vcov(object)`, so they reflect the VCE used at estimation time
  (IID, robust, cluster, HAC, etc.). Not available after partialling
  (matching Stata's `predict, stdp`).

- na.action:

  Function for handling `NA`s in `newdata`.

- ...:

  Additional arguments (ignored).

## Value

When `se.fit = FALSE` (default), a numeric vector of predicted values.
When `se.fit = TRUE`, a list with components `fit` (predicted values)
and `se.fit` (standard errors of prediction).

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
[`print.ivreg2()`](https://restatr.com/ivreg2r/reference/print.ivreg2.md),
[`print.summary.ivreg2()`](https://restatr.com/ivreg2r/reference/print.summary.ivreg2.md),
[`residuals.ivreg2()`](https://restatr.com/ivreg2r/reference/residuals.ivreg2.md),
[`summary.ivreg2()`](https://restatr.com/ivreg2r/reference/summary.ivreg2.md),
[`terms.ivreg2()`](https://restatr.com/ivreg2r/reference/terms.ivreg2.md),
[`update.ivreg2()`](https://restatr.com/ivreg2r/reference/update.ivreg2.md),
[`vcov.ivreg2()`](https://restatr.com/ivreg2r/reference/vcov.ivreg2.md)

## Examples

``` r
fit <- ivreg2(lwage ~ exper | educ | nearc4, data = card)

# Fitted values on the estimation sample
head(predict(fit))
#>        1        2        3        4        5        6 
#> 5.420130 5.946853 6.730347 5.796737 6.730347 5.834925 

# Predictions with standard errors for new data
nd <- data.frame(exper = c(5, 10), educ = c(12, 16), nearc4 = c(0, 1))
predict(fit, newdata = nd, se.fit = TRUE)
#> $fit
#>        1        2 
#> 5.499142 7.106954 
#> 
#> $se.fit
#> [1] 0.1005463 0.1114562
#> 
```
