# Extract IV diagnostic tests as a tibble

Flattens the IV specification tests stored on a fitted model into a
single tibble, one row per test, so they can be filtered, joined, or
printed like any other tidyverse output instead of navigated as a nested
list (`x$diagnostics$overid$stat` and friends).

## Usage

``` r
diagnostics(x, ...)

# S3 method for class 'ivreg2'
diagnostics(x, ...)
```

## Arguments

- x:

  A fitted model object.

- ...:

  Additional arguments (ignored).

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with one row per diagnostic test present on the model and the following
columns:

- `test`:

  Character. A stable machine-readable key, meant for
  `filter(test == ...)`-style selection (e.g. `"overid"`,
  `"weak_id_robust"`, `"endogeneity"`).

- `test_name`:

  Character. The display name printed by
  [`summary()`](https://rdrr.io/r/base/summary.html).

- `statistic`:

  Double. The test statistic, or — for the Stock-Yogo rows — the
  critical value rather than a statistic.

- `df`:

  Integer. Degrees of freedom, `NA` when not applicable.

- `df2`:

  Integer. Denominator degrees of freedom, `NA` unless the statistic is
  F-distributed.

- `p_value`:

  Double. `NA` when not defined (e.g. Stock-Yogo critical values, or an
  exactly identified overidentification row).

- `tested_vars`:

  Character. The variable(s) tested; populated only for the endogeneity,
  orthogonality, and redundancy rows.

Rows appear only for the tests actually computed on the model: an IV fit
reports underidentification, weak identification, and overidentification
by default, while `endog =`, `orthog =`, and `redundant =` add the
corresponding endogeneity, orthogonality, and redundancy rows. The
Stock-Yogo rows (keys of the form `"sy_iv_size_10"`) report critical
values, not test statistics, so their `p_value` is always `NA`. An OLS
fit (single-part formula) has no diagnostics at all and returns a
zero-row tibble with the columns above.

## See also

[`first_stage()`](https://restatr.com/ivreg2r/reference/first_stage.md),
[`glance.ivreg2()`](https://restatr.com/ivreg2r/reference/glance.ivreg2.md)

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
[`print.summary.ivreg2()`](https://restatr.com/ivreg2r/reference/print.summary.ivreg2.md),
[`residuals.ivreg2()`](https://restatr.com/ivreg2r/reference/residuals.ivreg2.md),
[`summary.ivreg2()`](https://restatr.com/ivreg2r/reference/summary.ivreg2.md),
[`terms.ivreg2()`](https://restatr.com/ivreg2r/reference/terms.ivreg2.md),
[`update.ivreg2()`](https://restatr.com/ivreg2r/reference/update.ivreg2.md),
[`vcov.ivreg2()`](https://restatr.com/ivreg2r/reference/vcov.ivreg2.md)

## Examples

``` r
data(card)
fit <- ivreg2(
  lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
  data = card, vcov = "robust"
)
diagnostics(fit)
#> # A tibble: 12 × 7
#>    test                test_name     statistic    df   df2   p_value tested_vars
#>    <chr>               <chr>             <dbl> <int> <int>     <dbl> <chr>      
#>  1 underid             Kleibergen-P…     40.1      2    NA  1.98e- 9 NA         
#>  2 weak_id             Cragg-Donald…     19.7     NA    NA NA        NA         
#>  3 weak_id_robust      Kleibergen-P…     20.4     NA    NA NA        NA         
#>  4 sy_iv_size_10       Stock-Yogo c…     19.9     NA    NA NA        NA         
#>  5 sy_iv_size_15       Stock-Yogo c…     11.6     NA    NA NA        NA         
#>  6 sy_iv_size_20       Stock-Yogo c…      8.75    NA    NA NA        NA         
#>  7 sy_iv_size_25       Stock-Yogo c…      7.25    NA    NA NA        NA         
#>  8 overid              Hansen J           1.86     1    NA  1.72e- 1 NA         
#>  9 anderson_rubin_f    Anderson-Rub…     27.9      2  3003  1.04e-12 NA         
#> 10 anderson_rubin_chi2 Anderson-Rub…     55.8      2    NA  7.52e-13 NA         
#> 11 stock_wright        Stock-Wright…     54.0      2    NA  1.92e-12 NA         
#> 12 endogeneity         Endogeneity       27.0      1    NA  2.04e- 7 educ       
```
