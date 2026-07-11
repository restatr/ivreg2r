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

- `note`:

  Character. A Stata-parity disclosure, or `NA` when none applies. A
  note appears when an explicit option the user set is silently not
  honored inside an automatically computed diagnostic, matching Stata's
  `ivreg2`: a non-Bartlett `kernel` on the identification and redundancy
  tests (Stata's `ranktest` hard-codes Bartlett); `kiefer` on the
  identification tests (the Kiefer VCE structure is dropped); `psd` on
  the identification and redundancy tests (`ranktest` never receives
  it); `center` on the endogeneity test; and `method = "cue"` (or
  `"liml"`, for orthogonality) on the endogeneity and orthogonality
  C-statistics, which come from a recursive re-estimation that does not
  use those estimators. When more than one disclosure applies to the
  same test, the note column joins the sentences with a space; each note
  is a complete sentence, so splitting on sentence boundaries recovers
  the individual disclosures programmatically.

Rows appear only for the tests actually computed on the model: an IV fit
reports underidentification, weak identification, overidentification,
and — for every endogenous regressor, unless narrowed to a subset via
`endog =` — the endogeneity test by default, while `orthog =` and
`redundant =` add the corresponding orthogonality and redundancy rows
only when the user requests them. The Stock-Yogo rows (keys of the form
`"sy_iv_size_10"`) report critical values, not test statistics, so their
`p_value` is always `NA`. An OLS fit (single-part formula) has no
diagnostics at all and returns a zero-row tibble with the columns above.

The `note` column is unique to this accessor:
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) keep
their numeric-only tidy semantics and carry no disclosure text, so
`diagnostics()` is the surface that exposes a note programmatically. The
same notes are printed by
[`summary()`](https://rdrr.io/r/base/summary.html) as footnotes below
the diagnostics block.

## Small-sample correction

the overidentification row (Sargan/Hansen J) and the endogeneity row
(C-statistic) are never small-sample corrected, regardless of the main
model's `small` argument. The first-stage F-statistics are always
small-sample corrected, regardless of `small`; they are not columns of
this tibble but are reached via
[`first_stage()`](https://restatr.com/ivreg2r/reference/first_stage.md).
All three match Stata's `ivreg2`. See
[ivreg2r-conventions](https://restatr.com/ivreg2r/reference/ivreg2r-conventions.md)
for the full statement of which corrections `small` controls.

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
#> # A tibble: 12 × 8
#>    test              test_name statistic    df   df2   p_value tested_vars note 
#>    <chr>             <chr>         <dbl> <int> <int>     <dbl> <chr>       <chr>
#>  1 underid           Kleiberg…     40.1      2    NA  1.98e- 9 NA          NA   
#>  2 weak_id           Cragg-Do…     19.7     NA    NA NA        NA          NA   
#>  3 weak_id_robust    Kleiberg…     20.4     NA    NA NA        NA          NA   
#>  4 sy_iv_size_10     Stock-Yo…     19.9     NA    NA NA        NA          NA   
#>  5 sy_iv_size_15     Stock-Yo…     11.6     NA    NA NA        NA          NA   
#>  6 sy_iv_size_20     Stock-Yo…      8.75    NA    NA NA        NA          NA   
#>  7 sy_iv_size_25     Stock-Yo…      7.25    NA    NA NA        NA          NA   
#>  8 overid            Hansen J       1.86     1    NA  1.72e- 1 NA          NA   
#>  9 anderson_rubin_f  Anderson…     27.9      2  3003  1.04e-12 NA          NA   
#> 10 anderson_rubin_c… Anderson…     55.8      2    NA  7.52e-13 NA          NA   
#> 11 stock_wright      Stock-Wr…     54.0      2    NA  1.92e-12 NA          NA   
#> 12 endogeneity       Endogene…     27.0      1    NA  2.04e- 7 educ        NA   
```
