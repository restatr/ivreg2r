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
