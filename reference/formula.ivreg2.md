# Extract formula from an ivreg2 object

Extract formula from an ivreg2 object

## Usage

``` r
# S3 method for class 'ivreg2'
formula(x, ...)
```

## Arguments

- x:

  An object of class `"ivreg2"`.

- ...:

  Additional arguments passed to
  [`Formula::formula.Formula()`](https://rdrr.io/pkg/Formula/man/Formula.html)
  (e.g., `rhs`, `lhs`, `collapse`).

## Value

The original model formula.
