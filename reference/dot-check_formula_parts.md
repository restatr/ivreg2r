# Validate the number of RHS parts in a Formula object

Validate the number of RHS parts in a Formula object

## Usage

``` r
.check_formula_parts(formula)
```

## Arguments

- formula:

  A `Formula` object.

## Value

Integer: 1 (OLS) or 3 (IV). Stops on invalid counts.
