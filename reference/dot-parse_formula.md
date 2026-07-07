# Parse a three-part IV formula and build model matrices

Parse a three-part IV formula and build model matrices

## Usage

``` r
.parse_formula(
  formula,
  data,
  weights = NULL,
  subset = NULL,
  na.action = na.omit,
  tvar = NULL,
  ivar = NULL
)
```

## Arguments

- formula:

  A formula: `y ~ exog | endo | instruments` (3-part IV) or `y ~ exog`
  (1-part OLS).

- data:

  A data frame.

- weights:

  Optional weights expression (evaluated in `data`).

- subset:

  Optional subset expression (evaluated in `data`).

- na.action:

  Function for handling `NA`s (default `na.omit`).

- tvar, ivar:

  Optional time/panel variable names (single character), required when
  the formula contains time-series operators
  [`l()`](https://restatr.com/ivreg2r/reference/ts-operators.md)/[`d()`](https://restatr.com/ivreg2r/reference/ts-operators.md).

## Value

A named list; see Details.
