# Fit OLS via lm.fit / lm.wfit (QR-based least squares)

When weights are present, uses
[`lm.wfit()`](https://rdrr.io/r/stats/lmfit.html) which internally
transforms by `sqrt(w)` and returns unweighted residuals/fitted values
with a QR of the weighted design matrix. The bread is `(X'WX)^{-1}`
automatically.

## Usage

``` r
.fit_ols(parsed, small = FALSE, dofminus = 0L, sdofminus = 0L)
```

## Arguments

- parsed:

  A `parsed_formula` object from
  [`.parse_formula()`](https://restatr.com/ivreg2r/reference/dot-parse_formula.md).

- small:

  Logical: if `TRUE`, use `N-K` denominator for sigma; if `FALSE`, use
  `N`.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

## Value

A named list with: `coefficients`, `residuals`, `fitted.values`, `vcov`,
`sigma`, `df.residual`, `rank`, `r.squared`, `adj.r.squared`, `rss`,
`bread`.
