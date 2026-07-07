# Fit 2SLS via two calls to lm.fit / lm.wfit

Stage 1 projects regressors X onto the instrument space Z. Stage 2
regresses y on the projected regressors X_hat. Residuals and fitted
values use the *original* X, not X_hat.

## Usage

``` r
.fit_2sls(parsed, small = FALSE, dofminus = 0L, sdofminus = 0L)
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
`bread`, `X_hat`, `proj_coef` (the L x K first-stage coefficient matrix
A with `X_hat = Z %*% A`, used for psd-corrected VCV assembly).

## Details

When weights are present, uses
[`lm.wfit()`](https://rdrr.io/r/stats/lmfit.html) for both stages.
`lm.wfit` returns fitted values on the original (unweighted) scale and a
QR of the weighted design matrix, so the bread is `(X_hat'WX_hat)^{-1}`.
