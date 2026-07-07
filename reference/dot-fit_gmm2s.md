# Fit two-step efficient GMM

Computes two-step efficient GMM estimates using the inverse of the
moment covariance matrix (Omega) as the optimal weighting matrix.
Mirrors Stata's `s_egmm()` (ivreg2.ado:5384-5490).

## Usage

``` r
.fit_gmm2s(
  parsed,
  small = FALSE,
  dofminus = 0L,
  sdofminus = 0L,
  omega_fn,
  omega_rank_bound = Inf
)
```

## Arguments

- parsed:

  A `parsed_formula` object from
  [`.parse_formula()`](https://restatr.com/ivreg2r/reference/dot-parse_formula.md).

- small:

  Logical: if `TRUE`, use `N-K` denominator for sigma.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

- omega_fn:

  Function: closure that takes a residual vector and returns the L x L
  moment covariance matrix Omega. Built in
  [`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md).

## Value

A named list with: `coefficients`, `residuals`, `fitted.values`, `vcov`,
`sigma`, `df.residual`, `rank`, `r.squared`, `adj.r.squared`, `rss`,
`r2u`, `r2c`, `mss`, `bread`, `bread_gmm`, `X_hat`, `j_stat`, `j_df`,
`j_p`, `omega`, `method`.

## Details

Algorithm:

1.  First step: call
    [`.fit_2sls()`](https://restatr.com/ivreg2r/reference/dot-fit_2sls.md)
    to get IV coefficients and residuals.

2.  Compute Omega from step-1 residuals via `omega_fn`.

3.  Rank-check Omega; error if singular.

4.  Efficient GMM re-estimation using optimal weighting matrix W = Omega
    inverse.

5.  Compute J statistic, R-squared, sigma, etc.

The efficient GMM VCV is V = N \* (X'Z Omega_inv Z'X)\_inv — the
sandwich collapses, so no separate meat computation is needed.
