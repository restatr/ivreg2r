# Fit GMM with user-supplied weighting matrix (inefficient GMM)

Computes GMM estimates using user W as the weighting matrix. The VCV is
the full sandwich form (not the collapsed efficient form), and the J
statistic is computed via efficient 2-step re-estimation using the
estimated Omega (not the wmatrix-weighted residuals).

## Usage

``` r
.fit_gmm_wmatrix(
  parsed,
  small = FALSE,
  dofminus = 0L,
  sdofminus = 0L,
  W,
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

- W:

  L x L user-supplied weighting matrix.

- omega_fn:

  Function: closure that takes a residual vector and returns the L x L
  moment covariance matrix Omega.

## Value

A named list with the same fields as
[`.fit_gmm2s()`](https://restatr.com/ivreg2r/reference/dot-fit_gmm2s.md)
plus `method = "gmmw"`.

## Details

Mirrors Stata's `s_iegmm()` (ivreg2.ado:5495-5563).
