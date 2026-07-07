# Fit Continuously Updated GMM Estimator (CUE)

Minimizes the CUE objective
`J(b) = N * gbar(b)' * Omega(b)^{-1} * gbar(b)` where both moment
conditions `gbar` and moment covariance `Omega` are re-evaluated at each
candidate beta. Uses BFGS + Nelder-Mead for optimization. Starting
values are 2SLS (IID) or two-step efficient GMM (non-IID), matching
Stata ivreg2.ado lines 5918-5948.

## Usage

``` r
.fit_cue(
  parsed,
  small = FALSE,
  dofminus = 0L,
  sdofminus = 0L,
  omega_fn,
  b0 = NULL,
  iid = TRUE,
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

- b0:

  Optional numeric vector: fixed parameter vector for J evaluation (no
  optimization). NULL for standard CUE optimization.

- iid:

  Logical: if `TRUE` (default), use 2SLS starting values. If `FALSE`
  (non-IID VCE: robust, cluster, HAC, AC), use two-step efficient GMM
  starting values, matching Stata's behavior (ivreg2.ado changelog line
  6596).

## Value

A named list with the same fields as
[`.fit_gmm2s()`](https://restatr.com/ivreg2r/reference/dot-fit_gmm2s.md)
plus `convergence` (integer) and `cue_message` (string).
`method = "cue"`.

## Details

When `b0` is supplied, evaluates the CUE objective at the fixed
parameter vector without optimization (Stata's `b0()` option).
