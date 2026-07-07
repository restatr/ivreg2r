# Compute residuals from W-weighted GMM first step

Uses a user-supplied weighting matrix W to estimate beta via GMM, then
returns the residuals. Used in Path 3 (wmatrix + gmm2s) to feed the
omega_fn closure with W-step residuals for the efficient second step.

## Usage

``` r
.wmatrix_first_step_resid(parsed, W)
```

## Arguments

- parsed:

  A `parsed_formula` object from
  [`.parse_formula()`](https://restatr.com/ivreg2r/reference/dot-parse_formula.md).

- W:

  L x L user-supplied weighting matrix.

## Value

Numeric vector of residuals from the W-weighted first step.
