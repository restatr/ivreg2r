# Compute HC meat matrix with weight-type dispatch

For aweight/pweight: meat = (w \* X \* e)' (w \* X \* e) = X' diag(w^2
e^2) X. For fweight: meat = X' diag(w \* e^2) X (linear in weights, not
quadratic). For unweighted: meat = X' diag(e^2) X.

## Usage

``` r
.hc_meat(basis, resid, weights = NULL, weight_type = "aweight", center = FALSE)
```

## Arguments

- basis:

  N x K matrix (X_hat for IV, X for OLS).

- resid:

  N-vector of residuals.

- weights:

  Normalized weights or NULL.

- weight_type:

  Character: `"aweight"`, `"fweight"`, or `"pweight"`.

- center:

  Logical: if `TRUE`, subtract the mean of moment conditions before
  computing the outer product. The mean computation depends on weight
  type, matching Stata's `m_omega()` (livreg2.do lines 179-188).

## Value

K x K symmetric meat matrix.
