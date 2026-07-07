# FWL: partial out exogenous regressors from endogenous and excluded IVs

Returns the residuals from projecting X1 and Z1 onto X2, where X2 is the
exogenous regressor matrix (including intercept if present).

## Usage

``` r
.partial_fwl(X1, Z1, X2, weights)
```

## Arguments

- X1:

  N x K1 endogenous regressor matrix.

- Z1:

  N x L1 excluded instrument matrix.

- X2:

  N x K2 exogenous regressor matrix (including intercept).

- weights:

  Normalized weights (sum to N), or NULL.

## Value

List with `X1_perp` and `Z1_perp`.
