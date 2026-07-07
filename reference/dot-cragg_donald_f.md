# Cragg-Donald Wald F statistic (weak identification)

CD F = min(eval/(1-eval)) \* (N - L) / L1 where eval are squared
canonical correlations. This equals the minimum eigenvalue of the CD
matrix scaled by the appropriate degrees of freedom.

## Usage

``` r
.cragg_donald_f(cc_result, N, L, L1, dofminus = 0L, sdofminus = 0L)
```

## Arguments

- cc_result:

  Output from
  [`.canonical_correlations()`](https://restatr.com/ivreg2r/reference/dot-canonical_correlations.md).

- N:

  Number of observations.

- L:

  Total number of instruments (including exogenous + intercept).

- L1:

  Number of excluded instruments.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

## Value

Named list with stat, test_name.
