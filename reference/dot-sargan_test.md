# Sargan overidentification test (IID path)

Computes the Sargan statistic as a quadratic form in the
instrument-space score, normalized by the large-sample sigma^2 (rss/N,
always).

## Usage

``` r
.sargan_test(Z, residuals, rss, N, overid_df, weights, dofminus = 0L)
```

## Arguments

- Z:

  N x L instrument matrix.

- residuals:

  N x 1 residual vector.

- rss:

  Residual sum of squares (already weighted if applicable).

- N:

  Number of observations.

- overid_df:

  Degrees of freedom (L - K).

- weights:

  Normalized weights or NULL.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

## Value

Named list with `stat`, `p`, `df`, `test_name`.
