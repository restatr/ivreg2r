# Dispatcher for overidentification tests

Returns Sargan (IID) or Hansen J (robust/cluster), or a zero-stat
placeholder for exactly identified models. Returns NULL for OLS.

## Usage

``` r
.compute_overid_test(
  Z,
  X,
  y,
  residuals,
  rss,
  weights,
  cluster_vec,
  vcov_type,
  is_iv,
  N,
  K,
  L,
  overid_df,
  dofminus = 0L,
  weight_type = "aweight",
  kernel = NULL,
  bw = NULL,
  time_index = NULL,
  center = FALSE,
  psd = NULL,
  sw = FALSE,
  ivar_vec = NULL
)
```

## Arguments

- Z:

  N x L instrument matrix.

- X:

  N x K regressor matrix.

- y:

  N x 1 response vector.

- residuals:

  N x 1 residual vector.

- rss:

  Residual sum of squares.

- weights:

  Normalized weights or NULL.

- cluster_vec:

  Cluster membership vector, or NULL.

- vcov_type:

  Character: "iid", "robust", "HAC", "AC", or "CL".

- is_iv:

  Logical: TRUE if this is an IV model.

- N:

  Number of observations.

- K:

  Number of regressors.

- L:

  Number of instruments.

- overid_df:

  Degrees of freedom (L - K).

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

## Value

Named list with `stat`, `p`, `df`, `test_name`, or NULL.

## Note

The Sargan statistic is normalized by the large-sample sigma-squared
`e'e/(N-dofminus)`. No small-sample correction is applied even when
`small = TRUE`. This matches Stata's `ivreg2`.
