# Hansen J overidentification test (robust/cluster path)

Computes Omega from first-step residuals and delegates to
[`.compute_j_with_omega()`](https://restatr.com/ivreg2r/reference/dot-compute_j_with_omega.md)
for 2-step GMM re-estimation and J statistic.

## Usage

``` r
.hansen_j_test(
  Z,
  X,
  y,
  residuals,
  weights,
  cluster_vec,
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

  N x 1 first-step residual vector.

- weights:

  Normalized weights or NULL.

- cluster_vec:

  Cluster membership vector, or NULL.

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

Named list with `stat`, `p`, `df`, `test_name`.
