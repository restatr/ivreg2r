# Anderson-Rubin test for weak-instrument-robust inference

Computes the Anderson-Rubin F and chi-squared statistics for testing H0:
all endogenous variable coefficients are zero in the structural
equation. The test is valid regardless of instrument strength.

## Usage

``` r
.compute_anderson_rubin(
  Z,
  X,
  y,
  weights,
  cluster_vec,
  vcov_type,
  N,
  K,
  L,
  K1,
  L1,
  M,
  endo_names,
  excluded_names,
  dofminus = 0L,
  sdofminus = 0L,
  weight_type = "aweight",
  kernel = NULL,
  bw = NULL,
  time_index = NULL,
  center = FALSE,
  sw = FALSE,
  ivar_vec = NULL,
  psd = NULL
)
```

## Arguments

- Z:

  N x L instrument matrix.

- X:

  N x K regressor matrix.

- y:

  N-vector of responses.

- weights:

  Normalized weights (sum to N), or NULL.

- cluster_vec:

  Cluster membership vector, or NULL.

- vcov_type:

  Character: `"iid"`, `"robust"`, `"HAC"`, `"AC"`, or `"CL"`.

- N, K, L, K1, L1:

  Integer dimensions.

- M:

  Number of clusters (or NULL).

- endo_names:

  Character vector of endogenous variable names.

- excluded_names:

  Character vector of excluded instrument names.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

## Value

Named list: `f_stat`, `f_p`, `f_df1`, `f_df2`, `chi2_stat`, `chi2_p`,
`chi2_df`.
