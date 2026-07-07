# Instrument redundancy test (LM)

Computes the KP rk LM statistic testing H0: rank=0 for the first-stage
relationship between endogenous regressors and tested excluded
instruments, conditional on maintained instruments.

## Usage

``` r
.compute_redundancy_test(
  X,
  Z,
  weights,
  cluster_vec,
  vcov_type,
  N,
  K1,
  endo_colnames,
  excluded_colnames,
  redundant_vars,
  dofminus = 0L,
  weight_type = "aweight",
  kernel = NULL,
  bw = NULL,
  time_index = NULL,
  center = FALSE
)
```

## Arguments

- X:

  N x K regressor matrix.

- Z:

  N x L instrument matrix.

- weights:

  Normalized weights (sum to N), or NULL.

- cluster_vec:

  Cluster membership vector, or NULL.

- vcov_type:

  Character: "iid", "robust", "HAC", "AC", or "CL".

- N:

  Number of observations.

- K1:

  Number of endogenous regressors.

- endo_colnames:

  Character vector of endogenous regressor column names.

- excluded_colnames:

  Character vector of all excluded instrument column names.

- redundant_vars:

  Character vector of excluded instrument names to test.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- weight_type:

  Character: weight type.

- kernel:

  Canonical kernel name or NULL.

- bw:

  Numeric bandwidth or NULL.

- time_index:

  Time index list or NULL.

- center:

  Logical: center scores (default FALSE).

## Value

Named list with `stat`, `p`, `df`, `test_name`, `tested_vars`.
