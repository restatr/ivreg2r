# Dispatcher for underidentification and weak identification tests

Computes Anderson LM + Cragg-Donald F (always), plus Kleibergen-Paap rk
LM and rk Wald F when VCE is robust or clustered.

## Usage

``` r
.compute_id_tests(
  X,
  Z,
  y,
  residuals,
  weights,
  cluster_vec,
  vcov_type,
  N,
  K,
  L,
  K1,
  L1,
  M = NULL,
  endo_names,
  excluded_names,
  has_intercept,
  dofminus = 0L,
  sdofminus = 0L,
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

- y:

  N x 1 response vector.

- residuals:

  N x 1 residual vector from 2SLS.

- weights:

  Normalized weights or NULL.

- cluster_vec:

  Cluster vector or NULL.

- vcov_type:

  Character: "iid", "robust", "HAC", "AC", or "CL".

- N, K, L, K1, L1:

  Integer dimensions.

- endo_names:

  Character vector of endogenous variable names.

- excluded_names:

  Character vector of excluded instrument names.

- has_intercept:

  Logical.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

## Value

List with underid, weak_id, weak_id_robust (or NULL).
