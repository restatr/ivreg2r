# Compute reduced-form regression(s) for IV models

In `"rf"` mode, fits y ~ Z and returns coefficients, VCV, residuals,
fitted values, F-test of excluded instruments, and RMSE. In `"system"`
mode, additionally fits each endogenous X_j ~ Z and constructs the
cross-equation VCV.

## Usage

``` r
.compute_reduced_form(
  mode,
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
  depvar_name,
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

- mode:

  Character: `"rf"` or `"system"`.

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

- depvar_name:

  Character: name of the dependent variable.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

## Value

A list with `mode` and per-equation (rf) or multi-equation (system)
regression results.
