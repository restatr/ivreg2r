# Stock-Wright LM S statistic

Weak-instrument-robust LM test of H0: all endogenous coefficients = 0
and orthogonality conditions are valid. This is the LM counterpart of
the Anderson-Rubin Wald test.

## Usage

``` r
.compute_stock_wright(
  Z,
  X,
  y,
  weights,
  cluster_vec,
  vcov_type,
  N,
  K1,
  L1,
  endo_names,
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

- weights:

  Normalized weights (sum to N), or NULL.

- cluster_vec:

  Cluster membership vector, or NULL.

- vcov_type:

  Character: `"iid"`, `"robust"`, `"HAC"`, `"AC"`, or `"CL"`.

- N:

  Number of observations.

- K1:

  Number of endogenous regressors.

- L1:

  Number of excluded instruments.

- endo_names:

  Character vector of endogenous variable names.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

## Value

Named list with `stat`, `p`, `df`.

## Details

Unlike Sargan/Hansen J, the S statistic uses residuals from regressing
`y` on the included exogenous regressors only (constraining endogenous
coefficients to zero). For IID, a homoskedastic omega is used
(`sigma^2 * Z'WZ / N`); for robust/cluster, the HC/CL omega is used.
