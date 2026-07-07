# Endogeneity test (C-statistic / difference-in-J)

Tests H0: specified endogenous regressors are actually exogenous.
Computes a difference-of-J (robust/cluster) or difference-of-Sargan
(IID) statistic by re-estimating a restricted model where the tested
regressors are treated as exogenous (and hence become their own
instruments).

## Usage

``` r
.compute_endogeneity_test(
  Z,
  X,
  y,
  residuals,
  rss,
  weights,
  cluster_vec,
  vcov_type,
  N,
  K,
  L,
  K1,
  endo_names,
  endog_vars,
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

  N x L instrument matrix (full model).

- X:

  N x K regressor matrix.

- y:

  N x 1 response vector.

- residuals:

  N x 1 residual vector from the full model (unused directly, but
  included for API consistency).

- rss:

  Residual sum of squares from the full model (unused directly).

- weights:

  Normalized weights (sum to N), or NULL.

- cluster_vec:

  Cluster membership vector, or NULL.

- vcov_type:

  Character: "iid", "robust", "HAC", "AC", or "CL".

- N:

  Number of observations.

- K:

  Number of regressors.

- L:

  Number of instruments.

- K1:

  Number of endogenous regressors.

- endo_names:

  Character vector of endogenous regressor names.

- endog_vars:

  Character vector of variables to test (subset of `endo_names`), or
  NULL to test all.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

## Value

Named list with `stat`, `p`, `df`, `test_name`, `tested_vars`, or NULL
if this is not an IV model.

## Details

The same-S-matrix constraint guarantees C \>= 0: the S matrix comes from
the restricted model (larger instrument set), and the full model's J is
re-estimated using the appropriate submatrix of that S.

Unlike `.compute_orthog_test`, this test deliberately has NO `j_full`
parameter: Stata's endog test comes from a recursive call that
re-estimates BOTH the full and restricted models fresh (forwarding
`gmm2s`/`liml` but not `cue`, and passing no smatrix from the outer
model, ado:1576-1601), so the outer model's own J never enters
`e(estat)`; the internal fixed-omega computations here mirror that
exactly, and the sw_cue endog fixture verifies the CUE case to 8 digits.

## Note

The C-statistic for IID models is computed using large-sample
sigma-squared `e'e/(N-dofminus)` for both models. No small-sample
correction is applied even when `small = TRUE`. This matches Stata's
`ivreg2`.
