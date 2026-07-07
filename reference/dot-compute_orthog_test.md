# Instrument orthogonality test (C-statistic / difference-in-J)

Tests H0: specified instruments satisfy orthogonality conditions.
Computes a difference-of-J (robust/cluster) or difference-of-Sargan
(IID) statistic by comparing the full model to a restricted model where
the tested instruments are removed.

## Usage

``` r
.compute_orthog_test(
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
  orthog_vars,
  dofminus = 0L,
  weight_type = "aweight",
  kernel = NULL,
  bw = NULL,
  time_index = NULL,
  center = FALSE,
  psd = NULL,
  omega = NULL,
  j_full = NULL,
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

  N x 1 residual vector from the full model.

- rss:

  Residual sum of squares from the full model.

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

- orthog_vars:

  Character vector of instrument names to test.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- omega:

  Pre-computed Omega (S) matrix from the full model, or NULL. When
  non-NULL (GMM2S/CUE paths), used directly instead of recomputing from
  residuals. This ensures the C-statistic uses the same first-step Omega
  as the main model's J statistic (Hayashi 2000, p. 220).

- j_full:

  The full model's own overidentification statistic (Stata's
  `e(j)`/`e(sargan)`), or NULL. Stata computes the C-statistic as
  `cstat = j - cj` (ivreg2.ado:1547) where `j` is the main model's
  reported statistic, not a re-minimization of the fixed-S quadratic;
  the two coincide for 2SLS/GMM2S but differ for CUE (optimum of the
  self-consistent objective) and LIML (Sargan at the LIML estimate),
  since the recursive orthog call forwards `gmm2s` but not `cue`/`liml`
  (ado:1521-1538). When non-NULL and finite, used as J_full directly.

## Value

Named list with `stat`, `p`, `df`, `test_name`, `tested_vars`, or NULL
if this is not an IV model or orthog_vars is NULL.

## Details

This is the mirror image of the endogeneity test: the **full** model
(more instruments) provides the S (Omega) matrix, and the restricted
model (fewer instruments) is the constrained model. The same-S-matrix
constraint guarantees C \>= 0 (Hayashi 2000, p. 220).
