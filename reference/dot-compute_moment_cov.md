# Compute L x L moment covariance matrix Omega

Computes the instrument-space score covariance for the Hansen J test.
This is a *new* computation in Z-space, not the K x K meat from
vcov-robust.R.

## Usage

``` r
.compute_moment_cov(
  Z,
  residuals,
  weights,
  cluster_vec,
  N,
  iid,
  dofminus = 0L,
  weight_type = "aweight",
  kernel = NULL,
  bw = NULL,
  time_index = NULL,
  center = FALSE,
  psd = NULL,
  vcov_type = "HAC",
  ZwZ = NULL,
  sw = FALSE,
  ivar_vec = NULL
)
```

## Arguments

- Z:

  N x L instrument matrix.

- residuals:

  N x 1 residual vector.

- weights:

  Normalized weights (sum to N), or NULL.

- cluster_vec:

  Cluster membership vector, or NULL.

- N:

  Number of observations.

- iid:

  Logical: TRUE when no cluster, no kernel, and `vcov = "iid"`.

- dofminus:

  Integer: large-sample DoF adjustment (default 0). HC path divides by
  `N - dofminus` (Stata livreg2.do line 326); cluster path divides by
  `N` (line 545, no dofminus adjustment).

- kernel:

  Canonical kernel name, or NULL for non-HAC.

- bw:

  Numeric bandwidth, or NULL.

- time_index:

  List from
  [`.build_time_index()`](https://restatr.com/ivreg2r/reference/dot-build_time_index.md),
  or NULL.

## Value

L x L symmetric matrix S (psd-corrected when `psd` is set).

## Details

Moment-condition covariance S for any VCE type (Stata's m_omega)

Single entry point for the matrix Stata posts as `e(S)`: the estimated
covariance of the orthogonality conditions, in Stata's normalization.
The iid branch computes `sigma^2 * (Z'WZ)/N` with
`sigma^2 = RSS/(N - dofminus)` (livreg2.do:197-198, 235); every other
VCE type delegates to
[`.compute_omega()`](https://restatr.com/ivreg2r/reference/dot-compute_omega.md),
whose parameters this function shares. This is the former body of the
`omega_fn` closure in
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md), extracted
so the same formula serves both the GMM/J machinery and the stored
`fit$S`.
