# Assemble a sandwich VCV from the L x L moment covariance Omega

Stata assembles every plain (non-GMM) robust-family VCV from the L x L
instrument-space moment covariance returned by `m_omega`, by congruence
with the first-stage map A (`X_hat = Z A`):
`V = 1/N * aux5 * omega * aux5'` (s_iegmm, ivreg2.ado:5556-5558). In
unscaled R quantities that is `V = N * bread * (A' Omega A) * bread`,
followed by the small-sample factor. This route is used when `psd` is
set, so the psd correction is applied to Omega *before* the congruence —
`psd(A' S A) != A' psd(S) A`, which is why correcting the K x K meat (or
the final VCV) diverges from Stata whenever the correction binds.

## Usage

``` r
.vcov_from_omega(
  bread,
  A,
  omega,
  N,
  K,
  M = NULL,
  small = FALSE,
  dofminus = 0L,
  sdofminus = 0L,
  cluster = FALSE
)
```

## Arguments

- bread:

  K x K bread matrix (k-class bread for LIML/k-class).

- A:

  L x K first-stage coefficient matrix, or `NULL` for OLS (where Z = X
  and the congruence is the identity).

- omega:

  L x L moment covariance from
  [`.compute_moment_cov()`](https://restatr.com/ivreg2r/reference/dot-compute_moment_cov.md)
  (psd-corrected there when `psd` is set).

- N:

  Integer: number of observations.

- K:

  Integer: number of regressors.

- M:

  Integer: number of clusters (required when `cluster = TRUE`).

- small:

  Logical: apply the finite-sample correction.

- dofminus, sdofminus:

  Integer DoF adjustments.

- cluster:

  Logical: `TRUE` for the cluster-family small-sample factor
  `(N-1)/(N-K-sdofminus) * M/(M-1)`; `FALSE` for the robust-family
  factor `(N-dofminus)/(N-K-dofminus-sdofminus)`.

## Value

K x K variance-covariance matrix.

## Details

The single formula reproduces every plain path exactly because each VCE
type's Omega normalization cancels its sandwich scaling: HC and HAC use
`Omega = meat/(N - dofminus)` with
`V = bread meat bread * N/(N - dofminus)`; cluster, two-way, DK, and
cluster+kernel use `Omega = meat/N` with `V = bread meat bread`; AC and
Stock-Watson Omegas are already per-observation with
`V = bread Omega bread * N`.
