# Compute Stock-Watson (2008) panel-robust meat matrix

Implements the bias-corrected per-panel estimator from Stock & Watson
(2008, Econometrica, eq. 6). Panels with T \<= 2 are skipped (SW is
undefined). The meat has its own internal normalization: divided by
`(N - N_panels)`, NOT by `N` or `(N - dofminus)`.

## Usage

``` r
.sw_meat(
  basis,
  resid,
  ivar_vec,
  N,
  weights = NULL,
  weight_type = "aweight",
  center = FALSE,
  eZmean = NULL
)
```

## Arguments

- basis:

  N x P matrix (Z for diagnostics, X_hat for VCV).

- resid:

  N-vector of residuals.

- ivar_vec:

  N-vector of panel identifiers.

- N:

  Integer: number of observations.

- weights:

  Normalized weights or NULL.

- weight_type:

  Character: `"aweight"` or `"pweight"`.

- center:

  Logical: if `TRUE`, subtract the global mean of moment conditions
  before computing the outer product.

- eZmean:

  P-vector: pre-computed global mean of scores (used when
  `center = TRUE`), or NULL.

## Value

P x P symmetric meat matrix (already normalized by `1/(N - N_panels)`),
or a zero matrix with a warning if all panels have T \<= 2.
