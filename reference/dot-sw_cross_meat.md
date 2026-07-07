# Compute cross-equation Stock-Watson meat block

Used by the reduced-form system VCV to compute the (i,j) block of the
stacked meat matrix. Same panel-loop algorithm as
[`.sw_meat()`](https://restatr.com/ivreg2r/reference/dot-sw_meat.md) but
with cross-equation residuals e_i, e_j.

## Usage

``` r
.sw_cross_meat(
  basis,
  resid_i,
  resid_j,
  panels,
  N,
  weights = NULL,
  weight_type = "aweight",
  center = FALSE,
  eZmean_i = NULL,
  eZmean_j = NULL
)
```

## Arguments

- basis:

  N x P matrix.

- resid_i, resid_j:

  N-vectors of residuals from equations i and j.

- panels:

  List of index vectors (output of `split(seq_len(N), ivar_vec)`).

- N:

  Integer.

- weights:

  Normalized weights or NULL.

- weight_type:

  Character.

- center:

  Logical.

- eZmean_i, eZmean_j:

  P-vectors or NULL.

## Value

P x P matrix (already normalized by `1/(N - N_panels)`).
