# Compute HAC meat matrix (heteroskedasticity and autocorrelation consistent)

Returns the **unscaled** P x P meat matrix. The caller divides by
`(N - dofminus)` to get the HAC omega (Stata livreg2.do line 326).

## Usage

``` r
.hac_meat(
  basis,
  resid,
  time_index,
  kernel,
  bw,
  weights = NULL,
  weight_type = "aweight",
  center = FALSE
)
```

## Arguments

- basis:

  N x P matrix (X_hat for VCV, Z for diagnostics). Must be in sorted
  time order (matching `time_index`).

- resid:

  N-vector of residuals, in sorted time order.

- time_index:

  List from
  [`.build_time_index()`](https://restatr.com/ivreg2r/reference/dot-build_time_index.md).

- kernel:

  Canonical kernel name (from
  [`.validate_kernel()`](https://restatr.com/ivreg2r/reference/dot-validate_kernel.md)).

- bw:

  Numeric bandwidth.

- weights:

  Normalized weights (sorted) or NULL.

- weight_type:

  Character: `"aweight"` or `"pweight"` (fweight blocked).

## Value

P x P symmetric meat matrix (unscaled).
