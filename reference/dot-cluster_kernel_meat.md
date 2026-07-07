# Compute cluster+kernel meat matrix (Driscoll-Kraay / Thompson time dimension)

Aggregates observation-level scores by time period, then accumulates
kernel-weighted cross-lag products of time-clustered scores. This is the
time-dimension component used by both Driscoll-Kraay (DK) and Thompson
(two-way cluster+kernel) VCE.

## Usage

``` r
.cluster_kernel_meat(
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

  Canonical kernel name.

- bw:

  Numeric bandwidth.

- weights:

  Normalized weights (sorted) or NULL.

- weight_type:

  Character: `"aweight"` or `"pweight"`.

## Value

P x P symmetric meat matrix (unscaled).

## Details

Algorithm matches livreg2.do lines 444–535.
