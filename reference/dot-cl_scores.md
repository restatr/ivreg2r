# Compute cluster scores (weight-type-agnostic), optionally centered

Cluster scores are `weights * basis * resid` for all weight types. The
definition of `weights` (normalized for aweight/pweight, raw for
fweight) makes this expression correct for all types.

## Usage

``` r
.cl_scores(
  basis,
  resid,
  weights = NULL,
  center = FALSE,
  weight_type = "aweight"
)
```

## Arguments

- basis:

  N x K matrix.

- resid:

  N-vector of residuals.

- weights:

  Normalized weights or NULL.

- center:

  Logical: if `TRUE`, center scores by subtracting their mean.

- weight_type:

  Character: `"aweight"`, `"fweight"`, or `"pweight"`.

## Value

N x K score matrix.

## Details

When `center = TRUE`, subtracts the (weighted) mean of scores before
returning, so all downstream cluster/kernel aggregation operates on
already-centered scores.
