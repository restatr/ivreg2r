# Compute clustered meat matrix (one-way or two-way)

Shared helper for the `rowsum() -> crossprod()` pattern used in VCV
computation and all diagnostic score covariance computations.

## Usage

``` r
.cluster_meat(scores, cluster_vec)
```

## Arguments

- scores:

  N x P score matrix (P = K for VCV, P = L for diagnostics).

- cluster_vec:

  Length-N vector (one-way) or list of 2 vectors (two-way) of cluster
  identifiers.

## Value

P x P symmetric meat matrix (unscaled).

## Details

For one-way clustering, computes
`crossprod(rowsum(scores, cluster_vec))`. For two-way clustering
(Cameron-Gelbach-Miller 2006), computes
`meat1 + meat2 - meat_intersection` where the intersection uses
`interaction(cv1, cv2, drop = TRUE)`.
