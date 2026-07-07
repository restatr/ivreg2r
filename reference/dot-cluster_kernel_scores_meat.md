# Compute cluster+kernel meat from pre-formed score vectors

Like
[`.cluster_kernel_meat()`](https://restatr.com/ivreg2r/reference/dot-cluster_kernel_meat.md)
but takes pre-formed N x P scores directly. Used by the Kleibergen-Paap
path where scores are Kronecker products.

## Usage

``` r
.cluster_kernel_scores_meat(scores, time_index, kernel, bw)
```

## Arguments

- scores:

  N x P score matrix, in sorted time order.

- time_index:

  List from
  [`.build_time_index()`](https://restatr.com/ivreg2r/reference/dot-build_time_index.md).

- kernel:

  Canonical kernel name.

- bw:

  Numeric bandwidth.

## Value

P x P symmetric meat matrix (unscaled).
