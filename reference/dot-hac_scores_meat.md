# Compute HAC meat from pre-formed score vectors

Like
[`.hac_meat()`](https://restatr.com/ivreg2r/reference/dot-hac_meat.md)
but takes pre-formed N x P scores directly, without separate basis and
residual vectors. Used by the Kleibergen-Paap path where scores are
Kronecker products that cannot be factored into basis \* resid form.

## Usage

``` r
.hac_scores_meat(scores, time_index, kernel, bw, weights = NULL)
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

- weights:

  Normalized weights (sorted) or NULL. For the KP path, weights are
  already incorporated into the scores.

## Value

P x P symmetric meat matrix (unscaled).
