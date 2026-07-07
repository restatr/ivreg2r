# Compute global mean of scores for SW centering

Helper that computes the centering term for SW meat when
`center = TRUE`. Matches the centering logic in
[`.hc_meat()`](https://restatr.com/ivreg2r/reference/dot-hc_meat.md) for
the corresponding weight type.

## Usage

``` r
.sw_eZmean(basis, resid, N, weights = NULL, weight_type = "aweight")
```

## Arguments

- basis:

  N x P matrix.

- resid:

  N-vector.

- N:

  Integer.

- weights:

  Normalized weights or NULL.

- weight_type:

  Character.

## Value

P-vector of score means.
