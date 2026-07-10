# Re-apply estimation's FWL partialling to reconstructed matrices

Reproduces the projection
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md) ran so
that matrices rebuilt from the model frame match the post-partialling
matrices stored when `x = TRUE`. The FWL residuals are invariant to a
positive rescaling of the weights, so the raw weights stored on the
object reproduce estimation's projection even though estimation used
normalized weights.

## Usage

``` r
.reapply_partial(object, X, Z)
```

## Arguments

- object:

  An `ivreg2` object recording a partial specification.

- X:

  Reconstructed pre-partialling regressor matrix (with an `assign`
  attribute and an `(Intercept)` column when the model has one).

- Z:

  Reconstructed pre-partialling instrument matrix, or `NULL`.

## Value

A list with post-partialling `X` and `Z`.
