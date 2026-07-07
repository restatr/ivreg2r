# Update and re-fit an ivreg2 model

Updates the formula and/or other arguments of an `ivreg2` call and
(optionally) re-fits the model.

## Usage

``` r
# S3 method for class 'ivreg2'
update(object, formula., ..., evaluate = TRUE)
```

## Arguments

- object:

  An object of class `"ivreg2"`.

- formula.:

  A formula to update the model formula (see
  [update.formula](https://rdrr.io/r/stats/update.formula.html)).
  Multi-part formula updates are supported.

- ...:

  Additional arguments to update in the call (e.g., `vcov = "robust"`,
  `data = new_data`).

- evaluate:

  Logical: if `TRUE` (default), evaluate the updated call; if `FALSE`,
  return the unevaluated call.

## Value

If `evaluate = TRUE`, a new `ivreg2` object. If `evaluate = FALSE`, the
unevaluated call.
