# Reject offset() terms in any formula part

Stata's `ivreg2` has no offset concept, so there is no parity target to
implement against; a silently-ignored offset (the model frame never
calls [`model.offset()`](https://rdrr.io/r/stats/model.extract.html)) is
worse than an explicit rejection.
[`terms()`](https://rdrr.io/r/stats/terms.html) sets a non-`NULL`
`"offset"` attribute only for the `offset(...)` function form, so a
column merely named `offset` used as an ordinary variable is unaffected.

## Usage

``` r
.check_no_offset(formula, n_rhs)
```

## Arguments

- formula:

  A `Formula` object (already validated by
  [`.check_formula_parts()`](https://restatr.com/ivreg2r/reference/dot-check_formula_parts.md)).

- n_rhs:

  Integer: number of RHS parts (1 or 3).

## Value

Invisible `NULL`. Stops if any part contains an
[`offset()`](https://rdrr.io/r/stats/offset.html) term.
