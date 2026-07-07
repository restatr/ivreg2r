# Partial out exogenous regressors via FWL projection

Projects y, X, and Z onto the orthogonal complement of the specified
exogenous regressor columns (P). Uses three branches matching Stata:

## Usage

``` r
.partial_out(parsed, partial_colnames, partialcons)
```

## Arguments

- parsed:

  A `parsed_formula` object from
  [`.parse_formula()`](https://restatr.com/ivreg2r/reference/dot-parse_formula.md).

- partial_colnames:

  Character vector of X column names to partial out.

- partialcons:

  Logical: whether the constant is being partialled.

## Value

A modified copy of `parsed` with:

- `y`, `X`, `Z` projected onto M_P

- Partialled columns removed from X and Z

- `has_intercept` set to FALSE (intercept absorbed)

## Details

1.  **Variables + constant**: demean everything, then regress on P, take
    residuals.

2.  **Variables only (no constant)**: regress on P directly, take
    residuals.

3.  **Constant only**: demean (subtract weighted/unweighted column
    means).
