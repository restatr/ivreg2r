# Duplicate-check keys for one formula part

Terms containing a time-series operator are keyed by their canonical
term label, not their underlying variable, so `unem` (endogenous) and
`l(unem, 1)` (instrument) coexist — matching Stata, where `unem` and
`L.unem` are distinct varlist entries. All other terms contribute their
variable names, as before.

## Usage

``` r
.dup_check_keys(formula, rhs)
```

## Arguments

- formula:

  A 3-part `Formula` object.

- rhs:

  Which part.

## Value

Character vector of keys.
