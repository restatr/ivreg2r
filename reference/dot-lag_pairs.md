# Find observation pairs separated by tau time periods

For each observation at time t, finds the observation at time t -
tau\*tdelta. Uses [`match()`](https://rdrr.io/r/base/match.html) for
O(N) lookup within panels.

## Usage

``` r
.lag_pairs(time_index, tau)
```

## Arguments

- time_index:

  List from
  [`.build_time_index()`](https://restatr.com/ivreg2r/reference/dot-build_time_index.md).

- tau:

  Integer lag (positive).

## Value

Two-column integer matrix `[i_now, i_lag]` with row indices into the
**sorted** data. Zero rows if no matches exist.
