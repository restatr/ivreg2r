# Build time-index structure for HAC/AC estimation

Computes sort order, time span, minimum time gap, and gap count from
time and panel variables.

## Usage

``` r
.build_time_index(tvar_vec, ivar_vec = NULL)
```

## Arguments

- tvar_vec:

  Numeric vector of time values (aligned to model frame rows).

- ivar_vec:

  Factor/character/numeric panel identifier vector, or NULL for pure
  time series.

## Value

List with components:

- sort_order:

  Integer permutation that sorts data by (ivar, tvar).

- tvar_sorted:

  Time values in sorted order.

- ivar_sorted:

  Panel IDs in sorted order, or NULL.

- T_span:

  Integer: max(tvar) - min(tvar) + 1.

- tdelta:

  Numeric: minimum nonzero gap between consecutive sorted time values
  (within panels if panel data).

- n_gaps:

  Integer: number of missing time periods in the grid.

- panel_info:

  List of per-panel start/end indices in sorted data, or NULL for pure
  time series.
