# Time-series operators for ivreg2 formulas

`l()` (lag) and `d()` (difference) may be used inside an
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md) formula to
create lags and differences of a variable on the fly, resolved against
the time variable `tvar` (and panel variable `ivar` for panel data).
They mirror Stata's `tsset`-aware `L.` and `D.` operators:

## Usage

``` r
l(x, k = 1)

d(x, k = 1)
```

## Arguments

- x:

  A numeric variable in the model data.

- k:

  Lag order(s) for `l()`; difference order(s) for `d()`. Non-negative
  integer (vector allowed, expanded to one term per element).

## Value

When used inside an
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md) formula: a
numeric vector aligned to the model data. Calling these functions
anywhere else is an error.

## Details

|  |  |  |
|----|----|----|
| R syntax | Meaning | Stata equivalent |
| `l(x)`, `l(x, 1)` | lag 1 | `L.x` |
| `l(x, 1:3)` | lags 1, 2, 3 (three terms) | `L(1/3).x` |
| `l(x, 0)` | `x` itself | `L0.x` |
| `d(x)`, `d(x, 1)` | first difference `x - l(x, 1)` | `D.x` |
| `d(x, 2)` | second-order difference `x - 2 l(x, 1) + l(x, 2)` | `D2.x` |

## Lag semantics

Lags are computed by time **value**, not row position: `l(x, k)` at time
`t` is the value of `x` at time `t - k` within the same panel unit. If
no observation exists at `t - k` (start of the series, or a gap in the
time variable), the lag is `NA` and the row is dropped from the
estimation sample — exactly Stata's behavior. Lags never cross
panel-unit boundaries. The data does not need to be sorted.

## Differences from fixest

Following Stata (not fixest), `d(x, k)` is the k-th **order** (iterated)
difference `(1-L)^k x`, matching Stata's `Dk.x`; fixest's `d(x, k)` is
the distance-k difference `x - l(x, k)`. The two agree for `k = 1`.

## Usage constraints

Operators require `tvar` (and `ivar` for panel data) to be passed to
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md); `l()` and
`d()` are not standalone lag utilities and signal an error when called
outside an [`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md)
formula (use, e.g.,
[`dplyr::lag()`](https://dplyr.tidyverse.org/reference/lead-lag.html)
with a consecutive-periods check, or `fixest`/`collapse`/ `plm`, for
general-purpose panel lagging). The time variable must be integer-valued
(greater than \\-2^{53}\\, at most \\2^{53}\\), like Stata's `tsset` —
lag arithmetic on fractional or astronomically large time values is not
exact in floating point. The lag/difference order `k` must be a constant
non-negative integer (or integer vector for ranges); leads (negative
`k`), nested operators (`l(d(x))`), and ranges on the left-hand side
(the model has a single dependent variable) are not supported. When
operators appear among the *regressors*, `predict(fit, newdata)`
recomputes them within `newdata`, so `newdata` must contain the history
rows needed for the lags; rows whose lags are missing get `NA`
predictions (again matching Stata). Operators confined to the instrument
part impose no requirement on `newdata`.

## Examples

``` r
data(phillips)
# Stata: ivreg2 cinf (unem = l(1/3).unem), bw(3)
fit <- ivreg2(cinf ~ 1 | unem | l(unem, 1:3), data = phillips,
              tvar = "year", vcov = "AC", kernel = "bartlett", bw = 3)
coef(fit)
#> (Intercept)        unem 
#>   2.7901216  -0.4755749 
```
