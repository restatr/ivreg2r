# Compute kernel weights for HAC estimation

Takes raw lag `tau` and bandwidth `bw`, normalizes internally as
`x = tau / bw`. Vectorized over `tau`.

## Usage

``` r
.kernel_weights(tau, bw, kernel)
```

## Arguments

- tau:

  Numeric vector of lags (non-negative integers in practice).

- bw:

  Numeric scalar bandwidth (\> 0).

- kernel:

  Canonical kernel name (from
  [`.validate_kernel()`](https://restatr.com/ivreg2r/reference/dot-validate_kernel.md)).

## Value

Numeric vector of kernel weights, same length as `tau`.

## Details

Formulas match Stata's `m_calckw()` (livreg2.do:627–664) with added
singularity guards for numerical stability:

- QS near zero: exponential interpolation (from sandwich::kweights)

- Daniell at zero: limit = 1

- Tent at zero: limit = tau^2 (and 0 when tau = 0)

- Lag-window kernels: clipped to 0 for x \> 1
