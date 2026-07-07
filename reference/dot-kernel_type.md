# Classify kernel as lag-window or spectral-window

Matches `m_omega()` classification at livreg2.do:161–167.

## Usage

``` r
.kernel_type(kernel)
```

## Arguments

- kernel:

  Canonical kernel name (from
  [`.validate_kernel()`](https://restatr.com/ivreg2r/reference/dot-validate_kernel.md)).

## Value

`"lag"` or `"spectral"`.
