# Validate bandwidth for HAC estimation

Rules from `s_vkernel()` at livreg2.do:71–91 and ivreg2.ado:4746.

## Usage

``` r
.validate_bandwidth(bw, kernel)
```

## Arguments

- bw:

  Numeric \> 0, or the string `"auto"`.

- kernel:

  Canonical kernel name (from
  [`.validate_kernel()`](https://restatr.com/ivreg2r/reference/dot-validate_kernel.md)).

## Value

`bw` (invisibly), after validation. Issues warning for bw=1 with kernels
where this means zero lags.
