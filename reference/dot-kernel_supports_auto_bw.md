# Check whether a kernel supports automatic bandwidth selection

Only Bartlett, Parzen, and Quadratic Spectral support `bw = "auto"`.
Matches ivreg2.ado:4746.

## Usage

``` r
.kernel_supports_auto_bw(kernel)
```

## Arguments

- kernel:

  Canonical kernel name.

## Value

Logical scalar.
