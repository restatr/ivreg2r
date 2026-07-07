# Validate and canonicalize a kernel name

Matches Stata's abbreviation table from `s_vkernel()`
(livreg2.do:96–115), case-insensitive. Accepts Stata's "Danielle"
spelling but returns the standard mathematical name "Daniell".

## Usage

``` r
.validate_kernel(kernel)
```

## Arguments

- kernel:

  Character string: kernel name or abbreviation.

## Value

Canonical kernel name string.
