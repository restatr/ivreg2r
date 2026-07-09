# Symmetric matrix square root via eigendecomposition

Computes the symmetric square root of a PSD matrix via
eigendecomposition. Eigenvalues that are negative only by floating-point
round-off are clamped to zero silently; a warning is issued only when a
negative eigenvalue is large enough to indicate a genuinely non-PSD
input.

## Usage

``` r
.sym_sqrt(A)
```

## Arguments

- A:

  Symmetric matrix.

## Value

Symmetric square root matrix.
