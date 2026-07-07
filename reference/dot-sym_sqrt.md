# Symmetric matrix square root via eigendecomposition

Computes the symmetric square root of a PSD matrix via
eigendecomposition. Negative eigenvalues (numerical noise) are clamped
to zero with a warning.

## Usage

``` r
.sym_sqrt(A)
```

## Arguments

- A:

  Symmetric matrix.

## Value

Symmetric square root matrix.
