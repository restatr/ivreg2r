# Cholesky-first linear solver

Tries Cholesky factorization first (fast, stable for PD matrices), falls
back to QR if Cholesky fails. Matches Stata's `cholqrsolve()`.

## Usage

``` r
.chol_solve(A, b)
```

## Arguments

- A:

  Symmetric positive-definite matrix.

- b:

  Right-hand side vector or matrix.

## Value

Solution x such that A x = b.
