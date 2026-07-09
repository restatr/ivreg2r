# Generalized inverse via Gauss-Jordan elimination with partial pivoting

At each step the largest remaining diagonal element is selected as the
pivot (partial pivoting); rows and columns whose pivot falls below a
tolerance are zeroed out. For a rank-deficient matrix this yields a
generalized inverse that zeroes the collinear rows rather than the
Moore-Penrose pseudoinverse, matching the rank-deficiency handling of
Stata's `invsym()` and hence how Stata's `test` command forms Wald
statistics.

## Usage

``` r
.syminv_sweep(A)
```

## Arguments

- A:

  Symmetric n x n matrix (typically a VCV submatrix).

## Value

The generalized inverse matrix, or `NULL` if no pivots succeed.
