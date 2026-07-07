# Generalized inverse via Gauss-Jordan elimination with partial pivoting

Matches the behavior of Stata's `syminv()`: at each step, the largest
remaining diagonal element is selected as the pivot (partial pivoting).
Rows/columns whose pivot falls below a tolerance are zeroed out. This
differs from the Moore-Penrose pseudoinverse for rank-deficient matrices
and matches how Stata's `test` command computes Wald statistics.

## Usage

``` r
.syminv_sweep(A)
```

## Arguments

- A:

  Symmetric n x n matrix (typically a VCV submatrix).

## Value

The generalized inverse matrix, or `NULL` if no pivots succeed.
