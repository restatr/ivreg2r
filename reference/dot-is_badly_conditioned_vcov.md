# Scale-invariant singularity check for covariance matrices

Cholesky factorization can succeed on matrices that are numerically
singular in floating-point arithmetic. To avoid treating coefficient
scale differences as singularity, this check first rescales to
correlation form and then uses the reciprocal condition estimate.

## Usage

``` r
.is_badly_conditioned_vcov(V)
```

## Arguments

- V:

  Symmetric positive definite covariance matrix.

## Value

Logical scalar: `TRUE` when `V` is numerically singular.
