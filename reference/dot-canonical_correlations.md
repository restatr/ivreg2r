# SVD-based canonical correlations and intermediates

Computes the Cholesky-standardized reduced-form coefficient matrix theta
= R_zz pihat irQyy and its full SVD, returning all intermediates needed
for both Anderson/CD (IID) and Kleibergen-Paap (robust) tests.

## Usage

``` r
.canonical_correlations(X1_perp, Z1_perp, N, K1, L1, weights)
```

## Arguments

- X1_perp:

  N x K1 partialled endogenous regressors.

- Z1_perp:

  N x L1 partialled excluded instruments.

- N:

  Number of observations.

- K1:

  Number of endogenous regressors.

- L1:

  Number of excluded instruments.

- weights:

  Normalized weights or NULL.

## Value

List with theta, U, cc, V, eval, pihat, irQyy, irQzz, or NULL if
Cholesky fails.
