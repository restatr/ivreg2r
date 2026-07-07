# Kleibergen-Paap rk statistic (LM or Wald variant)

Computes the rk chi-squared statistic from the canonical correlation
intermediates and the Kronecker score covariance.

## Usage

``` r
.kp_rk_stat(cc_result, shat0, N, K1, L1)
```

## Arguments

- cc_result:

  Output from
  [`.canonical_correlations()`](https://restatr.com/ivreg2r/reference/dot-canonical_correlations.md).

- shat0:

  Score covariance from
  [`.kp_omega()`](https://restatr.com/ivreg2r/reference/dot-kp_omega.md).

- N:

  Number of observations.

- K1:

  Number of endogenous regressors.

- L1:

  Number of excluded instruments.

## Value

Named list with chi2, df.
