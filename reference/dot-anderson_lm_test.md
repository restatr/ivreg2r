# Anderson canonical correlation LM test (IID underidentification)

Anderson canonical correlation LM test (IID underidentification)

## Usage

``` r
.anderson_lm_test(cc_result, N, L1, K1, dofminus = 0L)
```

## Arguments

- cc_result:

  Output from
  [`.canonical_correlations()`](https://restatr.com/ivreg2r/reference/dot-canonical_correlations.md).

- N:

  Number of observations.

- L1:

  Number of excluded instruments.

- K1:

  Number of endogenous regressors.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

## Value

Named list with stat, p, df, test_name.
