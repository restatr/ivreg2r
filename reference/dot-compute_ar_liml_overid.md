# Anderson-Rubin LIML overidentification statistics

Computes the LR and linearized forms of the Anderson-Rubin
overidentification test for LIML estimation under IID errors.

## Usage

``` r
.compute_ar_liml_overid(lambda, N, overid_df, dofminus = 0L)
```

## Arguments

- lambda:

  Numeric: LIML eigenvalue from
  [`.fit_kclass()`](https://restatr.com/ivreg2r/reference/dot-fit_kclass.md).

- N:

  Integer: number of observations.

- overid_df:

  Integer: degree of overidentification (L - K).

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

## Value

Named list with `lr_stat`, `lr_p`, `lin_stat`, `lin_p`, `df`, or a
zero-stat placeholder when exactly identified (df == 0).

## Details

Both statistics are chi-squared(L - K) under the null that all
overidentifying restrictions are valid. The LR form uses \\(N -
\text{dofminus}) \ln(\lambda)\\ and the linearized form uses \\(N -
\text{dofminus})(\lambda - 1)\\.
