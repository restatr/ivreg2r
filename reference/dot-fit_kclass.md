# Fit k-class estimator (LIML, Fuller, user-supplied k)

Computes k-class IV estimates. When `method = "liml"`, the LIML
eigenvalue lambda is computed from the concentration matrices. When
`fuller > 0`, applies the Fuller (1977) bias correction
`k = lambda - fuller / (N - L)`. When `method = "kclass"`, uses the
user-supplied `kclass` value directly.

## Usage

``` r
.fit_kclass(
  parsed,
  method = "liml",
  kclass = NULL,
  fuller = 0,
  small = FALSE,
  dofminus = 0L,
  sdofminus = 0L
)
```

## Arguments

- parsed:

  A `parsed_formula` object from
  [`.parse_formula()`](https://restatr.com/ivreg2r/reference/dot-parse_formula.md).

- method:

  Character: `"liml"` or `"kclass"`.

- kclass:

  Numeric scalar: user-supplied k value (used only when
  `method = "kclass"`).

- fuller:

  Numeric scalar: Fuller modification parameter (default 0).

- small:

  Logical: if `TRUE`, use `N-K` denominator for sigma.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

## Value

A named list with `coefficients`, `residuals`, `fitted.values`, `vcov`,
`sigma`, `df.residual`, `rank`, `r.squared`, `adj.r.squared`, `rss`,
`r2u`, `r2c`, `mss`, `bread`, `bread_kclass`, `X_hat`, `proj_coef`,
`lambda`, `kclass_value`, `method`, `fuller_param`.

## Details

The 2SLS-style bread `(X_hat'WX_hat)^{-1}` is returned for diagnostics
and downstream HC/CL VCV computation (H2), alongside the k-class VCV
`sigma^2 * solve(XhXh)`.
