# Extract terms object from an ivreg2 model

Extract terms object from an ivreg2 model

## Usage

``` r
# S3 method for class 'ivreg2'
terms(x, component = c("regressors", "instruments", "full"), ...)
```

## Arguments

- x:

  An object of class `"ivreg2"`.

- component:

  Character: which terms object to return. `"regressors"` (default)
  returns terms for all regressors (exogenous + endogenous);
  `"instruments"` returns terms for excluded instruments (`NULL` for
  OLS); `"full"` returns terms for the complete formula.

- ...:

  Additional arguments (ignored).

## Value

A [terms](https://rdrr.io/r/stats/terms.html) object, or `NULL` if
`component = "instruments"` for an OLS model.
