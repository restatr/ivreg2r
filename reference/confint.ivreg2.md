# Confidence intervals for ivreg2 coefficients

Computes confidence intervals using the t distribution when
`small = TRUE` was used at estimation time, and the standard normal
otherwise.

## Usage

``` r
# S3 method for class 'ivreg2'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `"ivreg2"`.

- parm:

  A specification of which parameters to give intervals for, either a
  numeric vector of positions or a character vector of names. If
  missing, all parameters are included.

- level:

  The confidence level (default 0.95).

- ...:

  Additional arguments (ignored).

## Value

A matrix with columns for the lower and upper confidence limits.
