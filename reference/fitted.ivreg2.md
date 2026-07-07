# Extract fitted values from an ivreg2 object

Respects `na.action`: when the model was fit with `na.exclude`, `NA`s
are reinserted at the omitted row positions so the result aligns with
the original data frame (matching base R's `fitted.lm`).

## Usage

``` r
# S3 method for class 'ivreg2'
fitted(object, ...)
```

## Arguments

- object:

  An object of class `"ivreg2"`.

- ...:

  Additional arguments (ignored).

## Value

Numeric vector of fitted values.
