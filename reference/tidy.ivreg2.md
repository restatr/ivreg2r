# Tidy an ivreg2 object

Constructs a tibble summarizing coefficient estimates, standard errors,
test statistics, and p-values.

## Usage

``` r
# S3 method for class 'ivreg2'
tidy(x, conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE, ...)
```

## Arguments

- x:

  An object of class `"ivreg2"`.

- conf.int:

  Logical: include confidence intervals? Default `TRUE`.

- conf.level:

  Confidence level for intervals. Default `0.95`.

- exponentiate:

  Logical: exponentiate the coefficient estimates and confidence
  interval bounds? Default `FALSE`. Useful for log-linear models where
  `exp(estimate)` gives a multiplicative effect. Standard errors remain
  on the original (log) scale, following broom convention.

- ...:

  Additional arguments (ignored).

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with columns `term`, `estimate`, `std.error`, `statistic`, `p.value`,
and optionally `conf.low`, `conf.high`.

## Examples

``` r
data(mroz)
mroz_work <- subset(mroz, inlf == 1)
fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
              data = mroz_work, vcov = "robust")
tidy(fit)
#> # A tibble: 4 × 7
#>   term         estimate std.error statistic p.value conf.low conf.high
#>   <chr>           <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 (Intercept) -0.385     1.06        -0.363  0.717  -2.46    1.69     
#> 2 educ         0.0964    0.0865       1.11   0.265  -0.0731  0.266    
#> 3 exper        0.0422    0.0167       2.53   0.0113  0.00954 0.0748   
#> 4 expersq     -0.000832  0.000471    -1.77   0.0770 -0.00175 0.0000902
tidy(fit, conf.int = FALSE)
#> # A tibble: 4 × 5
#>   term         estimate std.error statistic p.value
#>   <chr>           <dbl>     <dbl>     <dbl>   <dbl>
#> 1 (Intercept) -0.385     1.06        -0.363  0.717 
#> 2 educ         0.0964    0.0865       1.11   0.265 
#> 3 exper        0.0422    0.0167       2.53   0.0113
#> 4 expersq     -0.000832  0.000471    -1.77   0.0770
tidy(fit, exponentiate = TRUE)
#> # A tibble: 4 × 7
#>   term        estimate std.error statistic p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.681  1.06        -0.363  0.717    0.0852      5.43
#> 2 educ           1.10   0.0865       1.11   0.265    0.930       1.30
#> 3 exper          1.04   0.0167       2.53   0.0113   1.01        1.08
#> 4 expersq        0.999  0.000471    -1.77   0.0770   0.998       1.00

# \donttest{
# Compare 2SLS and LIML side-by-side on the help-file baseline spec
# (weak first stage; see the LIML example in ?ivreg2 for the framing)
fit_liml <- ivreg2(lwage ~ exper + expersq | educ |
                     age + kidslt6 + kidsge6,
                   data = mroz_work, method = "liml")
comparison <- rbind(
  cbind(method = "2SLS", tidy(fit)),
  cbind(method = "LIML", tidy(fit_liml))
)
comparison[comparison$term == "educ", c("method", "estimate", "std.error")]
#>   method   estimate  std.error
#> 2   2SLS 0.09640024 0.08646259
#> 6   LIML 0.09575813 0.08369059
# }
```
