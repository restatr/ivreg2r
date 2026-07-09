# Augment data with ivreg2 model predictions and residuals

Adds `.fitted` and `.resid` columns to the model frame (or user-supplied
data).

## Usage

``` r
# S3 method for class 'ivreg2'
augment(x, data = NULL, ...)
```

## Arguments

- x:

  An object of class `"ivreg2"`.

- data:

  A data frame to augment. If `NULL` (default), uses the stored model
  frame (`x$model`). An error is raised if `model = FALSE` was used at
  estimation time and `data` is not supplied.

- ...:

  Additional arguments (ignored).

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with all original data columns plus `.fitted` and `.resid`.

## See also

[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md)

Other broom methods:
[`glance.ivreg2()`](https://restatr.com/ivreg2r/reference/glance.ivreg2.md),
[`tidy.ivreg2()`](https://restatr.com/ivreg2r/reference/tidy.ivreg2.md)

## Examples

``` r
data(mroz)
mroz_work <- subset(mroz, inlf == 1)
fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
              data = mroz_work, vcov = "robust")
augment(fit) |> head()
#> # A tibble: 6 × 9
#>    lwage exper expersq  educ   age kidslt6 kidsge6 .fitted  .resid
#>    <dbl> <int>   <int> <int> <int>   <int>   <int>   <dbl>   <dbl>
#> 1 1.21      14     196    12    32       1       0   1.20   0.0107
#> 2 0.329      5      25    12    30       0       2   0.962 -0.634 
#> 3 1.51      15     225    12    35       1       3   1.22   0.297 
#> 4 0.0921     6      36    12    34       0       3   0.995 -0.903 
#> 5 1.52       7      49    14    31       1       2   1.22   0.305 
#> 6 1.56      33    1089    12    54       0       0   1.26   0.299 

# Residuals vs fitted
aug <- augment(fit)
plot(aug$.fitted, aug$.resid, xlab = "Fitted", ylab = "Residuals")
abline(h = 0, lty = 2)

```
