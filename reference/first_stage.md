# Extract first-stage regression objects

Returns a named list of `ivreg2_first_stage` objects, one per endogenous
variable. Each object supports
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html), and
[`glance()`](https://generics.r-lib.org/reference/glance.html).

## Usage

``` r
first_stage(x, ...)

# S3 method for class 'ivreg2'
first_stage(x, ...)
```

## Arguments

- x:

  A fitted model object.

- ...:

  Additional arguments (ignored).

## Value

A named list of `ivreg2_first_stage` objects.

## Examples

``` r
# Mirrors the Stata `ivreg2` help-file example at line 1325, which
# demonstrates the equivalence of the Kleibergen-Paap rk Wald F and the
# first-stage F with a single endogenous regressor (and, in passing,
# the `savefirst` option that `first_stage = TRUE` corresponds to).
data(mroz)
mroz_work <- subset(mroz, inlf == 1)
fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
              data = mroz_work, vcov = "robust", first_stage = TRUE)
fs <- first_stage(fit)
coef(fs$educ)
#>  (Intercept)        exper      expersq          age      kidslt6      kidsge6 
#> 13.121940433  0.050141381 -0.001765715 -0.012015042  0.725442450 -0.221944688 
summary(fs$educ)
#> 
#> First-stage regression: educ ~ instruments
#> Observations: 428
#> VCV type:     robust
#> 
#> Coefficients:
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) 13.121940   0.818270  16.036  < 2e-16 ***
#> exper        0.050141   0.045075   1.112  0.26660    
#> expersq     -0.001766   0.001469  -1.202  0.23012    
#> age         -0.012015   0.017442  -0.689  0.49129    
#> kidslt6      0.725442   0.278648   2.603  0.00955 ** 
#> kidsge6     -0.221945   0.092030  -2.412  0.01631 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Root MSE:       2.2586 
#> F(3, 422) = 5.02 (p = 0.0020)
#> Partial R2:     0.0299 
#> Shea PR2:       0.0299 
#> 
# KP rk Wald F = the robust first-stage F (single endogenous regressor)
c(kp_rk_wald_F = fit$diagnostics$weak_id_robust$stat,
  first_stage_F = glance(fs$educ)$f_stat)
#>  kp_rk_wald_F first_stage_F 
#>      5.021219      5.021219 
```
