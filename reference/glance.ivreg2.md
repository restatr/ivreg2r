# Glance at an ivreg2 object

Returns a single-row tibble of model-level summary statistics and the
headline IV diagnostics. The column set is deliberately compact so that
table tools such as modelsummary render a sensible default.

## Usage

``` r
# S3 method for class 'ivreg2'
glance(x, diagnostics = TRUE, ...)
```

## Arguments

- x:

  An object of class `"ivreg2"`.

- diagnostics:

  Logical: include the headline IV diagnostic columns? Default `TRUE`.
  Set to `FALSE` for a goodness-of-fit summary without the test
  statistics. Follows the same convention as broom's `glance.ivreg()`.

- ...:

  Additional arguments (ignored).

## Value

A single-row
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html).

**Always present** (9 columns): `r.squared`, `adj.r.squared`, `sigma`,
`statistic` (model F or Wald chi-squared), `p.value`, `df` (model
numerator degrees of freedom), `df.residual`, `nobs`, `vcov_type`.

**When `diagnostics = TRUE`** (default, 6 additional columns): the
headline IV specification tests — `weak_id_stat` (Cragg-Donald Wald F),
`weak_id_robust_stat` (Kleibergen-Paap rk Wald F), `underid_stat` and
`underid_p` (underidentification), `overid_stat` and `overid_p`
(Sargan/Hansen J overidentification).

The remaining stored quantities and configuration flags are **not** in
[`glance()`](https://generics.r-lib.org/reference/glance.html) — this
keeps the goodness-of-fit block usable in rendered tables. They remain
available as named elements on the fitted object: the estimation
`method`, `lambda`/`kclass_value`/`fuller_parameter`, `coviv`, `center`,
`psd`, `kernel`/`bw`, `kiefer`, `dkraay`, `sw`, cluster counts,
`cue_convergence`, `partial_ct`, `small`, the cross-products `yy`/`yyc`,
ranks and condition numbers (`rank`, `rankzz`, `condxx`, `condzz`), the
log-likelihood `ll`, and the full diagnostic list `x$diagnostics`
(endogeneity, orthogonality, redundancy, Anderson-Rubin, Stock-Wright,
and the Cragg-Donald/Kleibergen-Paap eigenvalues).

## Details

[`glance()`](https://generics.r-lib.org/reference/glance.html) returns a
fixed set of columns for a given value of `diagnostics`, using `NA` for
metrics that do not apply to the fitted model. All diagnostic columns
are `NA` for OLS models (single-part formula). `overid_stat` and
`overid_p` are also `NA` when the model is exactly identified (the
number of excluded instruments equals the number of endogenous
regressors), and `weak_id_robust_stat` is `NA` under `vcov = "iid"` (the
Cragg-Donald F in `weak_id_stat` is reported instead of the
Kleibergen-Paap F).

Set `diagnostics = FALSE` for a compact goodness-of-fit summary without
the IV test columns.

## Examples

``` r
data(mroz)
mroz_work <- subset(mroz, inlf == 1)
fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
              data = mroz_work, vcov = "robust")

# Full output with diagnostics
glance(fit)
#> # A tibble: 1 × 15
#>   r.squared adj.r.squared sigma statistic  p.value    df df.residual  nobs
#>       <dbl>         <dbl> <dbl>     <dbl>    <dbl> <int>       <int> <int>
#> 1     0.156         0.150 0.664      6.02 0.000508     3         424   428
#> # ℹ 7 more variables: vcov_type <chr>, weak_id_stat <dbl>,
#> #   weak_id_robust_stat <dbl>, underid_stat <dbl>, underid_p <dbl>,
#> #   overid_stat <dbl>, overid_p <dbl>

# Compact output without diagnostics
glance(fit, diagnostics = FALSE)
#> # A tibble: 1 × 9
#>   r.squared adj.r.squared sigma statistic  p.value    df df.residual  nobs
#>       <dbl>         <dbl> <dbl>     <dbl>    <dbl> <int>       <int> <int>
#> 1     0.156         0.150 0.664      6.02 0.000508     3         424   428
#> # ℹ 1 more variable: vcov_type <chr>

# Extract specific diagnostics
glance(fit)[, c("overid_stat", "overid_p")]
#> # A tibble: 1 × 2
#>   overid_stat overid_p
#>         <dbl>    <dbl>
#> 1       0.514    0.773
glance(fit)[, c("weak_id_stat", "weak_id_robust_stat")]
#> # A tibble: 1 × 2
#>   weak_id_stat weak_id_robust_stat
#>          <dbl>               <dbl>
#> 1         4.34                5.02

# Diagnostics dropped from glance() remain on the fitted object
fit$diagnostics$endogeneity
#> $stat
#> [1] 0.001299792
#> 
#> $p
#> [1] 0.9712404
#> 
#> $df
#> [1] 1
#> 
#> $test_name
#> [1] "Endogeneity"
#> 
#> $tested_vars
#> [1] "educ"
#> 

# \donttest{
# Compare Sargan (IID) vs Hansen J (robust)
fit_iid <- ivreg2(lwage ~ exper + expersq | educ |
                    age + kidslt6 + kidsge6, data = mroz_work)
data.frame(
  vcov = c("iid", "robust"),
  overid = c(glance(fit_iid)$overid_stat, glance(fit)$overid_stat),
  overid_p = c(glance(fit_iid)$overid_p, glance(fit)$overid_p)
)
#>     vcov    overid  overid_p
#> 1    iid 0.7015124 0.7041554
#> 2 robust 0.5138488 0.7734267

# A compact modelsummary table built from the curated glance() columns
if (requireNamespace("modelsummary", quietly = TRUE)) {
  modelsummary::modelsummary(
    list("2SLS" = fit),
    statistic = "std.error",
    gof_map = c("nobs", "r.squared", "weak_id_stat", "overid_stat")
  )
}
#> 
#> +--------------+---------+
#> |              | 2SLS    |
#> +==============+=========+
#> | (Intercept)  | -0.385  |
#> +--------------+---------+
#> |              | (1.060) |
#> +--------------+---------+
#> | educ         | 0.096   |
#> +--------------+---------+
#> |              | (0.086) |
#> +--------------+---------+
#> | exper        | 0.042   |
#> +--------------+---------+
#> |              | (0.017) |
#> +--------------+---------+
#> | expersq      | -0.001  |
#> +--------------+---------+
#> |              | (0.000) |
#> +--------------+---------+
#> | Num.Obs.     | 428     |
#> +--------------+---------+
#> | R2           | 0.156   |
#> +--------------+---------+
#> | weak_id_stat | 4       |
#> +--------------+---------+
#> | overid_stat  | 1       |
#> +--------------+---------+ 
# }
```
