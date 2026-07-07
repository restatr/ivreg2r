# Glance at an ivreg2 object

Returns a single-row tibble of model-level summary statistics and
diagnostic test results.

## Usage

``` r
# S3 method for class 'ivreg2'
glance(x, diagnostics = TRUE, ...)
```

## Arguments

- x:

  An object of class `"ivreg2"`.

- diagnostics:

  Logical: include IV diagnostic test columns? Default `TRUE`. Set to
  `FALSE` for a compact summary without test statistics. Follows the
  same convention as broom's `glance.ivreg()`.

- ...:

  Additional arguments (ignored).

## Value

A single-row
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html).

**Always present** (33 columns): `r.squared`, `adj.r.squared`, `sigma`,
`statistic`, `p.value`, `df`, `df.residual`, `nobs`, `vcov_type`,
`small`, `weight_type`, `method`, `lambda`, `kclass_value`,
`fuller_parameter`, `coviv`, `center`, `psd`, `kernel`, `bw`, `kiefer`,
`dkraay`, `n_clusters1`, `n_clusters2`, `cue_convergence`, `partial_ct`,
`yy`, `yyc`, `rankxx`, `rankzz`, `condxx`, `condzz`, `ll`.

**When `diagnostics = TRUE`** (default, 24 additional columns):
`weak_id_stat`, `weak_id_robust_stat`, `underid_stat`, `underid_p`,
`overid_stat`, `overid_p`, `ar_overid_lr_stat`, `ar_overid_lr_p`,
`ar_overid_lin_stat`, `ar_overid_lin_p`, `ar_overid_df`,
`endogeneity_stat`, `endogeneity_p`, `stock_wright_stat`,
`stock_wright_p`, `stock_wright_df`, `orthog_stat`, `orthog_p`,
`redundancy_stat`, `redundancy_p`, `rf_f_stat`, `rf_f_p`, `ccev_min`,
`cdev_min`.

## Details

[`glance()`](https://generics.r-lib.org/reference/glance.html) always
returns the same columns for a given value of the `diagnostics`
argument, using `NA` for metrics that do not apply to the fitted model.

The `small` column indicates whether finite-sample corrections were
applied. When `small = TRUE`, test statistics are F-distributed; when
`small = FALSE`, they are chi-squared.

When `diagnostics = TRUE` (default), the output includes IV
specification tests. Columns that are conditionally `NA`:

- All diagnostic columns are `NA` for OLS models (1-part formula).

- `overid_stat`, `overid_p`: also `NA` when exactly identified (number
  of excluded instruments equals number of endogenous regressors).

- `weak_id_robust_stat`: `NA` when `vcov = "iid"` (Cragg-Donald F is
  used instead of Kleibergen-Paap).

- `ar_overid_*`: only non-`NA` for `method = "liml"` with
  `vcov = "iid"`.

- `orthog_*`, `redundancy_*`: `NA` unless `orthog` or `redundant` was
  specified.

- `rf_f_*`: `NA` unless `reduced_form = "rf"`.

Set `diagnostics = FALSE` for a compact summary without test statistics.

## Examples

``` r
data(mroz)
mroz_work <- subset(mroz, inlf == 1)
fit <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
              data = mroz_work, vcov = "robust")

# Full output with diagnostics
glance(fit)
#> # A tibble: 1 × 58
#>   r.squared adj.r.squared sigma statistic  p.value    df df.residual  nobs
#>       <dbl>         <dbl> <dbl>     <dbl>    <dbl> <int>       <int> <int>
#> 1     0.156         0.150 0.664      6.02 0.000508     3         424   428
#> # ℹ 50 more variables: vcov_type <chr>, small <lgl>, weight_type <chr>,
#> #   method <chr>, lambda <dbl>, kclass_value <dbl>, fuller_parameter <dbl>,
#> #   coviv <lgl>, center <lgl>, psd <chr>, kernel <chr>, bw <dbl>, kiefer <lgl>,
#> #   dkraay <dbl>, sw <lgl>, n_clusters1 <int>, n_clusters2 <int>,
#> #   cue_convergence <int>, partial_ct <int>, yy <dbl>, yyc <dbl>, rankxx <int>,
#> #   rankzz <int>, condxx <dbl>, condzz <dbl>, ll <dbl>, weak_id_stat <dbl>,
#> #   weak_id_robust_stat <dbl>, underid_stat <dbl>, underid_p <dbl>, …

# Compact output without diagnostics
glance(fit, diagnostics = FALSE)
#> # A tibble: 1 × 34
#>   r.squared adj.r.squared sigma statistic  p.value    df df.residual  nobs
#>       <dbl>         <dbl> <dbl>     <dbl>    <dbl> <int>       <int> <int>
#> 1     0.156         0.150 0.664      6.02 0.000508     3         424   428
#> # ℹ 26 more variables: vcov_type <chr>, small <lgl>, weight_type <chr>,
#> #   method <chr>, lambda <dbl>, kclass_value <dbl>, fuller_parameter <dbl>,
#> #   coviv <lgl>, center <lgl>, psd <chr>, kernel <chr>, bw <dbl>, kiefer <lgl>,
#> #   dkraay <dbl>, sw <lgl>, n_clusters1 <int>, n_clusters2 <int>,
#> #   cue_convergence <int>, partial_ct <int>, yy <dbl>, yyc <dbl>, rankxx <int>,
#> #   rankzz <int>, condxx <dbl>, condzz <dbl>, ll <dbl>

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
# }
```
