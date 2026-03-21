
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ivreg2r

Instrumental variables and GMM estimation for R with built-in
diagnostics, designed to reproduce Stata’s
[`ivreg2`](https://ideas.repec.org/c/boc/bocode/s425401.html) output
(Baum, Schaffer & Stillman, 2003, 2007).

> **Note:** This package is not yet on CRAN. It is an independent
> reimplementation, not affiliated with or endorsed by the original
> `ivreg2` authors. Version 0.1.0. Requires R \>= 4.4.0.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("restatr/ivreg2r")
```

## Example

``` r
library(ivreg2r)
data(card)

fit <- ivreg2(lwage ~ exper + expersq + black + south + smsa + married |
                educ | nearc2 + nearc4,
              data = card, vcov = "robust", small = TRUE)
summary(fit)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq + black + south + smsa + 
#>     married | educ | nearc2 + nearc4, data = card, vcov = "robust", 
#>     small = TRUE)
#> 
#> Observations: 3,003 
#> VCV type:     Robust, small-sample corrected 
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  3.2920057  0.8423075   3.908 9.50e-05 ***
#> educ         0.1655671  0.0491900   3.366 0.000773 ***
#> exper        0.1125311  0.0226412   4.970 7.06e-07 ***
#> expersq     -0.0020232  0.0003786  -5.343 9.82e-08 ***
#> black       -0.0795648  0.0501282  -1.587 0.112567    
#> south       -0.0947853  0.0237991  -3.983 6.98e-05 ***
#> smsa         0.1223566  0.0312819   3.911 9.38e-05 ***
#> married     -0.0264721  0.0053246  -4.972 7.01e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.1440 
#> Adj. R-squared: 0.1420 
#> F(7, 2995):     116.8 (p < 2.2e-16)
#> Root MSE:       0.4110 
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           9.30 
#>   Kleibergen-Paap rk Wald F:     9.51 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       19.93 
#>      15%  maximal IV size       11.59 
#>      20%  maximal IV size       8.75 
#>      25%  maximal IV size       7.25 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(2) = 18.74 (p = 0.0001)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(2,2994) = 8.41 (p = 0.0002)
#>   Anderson-Rubin Wald Chi-sq(2) = 16.88 (p = 0.0002)
#>   Stock-Wright LM S Chi-sq(2) = 16.75 (p = 0.0002)
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(1) = 3.81 (p = 0.0510)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 4.45 (p = 0.0349)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ                9.51    0.0001      0.0062    0.0062      9.51      9.51
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq, black, south, smsa, married 
#> Excluded instruments:  nearc2, nearc4
```

Extract results programmatically with broom:

``` r
tidy(fit)
#> # A tibble: 8 × 7
#>   term        estimate std.error statistic      p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>        <dbl>    <dbl>     <dbl>
#> 1 (Intercept)  3.29     0.842         3.91 0.0000950     1.64      4.94   
#> 2 educ         0.166    0.0492        3.37 0.000773      0.0691    0.262  
#> 3 exper        0.113    0.0226        4.97 0.000000706   0.0681    0.157  
#> 4 expersq     -0.00202  0.000379     -5.34 0.0000000982 -0.00277  -0.00128
#> 5 black       -0.0796   0.0501       -1.59 0.113        -0.178     0.0187 
#> 6 south       -0.0948   0.0238       -3.98 0.0000698    -0.141    -0.0481 
#> 7 smsa         0.122    0.0313        3.91 0.0000938     0.0610    0.184  
#> 8 married     -0.0265   0.00532      -4.97 0.000000701  -0.0369   -0.0160
```

## Formula syntax

Three-part formula: `y ~ exogenous | endogenous | excluded_instruments`.
Each variable appears once; exogenous regressors are automatically
included as instruments.

| Part    | Contents              | Example                           |
|---------|-----------------------|-----------------------------------|
| LHS     | Dependent variable    | `lwage`                           |
| 1st RHS | Exogenous regressors  | `exper + expersq + black + south` |
| 2nd RHS | Endogenous regressors | `educ`                            |
| 3rd RHS | Excluded instruments  | `nearc2 + nearc4`                 |

A one-part formula (`y ~ x1 + x2`) estimates OLS.

## Features

**Estimators:** OLS, 2SLS, LIML, Fuller, k-class, two-step efficient
GMM, CUE

**Standard errors:** classical, robust, one- and two-way clustering,
HAC/AC (8 kernels), Kiefer, Driscoll-Kraay, cluster+kernel

**Diagnostics:** weak identification (Kleibergen-Paap rk, Cragg-Donald,
Stock-Yogo critical values), underidentification, overidentification
(Sargan, Hansen J), endogeneity, first-stage F/partial R-squared,
Anderson-Rubin, Stock-Wright S, orthogonality, redundancy

**Integration:** `tidy()`, `glance()`, `augment()` via
[broom](https://broom.tidymodels.org/)

All estimators and diagnostics are tested against Stata `ivreg2`
benchmarks.

## Stata migration

Common translations (see `vignette("introduction")` for the [full
table](https://restatr.github.io/ivreg2r/articles/introduction.html#stata-migration-quick-reference)
and [intentional
differences](https://restatr.github.io/ivreg2r/articles/introduction.html#intentional-differences-from-stata)):

| Stata | R (`ivreg2r`) |
|----|----|
| `ivreg2 y x (endo = z1 z2)` | `ivreg2(y ~ x \| endo \| z1 + z2, data = d)` |
| `ivreg2 ..., robust` | `ivreg2(..., vcov = "robust")` |
| `ivreg2 ..., robust small` | `ivreg2(..., vcov = "robust", small = TRUE)` |
| `ivreg2 ..., cluster(g)` | `ivreg2(..., clusters = ~ g)` |
| `ivreg2 ..., cluster(g) small` | `ivreg2(..., clusters = ~ g, small = TRUE)` |
| `ivreg2 ..., liml` | `ivreg2(..., method = "liml")` |
| `ivreg2 ..., gmm2s` | `ivreg2(..., method = "gmm2s")` |

## Learn more

- `vignette("introduction")` — OLS, 2SLS, diagnostics, robust/cluster
  SEs, weights, broom integration, Stata migration reference
- `vignette("advanced-iv")` — LIML, Fuller, k-class, COVIV, two-way
  clustering, reduced-form, model comparison
- `vignette("time-series-gmm")` — HAC/AC, Kiefer, Driscoll-Kraay, GMM,
  CUE, partialling, user matrices

## Citation

``` r
citation("ivreg2r")
DiTraglia F (2026). _ivreg2r: Extended Instrumental Variables
Estimation with Diagnostics_. R package,
<https://github.com/restatr/ivreg2r>.

The Stata ivreg2 command that inspired this package is described in:

  Baum C, Schaffer M, Stillman S (2003). "Instrumental Variables and
  GMM: Estimation and Testing." _The Stata Journal_, *3*(1), 1-31.
  doi:10.1177/1536867X0300300101
  <https://doi.org/10.1177/1536867X0300300101>.

Baum C, Schaffer M, Stillman S (2007). "Enhanced Routines for
Instrumental Variables/Generalized Method of Moments Estimation and
Testing." _The Stata Journal_, *7*(4), 465-506.
doi:10.1177/1536867X0800700402
<https://doi.org/10.1177/1536867X0800700402>.

To see these entries in BibTeX format, use 'print(<citation>,
bibtex=TRUE)', 'toBibtex(.)', or set
'options(citation.bibtex.max=999)'.
```

## License

MIT
