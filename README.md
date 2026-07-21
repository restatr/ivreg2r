
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ivreg2r <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/restatr/ivreg2r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/restatr/ivreg2r/actions/workflows/R-CMD-check.yaml) <!-- badges: end -->

Instrumental variables and GMM estimation for R with built-in diagnostics, designed to reproduce Stata’s [`ivreg2`](https://ideas.repec.org/c/boc/bocode/s425401.html) output (Baum, Schaffer & Stillman, 2003, 2007).

## Installation

``` r
install.packages("ivreg2r")
```

The development version, containing any changes made since the last CRAN release, is on GitHub:

``` r
# install.packages("devtools")
devtools::install_github("restatr/ivreg2r")
```

## Example

``` r
library(ivreg2r)
data(card)

fit <- ivreg2(lwage ~ exper + expersq + black + smsa + south + smsa66 +
                reg662 + reg663 + reg664 + reg665 + reg666 + reg667 +
                reg668 + reg669 | educ | nearc2 + nearc4,
              data = card, vcov = "robust", small = TRUE)
tidy(fit)
#> # A tibble: 16 × 7
#>    term        estimate std.error statistic  p.value conf.low conf.high
#>    <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#>  1 (Intercept)  3.24     0.884        3.66  2.56e- 4  1.50      4.97   
#>  2 educ         0.157    0.0526       2.99  2.83e- 3  0.0540    0.260  
#>  3 exper        0.119    0.0230       5.18  2.41e- 7  0.0738    0.164  
#>  4 expersq     -0.00236  0.000368    -6.40  1.83e-10 -0.00308  -0.00163
#>  5 black       -0.123    0.0516      -2.39  1.70e- 2 -0.225    -0.0220 
#>  6 smsa         0.101    0.0314       3.20  1.37e- 3  0.0391    0.162  
#>  7 south       -0.143    0.0303      -4.73  2.34e- 6 -0.203    -0.0838 
#>  8 smsa66       0.0151   0.0212       0.712 4.77e- 1 -0.0264    0.0566 
#>  9 reg662       0.103    0.0379       2.71  6.79e- 3  0.0284    0.177  
#> 10 reg663       0.150    0.0371       4.04  5.51e- 5  0.0771    0.223  
#> 11 reg664       0.0476   0.0457       1.04  2.98e- 1 -0.0420    0.137  
#> 12 reg665       0.154    0.0509       3.03  2.44e- 3  0.0546    0.254  
#> 13 reg666       0.173    0.0535       3.23  1.25e- 3  0.0680    0.278  
#> 14 reg667       0.142    0.0524       2.71  6.75e- 3  0.0393    0.245  
#> 15 reg668      -0.0951   0.0588      -1.62  1.06e- 1 -0.210     0.0201 
#> 16 reg669       0.103    0.0427       2.41  1.59e- 2  0.0193    0.187
```

## Formula syntax

Three-part formula: `y ~ exogenous | endogenous | excluded_instruments`. Each variable appears once; exogenous regressors are automatically included as instruments.

| Part    | Contents              | Example                                |
|---------|-----------------------|----------------------------------------|
| LHS     | Dependent variable    | `lwage`                                |
| 1st RHS | Exogenous regressors  | `exper + expersq + black + smsa + ...` |
| 2nd RHS | Endogenous regressors | `educ`                                 |
| 3rd RHS | Excluded instruments  | `nearc2 + nearc4`                      |

A one-part formula (`y ~ x1 + x2`) estimates OLS.

## Features

**Estimators:** OLS, 2SLS, LIML, Fuller, k-class, two-step efficient GMM, CUE

**Standard errors:** classical, robust, one- and two-way clustering, HAC/AC (8 kernels), Kiefer, Driscoll-Kraay, cluster+kernel

**Diagnostics:** weak identification (Kleibergen-Paap rk, Cragg-Donald, Stock-Yogo critical values), underidentification, overidentification (Sargan, Hansen J), endogeneity, first-stage F/partial R-squared, Anderson-Rubin, Stock-Wright S, orthogonality, redundancy

**Integration:** `tidy()`, `glance()`, `augment()` via [broom](https://broom.tidymodels.org/)

## Full output

`summary()` prints the complete `ivreg2` output, including the automatic diagnostics described above, for the model fit in the example:

``` r
summary(fit)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq + black + smsa + south + 
#>     smsa66 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + 
#>     reg668 + reg669 | educ | nearc2 + nearc4, data = card, vcov = "robust", 
#>     small = TRUE)
#> 
#> Observations: 3,010 
#> VCV type:     Robust, small-sample corrected 
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  3.2367114  0.8842788   3.660 0.000256 ***
#> educ         0.1570593  0.0525525   2.989 0.002825 ** 
#> exper        0.1188149  0.0229516   5.177 2.41e-07 ***
#> expersq     -0.0023565  0.0003684  -6.397 1.83e-10 ***
#> black       -0.1232778  0.0516278  -2.388 0.017010 *  
#> smsa         0.1007530  0.0314458   3.204 0.001369 ** 
#> south       -0.1431945  0.0302678  -4.731 2.34e-06 ***
#> smsa66       0.0150626  0.0211683   0.712 0.476792    
#> reg662       0.1027473  0.0379324   2.709 0.006793 ** 
#> reg663       0.1499316  0.0371230   4.039 5.51e-05 ***
#> reg664       0.0475676  0.0456953   1.041 0.297972    
#> reg665       0.1544801  0.0509365   3.033 0.002444 ** 
#> reg666       0.1729728  0.0535430   3.231 0.001249 ** 
#> reg667       0.1420355  0.0523988   2.711 0.006753 ** 
#> reg668      -0.0950611  0.0587578  -1.618 0.105801    
#> reg669       0.1029760  0.0426891   2.412 0.015915 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.1702 
#> Adj. R-squared: 0.1660 
#> F(15, 2994):     51.4 (p < 2.2e-16)
#> Root MSE:       0.4053 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(2) = 16.37 (p = 0.0003)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           7.89 
#>   Kleibergen-Paap rk Wald F:     8.32 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       19.93 
#>      15%  maximal IV size       11.59 
#>      20%  maximal IV size       8.75 
#>      25%  maximal IV size       7.25 
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(1) = 1.27 (p = 0.2600)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(2,2993) = 5.28 (p = 0.0051)
#>   Anderson-Rubin Wald Chi-sq(2) = 10.63 (p = 0.0049)
#>   Stock-Wright LM S Chi-sq(2) = 10.61 (p = 0.0050)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 2.83 (p = 0.0925)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ                8.32    0.0002      0.0052    0.0052      8.32      8.32
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq, black, smsa, south, smsa66, reg662, reg663, reg664, reg665, reg666, reg667, reg668, reg669 
#> Excluded instruments:  nearc2, nearc4
```

## Validation

Every estimator, variance estimator, and diagnostic is tested against output generated by Stata’s `ivreg2` itself. Checked-in Stata scripts run each specification and write the results to machine-readable fixture files; the R test suite then reproduces them at tight, centrally documented tolerances. See [How ivreg2r is validated against Stata](https://restatr.com/ivreg2r/articles/validation.html) for more details.

## Stata migration

Common translations (see `vignette("introduction")` for the [full table](https://restatr.com/ivreg2r/articles/introduction.html#stata-migration-quick-reference) and [intentional differences](https://restatr.com/ivreg2r/articles/introduction.html#intentional-differences-from-stata)):

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

- `vignette("introduction")` — OLS, 2SLS, diagnostics, robust/cluster SEs, weights, broom integration, Stata migration reference
- `vignette("advanced-iv")` — LIML, Fuller, k-class, COVIV, two-way clustering, reduced-form, model comparison
- `vignette("time-series-gmm")` — HAC/AC, Kiefer, Driscoll-Kraay, GMM, CUE, partialling, user matrices

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

## Development

`ivreg2r` is part of [restatr](https://restatr.com/), a project that uses LLM assistance to port user-contributed Stata modules from the [SSC Archive](https://ideas.repec.org/s/boc/bocode.html) to R and Python. The package was developed with Claude Code, working from the original [GPL-3-licensed](https://www.gnu.org/licenses/gpl-3.0.txt) [`ivreg2` Stata source code](https://ideas.repec.org/c/boc/bocode/s425401.html), created by Baum, Schaffer, and Stillman. Every package feature has been numerically validated against the Stata original. For more details see [How ivreg2r is validated against Stata](https://restatr.com/ivreg2r/articles/validation.html).

## License

GPL-3. Adapted from the Stata package `ivreg2` by Christopher F Baum, Mark E Schaffer, and Steven Stillman, also distributed under GPL-3.
