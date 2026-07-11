# Stock and Watson U.S. Macroeconomic Dataset

Quarterly U.S. macroeconomic time series (1959Q1–2000Q4) used in Baum,
Schaffer & Stillman (2007) to illustrate HAC standard errors, two-step
GMM, and CUE estimation. Original source: Stock and Watson (2003).

## Usage

``` r
stockwatson
```

## Format

A data frame with 168 observations and 17 variables:

- year:

  Calendar year.

- quarter:

  Quarter (1–4).

- date:

  Stata quarterly date (0 = 1960Q1).

- UR:

  Unemployment rate (percent).

- CPI:

  Consumer Price Index.

- FFIR:

  Federal funds interest rate.

- TBILL:

  Treasury bill rate.

- TBON:

  Treasury bond rate.

- ER:

  Trade-weighted exchange rate.

- GDP:

  Nominal GDP (billions of dollars).

- inf:

  Annualized CPI inflation: `100 * log(CPI / L4.CPI)`. `NA` for the
  first 4 observations.

- ggdp:

  Annualized GDP growth: `100 * log(GDP / L4.GDP)`. `NA` for the first 4
  observations.

- dinf:

  First difference of inflation (`inf - L.inf`). `NA` for the first 5
  observations.

- ggdp_2:

  Second lag of GDP growth (`L2.ggdp`).

- TBILL_1:

  First lag of Treasury bill rate (`L.TBILL`).

- ER_1:

  First lag of exchange rate (`L.ER`).

- TBON_1:

  First lag of Treasury bond rate (`L.TBON`).

## Source

Stock, J.H. and Watson, M.W. (2003). *Introduction to Econometrics*.
Addison-Wesley.

Downloaded from
<http://fmwww.bc.edu/ec-p/data/stockwatson/macrodat.dta>.

Redistribution basis: `system.file("COPYRIGHTS", package = "ivreg2r")`.

## Details

The derived variables `inf`, `ggdp`, `dinf`, and the lagged instrument
columns are pre-computed following Baum, Schaffer & Stillman (2007, p.
474) so that users can replicate their examples directly without
time-series operators.

## References

Baum, C.F., Schaffer, M.E. and Stillman, S. (2007). Enhanced routines
for instrumental variables/generalized method of moments estimation and
testing. *The Stata Journal*, 7(4), 465–506.

## See also

Other ivreg2r datasets:
[`abdata`](https://restatr.com/ivreg2r/reference/abdata.md),
[`card`](https://restatr.com/ivreg2r/reference/card.md),
[`cigar`](https://restatr.com/ivreg2r/reference/cigar.md),
[`griliches`](https://restatr.com/ivreg2r/reference/griliches.md),
[`grunfeld`](https://restatr.com/ivreg2r/reference/grunfeld.md),
[`klein`](https://restatr.com/ivreg2r/reference/klein.md),
[`mroz`](https://restatr.com/ivreg2r/reference/mroz.md),
[`nlswork`](https://restatr.com/ivreg2r/reference/nlswork.md),
[`phillips`](https://restatr.com/ivreg2r/reference/phillips.md),
[`wagepan`](https://restatr.com/ivreg2r/reference/wagepan.md)

## Examples

``` r
data(stockwatson)

# Baum, Schaffer & Stillman (2007), p. 475: 2SLS with IID errors
fit_iv <- ivreg2(dinf ~ 1 | UR | ggdp_2 + TBILL_1 + ER_1 + TBON_1,
                 data = stockwatson)
summary(fit_iv)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = dinf ~ 1 | UR | ggdp_2 + TBILL_1 + ER_1 + TBON_1, 
#>     data = stockwatson)
#> 
#> Observations: 158 
#> VCV type:     Classical (iid) 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept)  0.93807    0.29420   3.189  0.00143 **
#> UR          -0.15501    0.04833  -3.208  0.00134 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.1914 
#> Adj. R-squared: 0.1862 
#> Wald chi2(1):  10.2 (p = 0.0017)
#> Root MSE:       0.5543 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(4) = 58.66 (p = 0.0000)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           22.58 
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       24.58 
#>      15%  maximal IV size       13.96 
#>      20%  maximal IV size       10.26 
#>      25%  maximal IV size       8.31 
#>   Stock-Yogo critical values (IV relative bias):
#>      5%  maximal IV relative bias   16.85 
#>      10%  maximal IV relative bias   10.27 
#>      20%  maximal IV relative bias   6.71 
#>      30%  maximal IV relative bias   5.34 
#> 
#> Overidentification test (Sargan):
#>   Chi-sq(3) = 5.85 (p = 0.1191)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(4,153) = 3.44 (p = 0.0100)
#>   Anderson-Rubin Wald Chi-sq(4) = 14.23 (p = 0.0066)
#>   Stock-Wright LM S Chi-sq(4) = 13.05 (p = 0.0110)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 0.50 (p = 0.4783)
#>   Tested: UR
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   UR                 22.58    0.0000      0.3712    0.3712     22.58     22.58
#> 
#> Instrumented:          UR 
#> Excluded instruments:  ggdp_2, TBILL_1, ER_1, TBON_1 
#> 

# \donttest{
# Baum, Schaffer & Stillman (2007), p. 476: two-step GMM with HAC (Bartlett, bw=5)
fit_gmm <- ivreg2(dinf ~ 1 | UR | ggdp_2 + TBILL_1 + ER_1 + TBON_1,
                  data = stockwatson, method = "gmm2s",
                  vcov = "robust", kernel = "bartlett", bw = 5,
                  tvar = "date")
summary(fit_gmm)
#> 
#> 2-Step GMM Estimation
#> 
#> Call:
#> ivreg2(formula = dinf ~ 1 | UR | ggdp_2 + TBILL_1 + ER_1 + TBON_1, 
#>     data = stockwatson, vcov = "robust", method = "gmm2s", kernel = "bartlett", 
#>     bw = 5, tvar = "date")
#> 
#> Observations: 158 
#> VCV type:     HAC (kernel=Bartlett; bandwidth=5) 
#>               Estimates efficient for arbitrary autocorrelation
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  0.58508    0.37240   1.571    0.116
#> UR          -0.10024    0.06346  -1.580    0.114
#> ---
#> R-squared:      0.1548 
#> Adj. R-squared: 0.1493 
#> Wald chi2(1):  2.5 (p = 0.1185)
#> Root MSE:       0.5668 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(4) = 7.95 (p = 0.0933)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           22.58 
#>   Kleibergen-Paap rk Wald F:     7.36 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       24.58 
#>      15%  maximal IV size       13.96 
#>      20%  maximal IV size       10.26 
#>      25%  maximal IV size       8.31 
#>   Stock-Yogo critical values (IV relative bias):
#>      5%  maximal IV relative bias   16.85 
#>      10%  maximal IV relative bias   10.27 
#>      20%  maximal IV relative bias   6.71 
#>      30%  maximal IV relative bias   5.34 
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(3) = 3.57 (p = 0.3119)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(4,153) = 1.41 (p = 0.2343)
#>   Anderson-Rubin Wald Chi-sq(4) = 5.81 (p = 0.2137)
#>   Stock-Wright LM S Chi-sq(4) = 3.18 (p = 0.5278)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 1.03 (p = 0.3112)
#>   Tested: UR
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   UR                  7.36    0.0000      0.3712    0.3712      7.36      7.36
#> 
#> Instrumented:          UR 
#> Excluded instruments:  ggdp_2, TBILL_1, ER_1, TBON_1 
#> 
# }
```
