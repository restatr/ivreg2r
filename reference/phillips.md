# Phillips Curve Macroeconomic Dataset

U.S. macroeconomic time series (1948–1996) with inflation and
unemployment rates plus lags and first differences. This is the same
dataset used in Stata's `ivreg2` help file for HAC and AC examples,
originally from Wooldridge (2019).

## Usage

``` r
phillips
```

## Format

A data frame with 49 observations and 11 variables:

- year:

  Calendar year (1948–1996).

- inf:

  Inflation rate (percent).

- unem:

  Unemployment rate (percent).

- cinf:

  Change in inflation (`inf - inf_1`).

- cunem:

  Change in unemployment (`unem - unem_1`).

- unem_1:

  Lagged unemployment (one year).

- inf_1:

  Lagged inflation (one year).

- unem_2:

  Lagged unemployment (two years).

- inf_2:

  Lagged inflation (two years).

- cinf_1:

  Lagged change in inflation.

- cunem_1:

  Lagged change in unemployment.

## Source

Wooldridge, J.M. (2020). *Introductory Econometrics: A Modern Approach*,
7th ed. Cengage Learning.

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
[`stockwatson`](https://restatr.com/ivreg2r/reference/stockwatson.md),
[`wagepan`](https://restatr.com/ivreg2r/reference/wagepan.md)

## Examples

``` r
data(phillips)
# OLS Phillips curve with AC standard errors: exact replication of the
# Stata ivreg2 help-file example at line 1501 (`ivreg2 cinf unem, bw(3)`;
# Bartlett is Stata's default kernel). Without vcov = "robust", a kernel
# gives AC (autocorrelation-consistent) inference, not HAC.
fit_ac <- ivreg2(cinf ~ unem, data = phillips,
                 kernel = "bartlett", bw = 3, tvar = "year")
summary(fit_ac)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = cinf ~ unem, data = phillips, kernel = "bartlett", 
#>     bw = 3, tvar = "year")
#> 
#> Observations: 48 
#> VCV type:     AC (kernel=Bartlett; bandwidth=3) 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept)   3.0306     1.2230   2.478  0.01321 * 
#> unem         -0.5426     0.2054  -2.642  0.00825 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.1078 
#> Adj. R-squared: 0.0884 
#> Wald chi2(1):  6.7 (p = 0.0129)
#> Root MSE:       2.3992 
#> 

# Adding vcov = "robust" gives HAC (Newey-West), replicating the
# help-file example at line 1511.
fit_hac <- ivreg2(cinf ~ unem, data = phillips, vcov = "robust",
                  kernel = "bartlett", bw = 3, tvar = "year",
                  small = TRUE)
summary(fit_hac)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = cinf ~ unem, data = phillips, vcov = "robust", 
#>     small = TRUE, kernel = "bartlett", bw = 3, tvar = "year")
#> 
#> Observations: 48 
#> VCV type:     HAC (kernel=Bartlett; bandwidth=3), small-sample corrected 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)  
#> (Intercept)   3.0306     1.3895   2.181   0.0343 *
#> unem         -0.5426     0.2215  -2.450   0.0182 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.1078 
#> Adj. R-squared: 0.0884 
#> F(1, 46):     6.0 (p = 0.0182)
#> Root MSE:       2.4508 
#> 
```
