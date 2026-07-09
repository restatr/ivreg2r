# Wooldridge Wage Panel Dataset

Balanced panel data from the National Longitudinal Survey of Youth
(NLSY), 1980–1987. Contains 4,360 observations on 545 young men observed
over 8 years. Useful for demonstrating panel methods, including two-way
clustering by individual and year.

## Usage

``` r
wagepan
```

## Format

A data frame with 4,360 observations and 11 variables:

- nr:

  Person identifier.

- year:

  Calendar year (1980–1987).

- lwage:

  Log hourly wage.

- educ:

  Years of education.

- black:

  Black (binary).

- hisp:

  Hispanic (binary).

- exper:

  Years of labor market experience.

- expersq:

  Experience squared (`exper^2`).

- married:

  Married (binary).

- union:

  Union member (binary).

- hours:

  Annual hours worked.

## Source

Wooldridge, J.M. (2010). *Econometric Analysis of Cross Section and
Panel Data*, 2nd ed. MIT Press.

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
[`stockwatson`](https://restatr.com/ivreg2r/reference/stockwatson.md)

## Examples

``` r
data(wagepan)
# OLS with two-way clustering by person and year
# The canonical two-way clustering example from the ivreg2 help file lives
# in ?nlswork; this block illustrates the same VCE on a different panel.
fit <- ivreg2(lwage ~ educ + black + hisp + exper + expersq + married + union,
              data = wagepan, clusters = ~ nr + year, small = TRUE)
summary(fit)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ educ + black + hisp + exper + expersq + 
#>     married + union, data = wagepan, clusters = ~nr + year, small = TRUE)
#> 
#> Observations: 4,360 
#> VCV type:     Cluster-robust, small-sample corrected 
#> Clusters:     545 (nr), 8 (year)
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -0.0347056  0.1179535  -0.294 0.777115    
#> educ         0.0993878  0.0086400  11.503 8.44e-06 ***
#> black       -0.1438417  0.0511370  -2.813 0.026038 *  
#> hisp         0.0156980  0.0378971   0.414 0.691106    
#> exper        0.0891791  0.0150777   5.915 0.000591 ***
#> expersq     -0.0028487  0.0009645  -2.953 0.021301 *  
#> married      0.1076656  0.0234920   4.583 0.002535 ** 
#> union        0.1800725  0.0288443   6.243 0.000427 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.1866 
#> Adj. R-squared: 0.1853 
#> F(7, 7):     52.6 (p = 0.0000)
#> Root MSE:       0.4807 
#> 
```
