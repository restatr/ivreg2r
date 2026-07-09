# Griliches (1976) Young Men's Wages Dataset

Wages and characteristics of young men from the National Longitudinal
Survey (NLS), originally analyzed in Griliches (1976). This dataset is a
standard IV example used in Hayashi (2000, Ch. 3) and Baum, Schaffer &
Stillman (2007, pp. 492–494) to illustrate weak identification tests,
Anderson–Rubin inference, Stock–Wright S statistics, and instrument
redundancy.

## Usage

``` r
griliches
```

## Format

A data frame with 758 observations and 20 variables:

- rns:

  Residency in the South (1 = yes).

- rns80:

  Residency in the South in 1980.

- mrt:

  Marital status (1 = married).

- mrt80:

  Marital status in 1980.

- smsa:

  Resides in metropolitan area (1 = urban).

- smsa80:

  Resides in metropolitan area in 1980.

- med:

  Mother's education (years).

- iq:

  IQ score.

- kww:

  Score on Knowledge of the World of Work test.

- year:

  Survey year (66–73 or 80).

- age:

  Age at survey.

- age80:

  Age in 1980.

- s:

  Completed years of schooling.

- s80:

  Completed years of schooling in 1980.

- expr:

  Work experience (years).

- expr80:

  Work experience in 1980 (years).

- tenure:

  Job tenure (years).

- tenure80:

  Job tenure in 1980 (years).

- lw:

  Log wage.

- lw80:

  Log wage in 1980.

## Source

Griliches, Z. (1976). Wages of very young men. *Journal of Political
Economy*, 84(4), S69–S85.

Downloaded from <http://fmwww.bc.edu/ec-p/data/hayashi/griliches76.dta>.

## Details

Two cross-sections are pooled (survey years 66–73 and 80). The Baum,
Schaffer & Stillman (2007) examples use year dummies via `factor(year)`.

## References

Hayashi, F. (2000). *Econometrics*. Princeton University Press. (Chapter
3, p. 255.)

Baum, C.F., Schaffer, M.E. and Stillman, S. (2007). Enhanced routines
for instrumental variables/generalized method of moments estimation and
testing. *The Stata Journal*, 7(4), 465–506. (pp. 492–494.)

## Examples

``` r
data(griliches)

# Baum, Schaffer & Stillman (2007), p. 493: IV with robust SEs, weak-ID
# diagnostics, redundancy test
fit <- ivreg2(lw ~ s + expr + tenure + rns + smsa + factor(year) |
                iq | age + mrt,
              data = griliches, vcov = "robust",
              redundant = "mrt")
summary(fit)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lw ~ s + expr + tenure + rns + smsa + factor(year) | 
#>     iq | age + mrt, data = griliches, vcov = "robust", redundant = "mrt")
#> 
#> Observations: 758 
#> VCV type:     Robust 
#> 
#> Coefficients:
#>                 Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)    10.550965   2.781762   3.793 0.000149 ***
#> iq             -0.094890   0.041890  -2.265 0.023500 *  
#> s               0.339712   0.118327   2.871 0.004092 ** 
#> expr           -0.006604   0.029255  -0.226 0.821406    
#> tenure          0.084885   0.030668   2.768 0.005643 ** 
#> rns            -0.376939   0.155997  -2.416 0.015678 *  
#> smsa            0.218119   0.103112   2.115 0.034399 *  
#> factor(year)67  0.007775   0.166325   0.047 0.962717    
#> factor(year)68  0.037799   0.152359   0.248 0.804061    
#> factor(year)69  0.334703   0.163799   2.043 0.041016 *  
#> factor(year)70  0.628643   0.246846   2.547 0.010875 *  
#> factor(year)71  0.444610   0.186188   2.388 0.016942 *  
#> factor(year)73  0.439027   0.166866   2.631 0.008513 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      -6.4195 
#> Adj. R-squared: -6.5390 
#> Wald chi2(12):  4.4 (p = 0.0000)
#> Root MSE:       1.1676 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(2) = 5.90 (p = 0.0524)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           2.72 
#>   Kleibergen-Paap rk Wald F:     2.93 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       19.93 
#>      15%  maximal IV size       11.59 
#>      20%  maximal IV size       8.75 
#>      25%  maximal IV size       7.25 
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(1) = 1.56 (p = 0.2111)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(2,744) = 46.95 (p < 2.2e-16)
#>   Anderson-Rubin Wald Chi-sq(2) = 95.66 (p < 2.2e-16)
#>   Stock-Wright LM S Chi-sq(2) = 71.33 (p = 0.0000)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 56.06 (p = 0.0000)
#>   Tested: iq
#> 
#> Redundancy test (LM):
#>   Chi-sq(1) = 0.00 (p = 0.9665)
#>   Tested: mrt
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   iq                  2.93    0.0539      0.0073    0.0073      2.93      2.93
#> 
#> Instrumented:          iq 
#> Included instruments:  s, expr, tenure, rns, smsa, factor(year)67, factor(year)68, factor(year)69, factor(year)70, factor(year)71, factor(year)73 
#> Excluded instruments:  age, mrt 
#> 
```
