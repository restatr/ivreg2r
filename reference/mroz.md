# Mroz (1987) Female Labor Supply Dataset

Cross-sectional data on 753 married white women from the Panel Study of
Income Dynamics (PSID), 1975. Of these, 428 were working (`inlf == 1`)
with observed wages. Used by Mroz (1987) to study married women's labor
force participation and hours of work. A classic IV application
instruments education with parental education (`motheduc`, `fatheduc`).

## Usage

``` r
mroz
```

## Format

A data frame with 753 observations and 22 variables:

- inlf:

  In the labor force in 1975 (binary).

- hours:

  Hours worked in 1975.

- kidslt6:

  Number of children younger than 6.

- kidsge6:

  Number of children aged 6–18.

- age:

  Age in years.

- educ:

  Years of education.

- wage:

  Estimated hourly wage (1975 dollars).

- repwage:

  Reported hourly wage at interview (1976).

- hushrs:

  Husband's hours worked in 1975.

- husage:

  Husband's age.

- huseduc:

  Husband's years of education.

- huswage:

  Husband's hourly wage (1975 dollars).

- faminc:

  Family income (1975 dollars).

- mtr:

  Federal marginal tax rate facing the woman.

- motheduc:

  Mother's years of education.

- fatheduc:

  Father's years of education.

- unem:

  Unemployment rate in county of residence.

- city:

  Lives in SMSA (binary).

- exper:

  Years of labor market experience.

- nwifeinc:

  Non-wife household income (`faminc - wage * hours`, in thousands of
  1975 dollars).

- lwage:

  Log estimated hourly wage.

- expersq:

  Experience squared (`exper^2`).

## Source

Mroz, T.A. (1987). "The Sensitivity of an Empirical Model of Married
Women's Hours of Work to Economic and Statistical Assumptions."
*Econometrica*, 55(4), 765–799.

## See also

Other ivreg2r datasets:
[`abdata`](https://restatr.com/ivreg2r/reference/abdata.md),
[`card`](https://restatr.com/ivreg2r/reference/card.md),
[`cigar`](https://restatr.com/ivreg2r/reference/cigar.md),
[`griliches`](https://restatr.com/ivreg2r/reference/griliches.md),
[`grunfeld`](https://restatr.com/ivreg2r/reference/grunfeld.md),
[`klein`](https://restatr.com/ivreg2r/reference/klein.md),
[`nlswork`](https://restatr.com/ivreg2r/reference/nlswork.md),
[`phillips`](https://restatr.com/ivreg2r/reference/phillips.md),
[`stockwatson`](https://restatr.com/ivreg2r/reference/stockwatson.md),
[`wagepan`](https://restatr.com/ivreg2r/reference/wagepan.md)

## Examples

``` r
data(mroz)
# Restrict to working women (observed wages)
mroz_work <- subset(mroz, inlf == 1)
# IV regression: instrument education with parental education (Example 15.5,
# "Return to Education for Working Women", in Wooldridge (2020),
# Introductory Econometrics).
# These instruments (both parents' education) differ from the specification
# in ?ivreg2, which instruments educ with fatheduc alone (just-identified)
# and with age + kidslt6 + kidsge6 (overidentified).
fit <- ivreg2(lwage ~ exper + expersq | educ | motheduc + fatheduc,
              data = mroz_work)
summary(fit)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq | educ | motheduc + 
#>     fatheduc, data = mroz_work)
#> 
#> Observations: 428 
#> VCV type:     Classical (iid) 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.0481003  0.3984530   0.121 0.903915    
#> educ         0.0613966  0.0312895   1.962 0.049737 *  
#> exper        0.0441704  0.0133696   3.304 0.000954 ***
#> expersq     -0.0008990  0.0003998  -2.249 0.024543 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.1357 
#> Adj. R-squared: 0.1296 
#> Wald chi2(3):  8.1 (p = 0.0000)
#> Root MSE:       0.6716 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(2) = 88.84 (p < 2.2e-16)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           55.40 
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       19.93 
#>      15%  maximal IV size       11.59 
#>      20%  maximal IV size       8.75 
#>      25%  maximal IV size       7.25 
#> 
#> Overidentification test (Sargan):
#>   Chi-sq(1) = 0.38 (p = 0.5386)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(2,423) = 1.90 (p = 0.1505)
#>   Anderson-Rubin Wald Chi-sq(2) = 3.85 (p = 0.1459)
#>   Stock-Wright LM S Chi-sq(2) = 3.81 (p = 0.1485)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 2.81 (p = 0.0938)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ               55.40    0.0000      0.2076    0.2076     55.40     55.40
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq 
#> Excluded instruments:  motheduc, fatheduc 
#> 
```
