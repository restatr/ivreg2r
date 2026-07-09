# NLS Young Women Panel (Stata `[XT]` Manual Extract)

Panel data from the U.S. Bureau of Labor Statistics' National
Longitudinal Survey of Young Women, 14–24 years of age in 1968,
interviewed in survey years 1968–1988. This is the extract distributed
with Stata's `[XT]` manual (`webuse nlswork`): 28,534 person-year
observations.

## Usage

``` r
nlswork
```

## Format

A data frame with 28,534 observations and 21 variables:

- idcode:

  Person identifier. Use as `ivar`.

- year:

  Interview year (two-digit, e.g. 70 = 1970). The time variable: use as
  `tvar`.

- birth_yr:

  Birth year (two-digit).

- age:

  Age in current year.

- race:

  Race: 1 = White, 2 = Black, 3 = Other.

- msp:

  1 if married, spouse present.

- nev_mar:

  1 if never married.

- grade:

  Current grade completed.

- collgrad:

  1 if college graduate.

- not_smsa:

  1 if not in an SMSA (metropolitan area).

- c_city:

  1 if central city.

- south:

  1 if in the South.

- ind_code:

  Industry of employment (code).

- occ_code:

  Occupation (code).

- union:

  1 if union member.

- wks_ue:

  Weeks unemployed, last year.

- ttl_exp:

  Total work experience (years).

- tenure:

  Job tenure (years).

- hours:

  Usual hours worked.

- wks_work:

  Weeks worked, last year.

- ln_wage:

  Log wage (deflated by the GNP deflator).

## Source

U.S. Bureau of Labor Statistics. National Longitudinal Survey of Young
Women, 14–24 years old in 1968. Center for Human Resource Research.

Distributed via Stata's `webuse nlswork`.

## Details

`race` is coded 1 = White, 2 = Black, 3 = Other (faithful to the
upstream numeric coding; not converted to a factor). Many columns
contain missing values, consistent with the source survey.

## See also

Other ivreg2r datasets:
[`abdata`](https://restatr.com/ivreg2r/reference/abdata.md),
[`card`](https://restatr.com/ivreg2r/reference/card.md),
[`cigar`](https://restatr.com/ivreg2r/reference/cigar.md),
[`griliches`](https://restatr.com/ivreg2r/reference/griliches.md),
[`grunfeld`](https://restatr.com/ivreg2r/reference/grunfeld.md),
[`klein`](https://restatr.com/ivreg2r/reference/klein.md),
[`mroz`](https://restatr.com/ivreg2r/reference/mroz.md),
[`phillips`](https://restatr.com/ivreg2r/reference/phillips.md),
[`stockwatson`](https://restatr.com/ivreg2r/reference/stockwatson.md),
[`wagepan`](https://restatr.com/ivreg2r/reference/wagepan.md)

## Examples

``` r
data(nlswork)

# ivreg2 help file line 1618: one-way cluster on person id
fit <- ivreg2(ln_wage ~ grade + age + ttl_exp + tenure, data = nlswork,
              clusters = ~idcode)
summary(fit)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = ln_wage ~ grade + age + ttl_exp + tenure, data = nlswork, 
#>     clusters = ~idcode)
#> 
#> Observations: 28,099 
#> VCV type:     Cluster-robust 
#> Clusters:    4,697 (idcode)
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.6518251  0.0337636  19.306  < 2e-16 ***
#> grade        0.0744171  0.0021613  34.432  < 2e-16 ***
#> age         -0.0052630  0.0009445  -5.572 2.52e-08 ***
#> ttl_exp      0.0295501  0.0018370  16.086  < 2e-16 ***
#> tenure       0.0195170  0.0016316  11.962  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.3170 
#> Adj. R-squared: 0.3169 
#> Wald chi2(4):  966.0 (p < 2.2e-16)
#> Root MSE:       0.3949 
#> 

# ivreg2 help file line 1628: two-way cluster on person id and year
fit2 <- ivreg2(ln_wage ~ grade + age + ttl_exp + tenure, data = nlswork,
               clusters = ~idcode + year)
summary(fit2)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = ln_wage ~ grade + age + ttl_exp + tenure, data = nlswork, 
#>     clusters = ~idcode + year)
#> 
#> Observations: 28,099 
#> VCV type:     Cluster-robust 
#> Clusters:     4,697 (idcode), 15 (year)
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.651825   0.043107  15.121  < 2e-16 ***
#> grade        0.074417   0.002654  28.035  < 2e-16 ***
#> age         -0.005263   0.001644  -3.201  0.00137 ** 
#> ttl_exp      0.029550   0.002678  11.035  < 2e-16 ***
#> tenure       0.019517   0.003052   6.395  1.6e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.3170 
#> Adj. R-squared: 0.3169 
#> Wald chi2(4):  638.7 (p = 0.0000)
#> Root MSE:       0.3949 
#> 
```
