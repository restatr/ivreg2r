# Card (1995) College Proximity Dataset

Cross-sectional data from the National Longitudinal Survey of Young Men
(1966–1981) used by Card (1995) to estimate the return to schooling
using college proximity as an instrument for education.

## Usage

``` r
card
```

## Format

A data frame with 3,010 observations and 33 variables:

- nearc2:

  Grew up near a 2-year college (binary).

- nearc4:

  Grew up near a 4-year college (binary).

- educ:

  Years of education.

- age:

  Age in years (1976).

- fatheduc:

  Father's years of education.

- motheduc:

  Mother's years of education.

- weight:

  NLS sampling weight.

- momdad14:

  Lived with both parents at age 14 (binary).

- sinmom14:

  Lived with single mother at age 14 (binary).

- step14:

  Lived with stepparent at age 14 (binary).

- reg661:

  Region dummy: New England (1966).

- reg662:

  Region dummy: Middle Atlantic (1966).

- reg663:

  Region dummy: East North Central (1966).

- reg664:

  Region dummy: West North Central (1966).

- reg665:

  Region dummy: South Atlantic (1966).

- reg666:

  Region dummy: East South Central (1966).

- reg667:

  Region dummy: West South Central (1966).

- reg668:

  Region dummy: Mountain (1966).

- reg669:

  Region dummy: Pacific (1966).

- south66:

  Lived in the South in 1966 (binary).

- black:

  Black (binary).

- smsa:

  Lives in SMSA (binary, 1976).

- south:

  Lives in the South (binary, 1976).

- smsa66:

  Lived in SMSA in 1966 (binary).

- wage:

  Hourly wage (cents, 1976).

- enroll:

  Enrolled in school in 1976 (binary).

- KWW:

  Knowledge of the World of Work test score.

- IQ:

  IQ score.

- married:

  Marital status (1976): coded 1–6 in the source data (1 = married; the
  remaining codes distinguish other marital statuses). Not a binary
  indicator; recode before use as a regressor.

- libcrd14:

  Had a library card at age 14 (binary).

- exper:

  Years of labor market experience (`age - educ - 6`).

- lwage:

  Log hourly wage.

- expersq:

  Experience squared (`exper^2`).

## Source

Card, D. (1995). "Using Geographic Variation in College Proximity to
Estimate the Return to Schooling." In L.N. Christofides, E.K. Grant, and
R. Swidinsky (Eds.), *Aspects of Labour Market Behaviour: Essays in
Honour of John Vanderkamp*. University of Toronto Press.

## Examples

``` r
data(card)
# IV regression: instrument education with college proximity.
# Control set from Card (1995, Table 2) / Wooldridge (2020), Example 15.4.
fit <- ivreg2(
  lwage ~ exper + expersq + black + smsa + south + smsa66 +
    reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669 |
    educ | nearc4,
  data = card
)
summary(fit)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq + black + smsa + south + 
#>     smsa66 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + 
#>     reg668 + reg669 | educ | nearc4, data = card)
#> 
#> Observations: 3,010 
#> VCV type:     Classical (iid) 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  3.6661519  0.9223682   3.975 7.05e-05 ***
#> educ         0.1315038  0.0548174   2.399 0.016442 *  
#> exper        0.1082711  0.0235956   4.589 4.46e-06 ***
#> expersq     -0.0023349  0.0003326  -7.020 2.22e-12 ***
#> black       -0.1467758  0.0537564  -2.730 0.006326 ** 
#> smsa         0.1118084  0.0315777   3.541 0.000399 ***
#> south       -0.1446715  0.0272120  -5.316 1.06e-07 ***
#> smsa66       0.0185311  0.0215511   0.860 0.389862    
#> reg662       0.1007677  0.0375854   2.681 0.007340 ** 
#> reg663       0.1482588  0.0367162   4.038 5.39e-05 ***
#> reg664       0.0498971  0.0436234   1.144 0.252701    
#> reg665       0.1462719  0.0469387   3.116 0.001832 ** 
#> reg666       0.1629029  0.0517714   3.147 0.001652 ** 
#> reg667       0.1345722  0.0492708   2.731 0.006309 ** 
#> reg668      -0.0830770  0.0591734  -1.404 0.160332    
#> reg669       0.1078142  0.0417024   2.585 0.009729 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.2382 
#> Adj. R-squared: 0.2343 
#> Wald chi2(15):  51.0 (p < 2.2e-16)
#> Root MSE:       0.3873 
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           13.26 
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       16.38 
#>      15%  maximal IV size       8.96 
#>      20%  maximal IV size       6.66 
#>      25%  maximal IV size       5.53 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(1) = 13.27 (p = 0.0003)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(1,2994) = 5.42 (p = 0.0200)
#>   Anderson-Rubin Wald Chi-sq(1) = 5.44 (p = 0.0196)
#>   Stock-Wright LM S Chi-sq(1) = 5.43 (p = 0.0197)
#> 
#> Overidentification test (Sargan):  excluded (exactly identified)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 1.17 (p = 0.2786)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ               13.26    0.0003      0.0044    0.0044     13.26     13.26
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq, black, smsa, south, smsa66, reg662, reg663, reg664, reg665, reg666, reg667, reg668, reg669 
#> Excluded instruments:  nearc4 
#> 
```
