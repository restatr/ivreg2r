# Klein (1950) Model I Macroeconomic Data

Annual U.S. national-accounts aggregates for 1920–1941, from Klein's
Model I of the U.S. economy: consumption, profits, wages, investment,
capital stock, government spending, and taxes. Klein, L.R. (1950).
*Economic Fluctuations in the United States, 1921–1941*. Cowles
Commission Monograph No. 11. Wiley.

## Usage

``` r
klein
```

## Format

A data frame with 22 observations and 12 variables:

- yr:

  Calendar year (1920–1941). The time variable: use as `tvar` with
  [`l()`](https://restatr.com/ivreg2r/reference/ts-operators.md)/[`d()`](https://restatr.com/ivreg2r/reference/ts-operators.md).

- consump:

  Consumption.

- profits:

  Private profits.

- wagepriv:

  Private wage bill.

- invest:

  Investment.

- capital1:

  Lagged value of capital stock (a primitive; no contemporaneous capital
  column exists to derive it from).

- totinc:

  Total income/demand.

- wagegovt:

  Government wage bill.

- govt:

  Government spending.

- taxnetx:

  Indirect business taxes plus net exports.

- wagetot:

  Total U.S. wage bill.

- year:

  Calendar year minus 1931 (a linear trend, -11 to 10), used as an
  instrument – not the time variable (see `yr`).

## Source

Klein, L.R. (1950). *Economic Fluctuations in the United States,
1921–1941*. Cowles Commission Monograph No. 11. New York: Wiley.

Distributed via Stata's `webuse klein`.

## Details

**Two columns both encode "year" – do not confuse them.** `yr` is the
calendar year (1920–1941) and is the time variable to pass as `tvar`
when using the time-series operators
[`l()`](https://restatr.com/ivreg2r/reference/ts-operators.md) /
[`d()`](https://restatr.com/ivreg2r/reference/ts-operators.md). `year`
is calendar year minus 1931 (a linear trend ranging from -11 to 10) and
is used as an *instrument* (a trend term) in the help-file examples
below – it is not a second copy of the time index.

The upstream Stata dataset (`webuse klein`) also ships two precomputed
one-period lags, `profits1` (= `L.profits`) and `totinc1` (=
`L.totinc`). Those columns are deliberately dropped here: this package's
[`l()`](https://restatr.com/ivreg2r/reference/ts-operators.md)
time-series operator computes the same lags directly from `profits` and
`totinc` given `tvar = "yr"`, so shipping precomputed lag columns would
be redundant and could drift out of sync with
[`l()`](https://restatr.com/ivreg2r/reference/ts-operators.md).
`capital1` (the lagged capital stock) is kept because it is a primitive
in this dataset – there is no contemporaneous `capital` column to lag it
from – and it appears directly in the instrument lists below.

## See also

Other ivreg2r datasets:
[`abdata`](https://restatr.com/ivreg2r/reference/abdata.md),
[`card`](https://restatr.com/ivreg2r/reference/card.md),
[`cigar`](https://restatr.com/ivreg2r/reference/cigar.md),
[`griliches`](https://restatr.com/ivreg2r/reference/griliches.md),
[`grunfeld`](https://restatr.com/ivreg2r/reference/grunfeld.md),
[`mroz`](https://restatr.com/ivreg2r/reference/mroz.md),
[`nlswork`](https://restatr.com/ivreg2r/reference/nlswork.md),
[`phillips`](https://restatr.com/ivreg2r/reference/phillips.md),
[`stockwatson`](https://restatr.com/ivreg2r/reference/stockwatson.md),
[`wagepan`](https://restatr.com/ivreg2r/reference/wagepan.md)

## Examples

``` r
data(klein)

# ivreg2 help file line 1462: LIML, consumption on lagged profits, with
# profits and wagetot treated as endogenous
fit <- ivreg2(
  consump ~ l(profits, 1) | profits + wagetot |
    govt + taxnetx + year + wagegovt + capital1 + l(totinc, 1),
  data = klein, tvar = "yr", method = "liml"
)
summary(fit)
#> 
#> LIML Estimation
#> 
#> Call:
#> ivreg2(formula = consump ~ l(profits, 1) | profits + wagetot | 
#>     govt + taxnetx + year + wagegovt + capital1 + l(totinc, 1), 
#>     data = klein, method = "liml", tvar = "yr")
#> 
#> Observations: 21 
#> VCV type:     Classical (iid) 
#> lambda:       1.498746 
#> kclass:       1.498746 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   17.14766    1.84030   9.318   <2e-16 ***
#> profits       -0.22251    0.20175  -1.103   0.2701    
#> wagetot        0.82256    0.05538  14.853   <2e-16 ***
#> l(profits, 1)  0.39603    0.17360   2.281   0.0225 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.9566 
#> Adj. R-squared: 0.9489 
#> Wald chi2(3):  118.4 (p = 0.0000)
#> Root MSE:       1.3953 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(5) = 12.01 (p = 0.0347)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           2.89 
#>   Stock-Yogo critical values (LIML size):
#>      10%  maximal LIML size      4.06 
#>      15%  maximal LIML size      2.95 
#>      20%  maximal LIML size      2.63 
#>      25%  maximal LIML size      2.46 
#> 
#> Overidentification test (Sargan):
#>   Chi-sq(4) = 6.99 (p = 0.1365)
#> 
#> Anderson-Rubin overidentification:
#>   LR Chi-sq(4) = 8.497 (p = 0.0750)
#>   Linearized Chi-sq(4) = 10.474 (p = 0.0332)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(6,13) = 18.01 (p = 0.0000)
#>   Anderson-Rubin Wald Chi-sq(6) = 174.61 (p < 2.2e-16)
#>   Stock-Wright LM S Chi-sq(6) = 18.75 (p = 0.0046)
#> 
#> Endogeneity test:
#>   Chi-sq(2) = 8.98 (p = 0.0112)
#>   Tested: profits, wagetot
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   profits             2.92    0.0497      0.5742    0.5926      3.65      3.14
#>   wagetot            38.92    0.0000      0.9473    0.9777    100.64     41.88
#> 
#> Instrumented:          profits, wagetot 
#> Included instruments:  l(profits, 1) 
#> Excluded instruments:  govt, taxnetx, year, wagegovt, capital1, l(totinc, 1) 
#> 
```
