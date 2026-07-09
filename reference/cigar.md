# Cigarette Demand Panel (Baltagi–Levin / Baltagi–Griffin–Xiong)

Annual per-capita cigarette sales for 46 U.S. states over 1963–1992 (a
balanced panel, 1,380 state-year observations), along with the price per
pack, population, consumer price index, per-capita disposable income,
and the minimum price per pack among neighboring states used to
instrument for price endogeneity. Baltagi, B. H. and Levin, D. (1992),
"Cigarette taxation: raising revenues and reducing consumption,"
*Structural Change and Economic Dynamics*, 3(2), 321–335. See also
Baltagi, B. H., Griffin, J. M. and Xiong, W. (2000), "To pool or not to
pool: homogeneous versus heterogeneous estimators applied to cigarette
demand," *Review of Economics and Statistics*, 82(1), 117–126.

## Usage

``` r
cigar
```

## Format

A data frame with 1,380 observations and 9 variables (46 states times 30
years, 1963–1992):

- state:

  State identifier code (46 distinct U.S. states). Use as `ivar`.

- year:

  Year, coded 63–92 (i.e., 1963–1992). The time variable: use as `tvar`.

- price:

  Price per pack of cigarettes (nominal).

- pop:

  Population.

- pop16:

  Population above the age of 16.

- cpi:

  Consumer price index (1983 = 100).

- ndi:

  Per-capita nominal disposable income.

- sales:

  Cigarette sales, in packs per capita.

- pimin:

  Minimum price per pack of cigarettes in adjoining states.

## Source

Baltagi, B. H. and Levin, D. (1992). Cigarette taxation: raising
revenues and reducing consumption. *Structural Change and Economic
Dynamics*, 3(2), 321–335.

Baltagi, B. H., Griffin, J. M. and Xiong, W. (2000). To pool or not to
pool: homogeneous versus heterogeneous estimators applied to cigarette
demand. *Review of Economics and Statistics*, 82(1), 117–126.

Distributed via R package `plm`, `data(Cigar)`.

## Details

`price` and `ndi` are nominal; deflate by `cpi` for real terms.

## See also

Other ivreg2r datasets:
[`abdata`](https://restatr.com/ivreg2r/reference/abdata.md),
[`card`](https://restatr.com/ivreg2r/reference/card.md),
[`griliches`](https://restatr.com/ivreg2r/reference/griliches.md),
[`grunfeld`](https://restatr.com/ivreg2r/reference/grunfeld.md),
[`klein`](https://restatr.com/ivreg2r/reference/klein.md),
[`mroz`](https://restatr.com/ivreg2r/reference/mroz.md),
[`nlswork`](https://restatr.com/ivreg2r/reference/nlswork.md),
[`phillips`](https://restatr.com/ivreg2r/reference/phillips.md),
[`stockwatson`](https://restatr.com/ivreg2r/reference/stockwatson.md),
[`wagepan`](https://restatr.com/ivreg2r/reference/wagepan.md)

## Examples

``` r
data(cigar)

# Price-endogeneity IV specification in real (CPI-deflated) terms, with
# the neighboring-state minimum price as the excluded instrument, as in
# the GFIC empirical example.
cigar_real <- transform(
  cigar,
  lsales  = log(sales),
  lrprice = log(price / cpi),
  lrndi   = log(ndi / cpi),
  lrpimin = log(pimin / cpi)
)
fit <- ivreg2(lsales ~ lrndi | lrprice | lrpimin, data = cigar_real)
summary(fit)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lsales ~ lrndi | lrprice | lrpimin, data = cigar_real)
#> 
#> Observations: 1,380 
#> VCV type:     Classical (iid) 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  3.58676    0.11626  30.852   <2e-16 ***
#> lrprice     -0.75698    0.04233 -17.882   <2e-16 ***
#> lrndi        0.24775    0.02521   9.827   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.3165 
#> Adj. R-squared: 0.3155 
#> Wald chi2(2):  168.3 (p < 2.2e-16)
#> Root MSE:       0.1856 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(1) = 901.42 (p < 2.2e-16)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           2593.62 
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       16.38 
#>      15%  maximal IV size       8.96 
#>      20%  maximal IV size       6.66 
#>      25%  maximal IV size       5.53 
#> 
#> Overidentification test (Sargan):  excluded (exactly identified)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(1,1377) = 261.83 (p < 2.2e-16)
#>   Anderson-Rubin Wald Chi-sq(1) = 262.40 (p < 2.2e-16)
#>   Stock-Wright LM S Chi-sq(1) = 220.48 (p < 2.2e-16)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 16.86 (p = 0.0000)
#>   Tested: lrprice
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   lrprice          2593.62    0.0000      0.6532    0.6532   2593.62   2593.62
#> 
#> Instrumented:          lrprice 
#> Included instruments:  lrndi 
#> Excluded instruments:  lrpimin 
#> 
```
