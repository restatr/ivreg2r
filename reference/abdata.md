# Arellano–Bond (1991) UK Employment Panel

Employment, wages, and capital stock for UK companies: an unbalanced
annual panel over 1976–1984 (140 companies, 1,031 firm-year
observations). This is the companion dataset to Arellano, M. and Bond,
S. (1991), "Some Tests of Specification for Panel Data: Monte Carlo
Evidence and an Application to Employment Equations," *Review of
Economic Studies*, 58(2), 277–297, shipped as distributed (levels plus
the log transforms and year dummies used throughout the Stata
literature).

## Usage

``` r
abdata
```

## Format

A data frame with 1,031 observations and 16 variables:

- ind:

  Industry code.

- year:

  Calendar year (1976–1984). The time variable: use as `tvar`.

- emp:

  Employment.

- wage:

  Real wage.

- cap:

  Gross capital stock.

- indoutpt:

  Industry output.

- n:

  Log employment, `log(emp)`.

- w:

  Log real wage, `log(wage)`.

- k:

  Log gross capital stock, `log(cap)`.

- ys:

  Log industry output, `log(indoutpt)`.

- yr1980:

  Year dummy: 1 if `year == 1980`.

- yr1981:

  Year dummy: 1 if `year == 1981`.

- yr1982:

  Year dummy: 1 if `year == 1982`.

- yr1983:

  Year dummy: 1 if `year == 1983`.

- yr1984:

  Year dummy: 1 if `year == 1984`.

- id:

  Firm identifier (140 distinct companies). Use as `ivar`.

## Source

Arellano, M. and Bond, S. (1991). Some tests of specification for panel
data: Monte Carlo evidence and an application to employment equations.
*Review of Economic Studies*, 58(2), 277–297.

Downloaded from <http://fmwww.bc.edu/ec-p/data/macro/abdata.dta>.

## Details

`n`, `w`, `k`, and `ys` are natural logs of `emp`, `wage`, `cap`, and
`indoutpt` respectively. `yr1980`–`yr1984` are year-dummy indicators (1
in the named year, 0 otherwise).

## See also

Other ivreg2r datasets:
[`card`](https://restatr.com/ivreg2r/reference/card.md),
[`cigar`](https://restatr.com/ivreg2r/reference/cigar.md),
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
data(abdata)

# ivreg2 help file line 1558: one-way cluster on firm id
fit <- ivreg2(n ~ w + k, data = abdata, clusters = ~id)
summary(fit)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = n ~ w + k, data = abdata, clusters = ~id)
#> 
#> Observations: 1,031 
#> VCV type:     Cluster-robust 
#> Clusters:    140 (id)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  2.55693    0.67565   3.784 0.000154 ***
#> w           -0.36363    0.21596  -1.684 0.092230 .  
#> k            0.81085    0.03248  24.963  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.8345 
#> Adj. R-squared: 0.8342 
#> Wald chi2(2):  322.5 (p < 2.2e-16)
#> Root MSE:       0.5455 
#> 
```
