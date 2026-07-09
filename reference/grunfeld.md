# Grunfeld (1958) Corporate Investment Panel

Corporate investment, market value, and capital stock for 10 large U.S.
firms observed annually over 1935–1954: a balanced panel of 10 firms
times 20 years (200 observations). Grunfeld, Y. (1958). *The
Determinants of Corporate Investment*. PhD dissertation, University of
Chicago. See Kleiber, C. and Zeileis, A. (2010) for the definitive
account of the dataset's provenance and its many circulating variants.

## Usage

``` r
grunfeld
```

## Format

A data frame with 200 observations and 6 variables (10 firms times 20
years, 1935–1954):

- company:

  Firm identifier (1–10). Use as `ivar`.

- year:

  Calendar year (1935–1954). The time variable: use as `tvar`.

- invest:

  Gross investment, current year.

- mvalue:

  Market value of the firm, prior year.

- kstock:

  Capital stock, prior year.

- time:

  Time index (1–20), a linear trend.

## Source

Grunfeld, Y. (1958). *The Determinants of Corporate Investment*. PhD
dissertation, University of Chicago.

Kleiber, C. and Zeileis, A. (2010). The Grunfeld data at 50. *German
Economic Review*, 11(4), 404–417.

Distributed via Stata's `webuse grunfeld`.

## See also

Other ivreg2r datasets:
[`abdata`](https://restatr.com/ivreg2r/reference/abdata.md),
[`card`](https://restatr.com/ivreg2r/reference/card.md),
[`cigar`](https://restatr.com/ivreg2r/reference/cigar.md),
[`griliches`](https://restatr.com/ivreg2r/reference/griliches.md),
[`klein`](https://restatr.com/ivreg2r/reference/klein.md),
[`mroz`](https://restatr.com/ivreg2r/reference/mroz.md),
[`nlswork`](https://restatr.com/ivreg2r/reference/nlswork.md),
[`phillips`](https://restatr.com/ivreg2r/reference/phillips.md),
[`stockwatson`](https://restatr.com/ivreg2r/reference/stockwatson.md),
[`wagepan`](https://restatr.com/ivreg2r/reference/wagepan.md)

## Examples

``` r
data(grunfeld)

# ivreg2 help file line 1595: Driscoll-Kraay standard errors (bw = 2),
# small-sample correction
fit <- ivreg2(
  invest ~ mvalue + kstock, data = grunfeld,
  dkraay = 2, small = TRUE, tvar = "year", ivar = "company"
)
summary(fit)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = invest ~ mvalue + kstock, data = grunfeld, small = TRUE, 
#>     tvar = "year", ivar = "company", dkraay = 2)
#> 
#> Observations: 200 
#> VCV type:     Driscoll-Kraay (kernel=Bartlett; bandwidth=2), small-sample corrected 
#> Clusters:    20 (year)
#> 
#> Coefficients:
#>              Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -42.71437   12.22308  -3.495  0.00243 ** 
#> mvalue        0.11556    0.01024  11.289 7.22e-10 ***
#> kstock        0.23068    0.04670   4.940 9.10e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> R-squared:      0.8124 
#> Adj. R-squared: 0.8105 
#> F(2, 19):     118.8 (p = 0.0000)
#> Root MSE:       94.4084 
#> 
```
