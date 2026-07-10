# Time Series IV, Panel VCE, and GMM

## Introduction

Classical and heteroskedasticity-robust standard errors both assume that
the errors are serially uncorrelated. In time-series and panel
instrumental-variables models that assumption rarely holds: shocks
persist across periods, and in panels they may also be correlated across
units in the same period. When it fails, the usual standard errors are
inconsistent, the reported weak-identification and overidentification
statistics use the wrong variance, and two-stage least squares is no
longer the efficient use of the moment conditions. This vignette covers
the tools `ivreg2r` provides for such settings: heteroskedasticity- and
autocorrelation-consistent (HAC) standard errors and their kernel
machinery, the panel variance estimators (Kiefer, Driscoll-Kraay,
Stock-Watson, and cluster-plus-kernel), and generalized method of
moments (GMM) estimation, including the continuously-updated estimator.
Each dataset is introduced where it is first used.

``` r

library(ivreg2r)
library(dplyr)
data(phillips)
data(stockwatson)
data(griliches)
data(abdata)
data(grunfeld)
data(nlswork)
data(wagepan)
```

## HAC and AC standard errors

The `phillips` dataset contains 49 annual U.S. observations (1948–1996)
on inflation and unemployment. A simple expectations-augmented Phillips
curve regresses the change in inflation on the unemployment rate. We
pass the data frame as-is:
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md) drops
incomplete cases per specification, exactly as Stata does, so each fit
uses the largest estimation sample its own variables allow. The dataset
ships several pre-built lag columns that most fits below do not use, so
running [`na.omit()`](https://rdrr.io/r/stats/na.fail.html) on the whole
data frame first would drop rows whose only missing values are in those
unused columns, silently changing the estimation sample.

Kernel-based inference is requested with `kernel`, `bw` (bandwidth), and
`tvar` (the time variable). The type of variance estimator follows the
`vcov` setting. With `vcov = "iid"` a kernel gives
autocorrelation-consistent (AC) standard errors, which assume
homoskedasticity but allow serial correlation. With `vcov = "robust"`
the same kernel gives HAC standard errors, which allow both. This
matches Stata’s `ivreg2`, where `bw`/`kernel` alone produces AC and
combining `robust` with a bandwidth produces HAC. The explicit values
`vcov = "AC"` and `vcov = "HAC"` request the same two estimators by
name, which is how the AC fit below is written. The following three fits
use the same Phillips curve and differ only in the variance estimator:

``` r

f_ph <- cinf ~ unem
fit_iid <- ivreg2(f_ph, data = phillips, small = TRUE)
fit_ac  <- ivreg2(f_ph, data = phillips, vcov = "AC",
                  kernel = "bartlett", bw = 3, tvar = "year", small = TRUE)
fit_hac <- ivreg2(f_ph, data = phillips, vcov = "robust",
                  kernel = "bartlett", bw = 3, tvar = "year", small = TRUE)
bind_rows(iid = tidy(fit_iid), AC = tidy(fit_ac), HAC = tidy(fit_hac),
          .id = "vce") |>
  select(vce, term, estimate, std.error)
#> # A tibble: 6 × 4
#>   vce   term        estimate std.error
#>   <chr> <chr>          <dbl>     <dbl>
#> 1 iid   (Intercept)    3.03      1.38 
#> 2 iid   unem          -0.543     0.230
#> 3 AC    (Intercept)    3.03      1.25 
#> 4 AC    unem          -0.543     0.210
#> 5 HAC   (Intercept)    3.03      1.39 
#> 6 HAC   unem          -0.543     0.221
```

The estimate column repeats each coefficient across the three `vce`
rows: the point estimates are identical, because the variance estimator
does not affect them, while the standard error on `unem` shifts as the
error assumptions change. The HAC fit is the specification from Stata’s
`ivreg2` help file, run on the full N = 48 sample. The Bartlett kernel
applies linearly declining weights to autocovariances up to lag
`bw - 1`; combined with `vcov = "robust"` this is the Newey-West (1987)
HAC estimator, and the help file pairs it with the equivalent
`newey cinf unem, lag(2)`.

### Choosing a kernel

`ivreg2r` supports eight kernels, matching Stata’s `ivreg2`:

| Kernel | Weights | Guaranteed PSD? | Typical use |
|----|----|----|----|
| Bartlett | Linear decline | Yes | Default choice; Newey-West HAC |
| Parzen | Piecewise cubic | Yes | Smoother than Bartlett |
| Quadratic Spectral | $`\sin/\cos`$ based | Yes | Andrews (1991) |
| Truncated | Uniform up to lag | No | All lags equally weighted; Kiefer equivalence |
| Tukey-Hanning | Cosine taper | No | Standard spectral window |
| Tukey-Hamming | Cosine taper | No | Variant of Hanning |
| Daniell | $`\sin(x)/x`$ | No | Spectral analysis |
| Tent | Tent-shaped spectral window | No | Spectral analysis |

Bartlett, Parzen, and Quadratic Spectral guarantee a positive
semi-definite covariance estimate (Andrews 1991). For the other kernels,
use the `psd` option if negative eigenvalues arise (see the PSD
corrections section below).

### Automatic bandwidth

Use `bw = "auto"` for the Newey-West (1994) plug-in bandwidth selector,
available for the Bartlett, Parzen, and Quadratic Spectral kernels.
Pairing it with the Quadratic Spectral kernel reproduces a help-file
example:

``` r

fit_auto <- ivreg2(cinf ~ unem, data = phillips,
                   kernel = "qs", bw = "auto", tvar = "year")
fit_auto$bw
#> [1] 15.66401
```

### Kernel-based VCE with IV

Kernel-based standard errors work with every estimation method.
Instrumenting unemployment with its own first three lags gives an
overidentified model with two overidentifying restrictions. The
formula’s first part is `~ 1`, which means no exogenous regressors
beyond the constant, so the whole specification is Stata’s
`ivreg2 cinf (unem = L.unem L2.unem L3.unem), ...`:

``` r

f_ph_iv <- cinf ~ 1 | unem | l(unem, 1) + l(unem, 2) + l(unem, 3)
fit_iv_hac <- ivreg2(f_ph_iv, data = phillips, tvar = "year",
                     vcov = "robust", kernel = "bartlett", bw = 3)
summary(fit_iv_hac)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = f_ph_iv, data = phillips, vcov = "robust", kernel = "bartlett", 
#>     bw = 3, tvar = "year")
#> 
#> Observations: 46 
#> VCV type:     HAC (kernel=Bartlett; bandwidth=3) 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)   2.7901     1.8519   1.507    0.132
#> unem         -0.4756     0.2987  -1.592    0.111
#> ---
#> R-squared:      0.1414 
#> Adj. R-squared: 0.1219 
#> Wald chi2(1):  2.4 (p = 0.1266)
#> Root MSE:       2.0145 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(3) = 7.23 (p = 0.0649)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           21.01 
#>   Kleibergen-Paap rk Wald F:     25.63 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       22.30 
#>      15%  maximal IV size       12.83 
#>      20%  maximal IV size       9.54 
#>      25%  maximal IV size       7.80 
#>   Stock-Yogo critical values (IV relative bias):
#>      5%  maximal IV relative bias   13.91 
#>      10%  maximal IV relative bias   9.08 
#>      20%  maximal IV relative bias   6.46 
#>      30%  maximal IV relative bias   5.39 
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(2) = 7.53 (p = 0.0231)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(3,42) = 6.92 (p = 0.0007)
#>   Anderson-Rubin Wald Chi-sq(3) = 22.75 (p = 0.0000)
#>   Stock-Wright LM S Chi-sq(3) = 6.88 (p = 0.0758)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 0.57 (p = 0.4500)
#>   Tested: unem
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   unem               25.63    0.0000      0.6001    0.6001     25.63     25.63
#> 
#> Instrumented:          unem 
#> Excluded instruments:  l(unem, 1), l(unem, 2), l(unem, 3)
```

The diagnostics adapt automatically to the kernel-based variance: the
identification tests are Kleibergen-Paap rk statistics and the
overidentification test is a kernel-based Hansen J, here 7.535 on 2
degrees of freedom (p = 0.023) — a rejection at the 5% level, and a
useful caution: a variable’s own lags are not automatically valid
instruments when the errors are serially correlated, since the same
persistence that motivates the HAC variance can reach the moment
conditions themselves. The Stock-Yogo critical values printed alongside
the weak-identification statistic are tabulated for the Cragg-Donald
statistic under i.i.d. errors; Baum, Schaffer & Stillman (2007, p. 490)
advise applying them to the Kleibergen-Paap rk Wald F only with caution.
The Kleibergen-Paap identification tests always use the Bartlett kernel
internally regardless of the kernel chosen for the variance, matching
Stata’s `ivreg2`.

### Lags and differences in the formula

The instruments above were written with
[`l()`](https://restatr.com/ivreg2r/reference/ts-operators.md), one of
two formula operators that build lags and differences on the fly against
the time variable, mirroring Stata’s `tsset`-aware `L.` and `D.`
operators. `l(x, k)` is the k-th lag (Stata’s `Lk.x`), and `d(x, k)` is
the k-th order iterated difference $`(1-L)^k x`$ (Stata’s `Dk.x`, so
`d(x, 2)` is $`x - 2\,l(x,1) + l(x,2)`$). A range such as `l(unem, 1:3)`
expands to three separate terms. See
[`?"ts-operators"`](https://restatr.com/ivreg2r/reference/ts-operators.md)
for the full reference.

Lags are computed by time *value*, not row position: `l(x, k)` at time
`t` is the value of `x` at time `t - k` within the same panel unit, and
if no observation exists at `t - k` — the start of a series, or a gap in
the time variable — the lag is `NA` and the row leaves the estimation
sample. This is why
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md) warns when
the time variable has gaps: the warning is informational, not an error,
and the lags remain correct at their true temporal distance. You can
skip this next step: the
[`l()`](https://restatr.com/ivreg2r/reference/ts-operators.md) operator
above already handles gaps, and nothing later in the vignette depends on
it. When it is more convenient to build lag columns before fitting, the
gap-safe idiom is to fill the calendar with
[`tidyr::complete()`](https://tidyr.tidyverse.org/reference/complete.html)
first, so that
[`dplyr::lag()`](https://dplyr.tidyverse.org/reference/lead-lag.html)
steps by one period rather than by one row:

``` r

library(tidyr)
phillips_lagged <- phillips |>
  complete(year = full_seq(year, 1)) |>
  arrange(year) |>
  mutate(unem_lag1 = lag(unem))
head(phillips_lagged[, c("year", "unem", "unem_lag1")], 3)
#> # A tibble: 3 × 3
#>    year  unem unem_lag1
#>   <dbl> <dbl>     <dbl>
#> 1  1948  3.80     NA   
#> 2  1949  5.90      3.80
#> 3  1950  5.30      5.90
```

On these annual data with no gaps, `l(unem, 1)` and this pre-computed
`unem_lag1` coincide; on panels with gaps, filling the calendar first is
what keeps a positional
[`lag()`](https://dplyr.tidyverse.org/reference/lead-lag.html) aligned
to calendar time. (`tidyr` is an optional Suggests, needed only for the
[`complete()`](https://tidyr.tidyverse.org/reference/complete.html) call
in this pre-computation idiom; the
[`l()`](https://restatr.com/ivreg2r/reference/ts-operators.md) and
[`d()`](https://restatr.com/ivreg2r/reference/ts-operators.md) formula
operators require neither it nor `dplyr`.)

## Panel variance estimators

Panel data introduce richer correlation structures than a single time
series. Beyond one-way and two-way clustering (covered in
[`vignette("introduction")`](https://restatr.com/ivreg2r/articles/introduction.md)
and
[`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md)),
`ivreg2r` provides three panel variance estimators, organized here by
the asymptotic regime each one suits. The help file groups its own panel
examples the same way.

### Large N, small T: the Arellano-Bond firm panel

The `abdata` panel follows 140 UK firms over 1976–1984 — many units, few
periods. The help file’s illustration for this regime estimates an
employment equation by two-step GMM with the regressors instrumented by
their own first and second differences (built here with
[`d()`](https://restatr.com/ivreg2r/reference/ts-operators.md)),
clustering on the firm with `clusters = ~ id`, which is Stata’s
`cluster(id)`:

``` r

fit_ab <- ivreg2(
  n ~ 1 | w + k + ys |
    d(w, 1) + d(k, 1) + d(ys, 1) + d(w, 2) + d(k, 2) + d(ys, 2),
  data = abdata, tvar = "year", ivar = "id",
  method = "gmm2s", clusters = ~ id
)
summary(fit_ab)
#> 
#> 2-Step GMM Estimation
#> 
#> Call:
#> ivreg2(formula = n ~ 1 | w + k + ys | d(w, 1) + d(k, 1) + d(ys, 
#>     1) + d(w, 2) + d(k, 2) + d(ys, 2), data = abdata, clusters = ~id, 
#>     method = "gmm2s", tvar = "year", ivar = "id")
#> 
#> Observations: 751 
#> VCV type:     Cluster-robust 
#>               Estimates efficient for arbitrary clustering
#> Clusters:    140 (id)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  2.69377    3.04158   0.886    0.376    
#> w           -0.33492    0.29047  -1.153    0.249    
#> k            0.71563    0.08636   8.287   <2e-16 ***
#> ys          -0.05380    0.48758  -0.110    0.912    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.8255 
#> Adj. R-squared: 0.8248 
#> Wald chi2(3):  61.1 (p < 2.2e-16)
#> Root MSE:       0.5602 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(4) = 29.79 (p = 0.0000)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           1.64 
#>   Kleibergen-Paap rk Wald F:     10.73 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV relative bias):
#>      5%  maximal IV relative bias   12.20 
#>      10%  maximal IV relative bias   7.77 
#>      20%  maximal IV relative bias   5.35 
#>      30%  maximal IV relative bias   4.40 
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(3) = 7.54 (p = 0.0565)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(6,139) = 12.42 (p = 0.0000)
#>   Anderson-Rubin Wald Chi-sq(6) = 75.69 (p = 0.0000)
#>   Stock-Wright LM S Chi-sq(6) = 51.94 (p = 0.0000)
#> 
#> Endogeneity test:
#>   Chi-sq(3) = 2.02 (p = 0.5690)
#>   Tested: w, k, ys
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   w                  54.24    0.0000      0.0628    0.0150     27.04     27.64
#>   k                  14.69    0.0000      0.0696    0.0413     11.43     26.78
#>   ys                 75.42    0.0000      0.4212    0.0801     16.94     18.94
#> 
#> Instrumented:          w, k, ys 
#> Excluded instruments:  d(w, 1), d(k, 1), d(ys, 1), d(w, 2), d(k, 2), d(ys, 2)
```

Clustering on the firm makes the moment covariance robust to arbitrary
serial correlation within a firm’s history, which is the correlation of
concern when the panel is short. The differenced instruments cost the
first two periods of each firm, so the estimation sample is N = 751
across 140 firms, and the Hansen J of 7.54 (p = 0.057) does not reject
the overidentifying restrictions. Two-step GMM is treated in full in the
GMM section below.

Kiefer (1980) offers a homoskedastic alternative for the same regime. It
assumes conditional homoskedasticity across all observations but allows
arbitrary within-unit serial correlation, and is equivalent to AC
inference with a truncated kernel spanning the full panel length $`T`$,
so that every available within-panel lag receives full weight. The
specification here drops down to the simpler uninstrumented `n ~ w + k`
to isolate the variance estimator on a plain regression; that
simplification is for exposition, not a requirement of `kiefer`, which
works with any specification:

``` r

fit_kiefer <- ivreg2(n ~ w + k, data = abdata, kiefer = TRUE,
                     tvar = "year", ivar = "id")
summary(fit_kiefer)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = n ~ w + k, data = abdata, tvar = "year", ivar = "id", 
#>     kiefer = TRUE)
#> 
#> Observations: 1,031 
#> VCV type:     Autocorrelation-consistent (Kiefer) 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  2.55693    0.43381   5.894 3.77e-09 ***
#> w           -0.36363    0.13703  -2.654  0.00796 ** 
#> k            0.81085    0.02435  33.297  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.8345 
#> Adj. R-squared: 0.8342 
#> Wald chi2(2):  553.3 (p < 2.2e-16)
#> Root MSE:       0.5455
```

The equivalence is exact. `kiefer = TRUE` sets the bandwidth to the full
span $`T`$ internally, just as Stata’s `kiefer` option does, while an
*explicitly* supplied bandwidth must be less than $`T`$ in both
programs, making $`T - 1 = 8`$ the largest explicit choice here. The two
nonetheless coincide, because the truncated kernel gives every lag up to
the bandwidth full weight and the largest lag a panel of length $`T`$
contains is $`T - 1`$ — so bandwidths $`T`$ and $`T - 1`$ select exactly
the same lags:

``` r

bw_max <- length(unique(abdata$year)) - 1
fit_tru8 <- ivreg2(n ~ w + k, data = abdata, kernel = "truncated", bw = bw_max,
                   tvar = "year", ivar = "id")
all.equal(vcov(fit_kiefer), vcov(fit_tru8))
#> [1] TRUE
```

Use Kiefer when conditional homoskedasticity is plausible but
within-unit serial correlation of arbitrary form cannot be ruled out; if
heteroskedasticity is also a concern, use Driscoll-Kraay or
cluster-plus-kernel instead. The heteroskedastic counterpart of that
equivalence also holds: clustering on the unit with robust standard
errors is the same as a robust truncated kernel at the maximum
bandwidth, since both weight every within-unit lag equally:

``` r

fit_cl_rob  <- ivreg2(n ~ w + k, data = abdata, clusters = ~ id,
                      vcov = "robust")
fit_tru_rob <- ivreg2(n ~ w + k, data = abdata, kernel = "truncated",
                      bw = bw_max, vcov = "robust",
                      tvar = "year", ivar = "id")
all.equal(vcov(fit_cl_rob), vcov(fit_tru_rob))
#> [1] TRUE
```

### Small N, large T: the Grunfeld investment panel

Driscoll and Kraay (1998) handle cross-sectional dependence together
with persistence. The estimator clusters on the time variable and
applies a kernel to scores summed within each period, so it is robust to
common shocks — business cycles, policy changes — that move all units at
once. It is a large-$`T`$ estimator and needs a reasonably long
time-series dimension, so we switch to the `grunfeld` investment panel
(10 firms over 1935–1954):

``` r

fit_dk <- ivreg2(invest ~ mvalue + kstock, data = grunfeld,
                 dkraay = 2, small = TRUE, tvar = "year", ivar = "company")
summary(fit_dk)
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
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.8124 
#> Adj. R-squared: 0.8105 
#> F(2, 19):     118.8 (p = 0.0000)
#> Root MSE:       94.4084
```

The `dkraay` argument sets the bandwidth: here Bartlett with bandwidth
2, i.e. autocovariances up to one lag. Driscoll-Kraay is equivalent to
clustering on the time variable with a kernel (Bartlett by default,
overridable through `kernel`), an equivalence the help file demonstrates
for exactly this model:

``` r

fit_dk_ck <- ivreg2(invest ~ mvalue + kstock, data = grunfeld,
                    clusters = ~ year, kernel = "bartlett", bw = 2,
                    small = TRUE, tvar = "year", ivar = "company")
all.equal(vcov(fit_dk), vcov(fit_dk_ck))
#> [1] TRUE
```

Even $`T = 20`$ is not very large, so Driscoll-Kraay results on panels
of this length should be read with caution.

### Large N, large T: the NLS young-women panel

Two-way cluster-plus-kernel accounts for both within-panel correlation
(clustering on the unit) and autocorrelated cross-panel shocks (a kernel
on time-aggregated scores). Specify both `clusters` and `kernel`. We use
the `nlswork` panel of 4711 women over 15 survey years:

``` r

fit_ck <- ivreg2(
  ln_wage ~ grade + age + ttl_exp + tenure, data = nlswork,
  clusters = ~ idcode + year, kernel = "truncated", bw = 2,
  tvar = "year", ivar = "idcode"
)
#> Warning: Time variable 'year' has 23927 gap(s) in relevant range.
summary(fit_ck)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = ln_wage ~ grade + age + ttl_exp + tenure, data = nlswork, 
#>     clusters = ~idcode + year, kernel = "truncated", bw = 2, 
#>     tvar = "year", ivar = "idcode")
#> 
#> Observations: 28,099 
#> VCV type:     Two-way cluster + kernel (Thompson; kernel=Truncated; bandwidth=2) 
#> Clusters:     4,697 (idcode), 15 (year)
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.651825   0.056445  11.548  < 2e-16 ***
#> grade        0.074417   0.003337  22.300  < 2e-16 ***
#> age         -0.005263   0.002538  -2.074   0.0381 *  
#> ttl_exp      0.029550   0.003996   7.395 1.41e-13 ***
#> tenure       0.019517   0.004359   4.477 7.57e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.3170 
#> Adj. R-squared: 0.3169 
#> Wald chi2(4):  693.3 (p = 0.0000)
#> Root MSE:       0.3949
```

The gap warning is expected: the NLS interviews were not conducted every
calendar year, so most women have gaps in their `year` sequence, and the
kernel steps by calendar time rather than by row. The variance
decomposition is $`V_{CK} = V_{cluster,i} + V_{kernel,t} - V_{HAC}`$,
where the third term is the observation-level HAC variance at the
individual-time intersection, following the two-way structure of
Cameron, Gelbach & Miller (2011). With the truncated kernel shown here,
this construction is exactly Thompson’s (2011) estimator: two-way
cluster-plus-kernel with a truncated kernel whose bandwidth is the
number of periods after which the common cross-panel shocks can be
ignored.

### Fixed-effects panels: the Stock-Watson correction

For a within-transformed fixed-effects panel regression, Stock & Watson
(2008) derive a bias-adjusted heteroskedasticity-robust variance that
corrects the degrees-of-freedom bias of the usual cluster-robust
estimator when $`T`$ is fixed and $`N`$ is large. Request it with
`sw = TRUE` and an `ivar`. Because this regime is $`N \to \infty`$ with
fixed $`T`$, a large-$`N`$, short-$`T`$ panel is the appropriate home;
we use `wagepan` (545 men over 8 years):

``` r

fit_sw <- ivreg2(lwage ~ exper + expersq + married + union,
                 data = wagepan, sw = TRUE, ivar = "nr", tvar = "year")
summary(fit_sw)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq + married + union, data = wagepan, 
#>     tvar = "year", ivar = "nr", sw = TRUE)
#> 
#> Observations: 4,360 
#> VCV type:     Stock-Watson heteroskedastic-robust 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.1177244  0.0376554  29.683   <2e-16 ***
#> exper        0.1140221  0.0107484  10.608   <2e-16 ***
#> expersq     -0.0063520  0.0006866  -9.251   <2e-16 ***
#> married      0.1584604  0.0168736   9.391   <2e-16 ***
#> union        0.1612066  0.0180438   8.934   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.0925 
#> Adj. R-squared: 0.0917 
#> Wald chi2(4):  99.2 (p < 2.2e-16)
#> Root MSE:       0.5073
```

The `sw` option forces robust standard errors and is incompatible with
clustering, kernels, Kiefer, and Driscoll-Kraay. It appears in neither
the help file nor the Baum, Schaffer & Stillman papers, and Stata’s
`ivreg2` labels it a beta feature; the correction itself is standard
(Stock & Watson 2008). It is unrelated to the Sanderson-Windmeijer “SW”
first-stage F statistic.

### When to use which

| VCE | Cross-sectional heteroskedasticity | Within-panel autocorrelation | Cross-panel dependence |
|----|----|----|----|
| Kiefer | No\* | Yes | No |
| Stock-Watson | Yes | No | No |
| Driscoll-Kraay | Yes | Yes (kernel) | Yes |
| Cluster-plus-kernel | Yes | Yes (clustering) | Yes (kernel) |
| Two-way cluster | Yes | Yes (within groups) | Yes (within groups) |

\*Kiefer assumes conditional homoskedasticity across all observations,
both across individuals and across time.

## Generalized method of moments

### When GMM helps

Two-step efficient GMM (`method = "gmm2s"`) uses an optimal weighting
matrix to exploit overidentifying restrictions more efficiently than
2SLS under non-i.i.d. errors. When the model is exactly identified, GMM
and 2SLS produce the same coefficients. Two finite-sample cautions apply
before reaching for it by default. The optimal weighting matrix is a
function of fourth moments, which are hard to estimate well, so
efficient-GMM Wald tests tend to over-reject in finite samples (Baum,
Schaffer & Stillman 2003, §2.7, citing Hayashi 2000, p. 215).
Performance also degrades as the instrument count grows relative to the
sample (Baum, Schaffer & Stillman 2007, p. 488). In small samples with
modest overidentification, robust 2SLS may be preferable.

### A HAC-GMM arc on the Stock-Watson data

We follow the worked example of Baum, Schaffer & Stillman (2007,
p. 476): a quarterly U.S. Phillips-curve equation from Stock & Watson’s
textbook dataset (168 observations, 1959–2000), regressing the change in
inflation on the unemployment rate, instrumented by lags of GDP growth,
the T-bill rate, the exchange rate, and the T-bond rate. The errors are
serially correlated, so the GMM weighting matrix is HAC (Bartlett,
bandwidth 5). In Stata terms the specification is
`ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1)`, with `~ 1` again
meaning no exogenous regressors beyond the constant. First the 2SLS
baseline, then two-step GMM:

``` r

sw_formula <- dinf ~ 1 | UR | ggdp_2 + TBILL_1 + ER_1 + TBON_1
fit_sw_2sls <- ivreg2(sw_formula, data = stockwatson)
ur_2sls <- tidy(fit_sw_2sls) |> filter(term == "UR")
ur_2sls
#> # A tibble: 1 × 7
#>   term  estimate std.error statistic p.value conf.low conf.high
#>   <chr>    <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 UR      -0.155    0.0483     -3.21 0.00134   -0.250   -0.0603
```

``` r

fit_gmm <- ivreg2(sw_formula, data = stockwatson, method = "gmm2s",
                  vcov = "robust", kernel = "bartlett", bw = 5, tvar = "date")
ur_gmm <- tidy(fit_gmm) |> filter(term == "UR")
summary(fit_gmm)
#> 
#> 2-Step GMM Estimation
#> 
#> Call:
#> ivreg2(formula = sw_formula, data = stockwatson, vcov = "robust", 
#>     method = "gmm2s", kernel = "bartlett", bw = 5, tvar = "date")
#> 
#> Observations: 158 
#> VCV type:     HAC (kernel=Bartlett; bandwidth=5) 
#>               Estimates efficient for arbitrary autocorrelation
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  0.58508    0.37240   1.571    0.116
#> UR          -0.10024    0.06346  -1.580    0.114
#> ---
#> R-squared:      0.1548 
#> Adj. R-squared: 0.1493 
#> Wald chi2(1):  2.5 (p = 0.1185)
#> Root MSE:       0.5668 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(4) = 7.95 (p = 0.0933)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           22.58 
#>   Kleibergen-Paap rk Wald F:     7.36 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       24.58 
#>      15%  maximal IV size       13.96 
#>      20%  maximal IV size       10.26 
#>      25%  maximal IV size       8.31 
#>   Stock-Yogo critical values (IV relative bias):
#>      5%  maximal IV relative bias   16.85 
#>      10%  maximal IV relative bias   10.27 
#>      20%  maximal IV relative bias   6.71 
#>      30%  maximal IV relative bias   5.34 
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(3) = 3.57 (p = 0.3119)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(4,153) = 1.41 (p = 0.2343)
#>   Anderson-Rubin Wald Chi-sq(4) = 5.81 (p = 0.2137)
#>   Stock-Wright LM S Chi-sq(4) = 3.18 (p = 0.5278)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 1.03 (p = 0.3112)
#>   Tested: UR
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   UR                  7.36    0.0000      0.3712    0.3712      7.36      7.36
#> 
#> Instrumented:          UR 
#> Excluded instruments:  ggdp_2, TBILL_1, ER_1, TBON_1
```

The first step estimates the model by 2SLS and forms the optimal
weighting matrix from the residuals; the second step re-estimates with
that weighting, and the Hansen (1982) J statistic tests the
overidentifying restrictions. Following the published narrative, the
i.i.d. 2SLS estimate of the unemployment coefficient is -0.155 and
clearly significant, while the HAC-GMM estimate falls to -0.100 (p =
0.114): once the serial correlation is accounted for, the apparent
significance disappears.

The continuously-updated estimator (CUE) completes the arc. CUE jointly
optimizes the GMM objective over the parameter vector and the weighting
matrix, and is the GMM generalization of LIML to non-i.i.d. errors.
Continuing the same equation:

``` r

fit_cue <- ivreg2(sw_formula, data = stockwatson, method = "cue",
                  vcov = "robust", kernel = "bartlett", bw = 5, tvar = "date")
ur_cue <- tidy(fit_cue) |> filter(term == "UR")
summary(fit_cue)
#> 
#> CUE Estimation
#> 
#> Call:
#> ivreg2(formula = sw_formula, data = stockwatson, vcov = "robust", 
#>     method = "cue", kernel = "bartlett", bw = 5, tvar = "date")
#> 
#> Observations: 158 
#> VCV type:     HAC (kernel=Bartlett; bandwidth=5) 
#>               Estimates efficient for arbitrary autocorrelation
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  0.29785    0.38046   0.783    0.434
#> UR          -0.04831    0.06447  -0.749    0.454
#> ---
#> R-squared:      0.0901 
#> Adj. R-squared: 0.0842 
#> Wald chi2(1):  0.6 (p = 0.4577)
#> Root MSE:       0.5881 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(4) = 7.95 (p = 0.0933)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           22.58 
#>   Kleibergen-Paap rk Wald F:     7.36 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (LIML size):
#>      10%  maximal LIML size      5.44 
#>      15%  maximal LIML size      3.87 
#>      20%  maximal LIML size      3.30 
#>      25%  maximal LIML size      2.98 
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(3) = 2.79 (p = 0.4246)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(4,153) = 1.41 (p = 0.2343)
#>   Anderson-Rubin Wald Chi-sq(4) = 5.81 (p = 0.2137)
#>   Stock-Wright LM S Chi-sq(4) = 3.18 (p = 0.5278)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 1.03 (p = 0.3112)
#>   Tested: UR
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   UR                  7.36    0.0000      0.3712    0.3712      7.36      7.36
#> 
#> Instrumented:          UR 
#> Excluded instruments:  ggdp_2, TBILL_1, ER_1, TBON_1
```

CUE moves the coefficient further in the same direction, to -0.048 (p =
0.454), with the Hansen J at 2.793 (p = 0.425) leaving the
overidentifying restrictions comfortably unrejected.

Like LIML, CUE has no finite moments of any order regardless of the
degree of overidentification, and its finite-sample distribution can
have very heavy tails (Hansen, Heaton & Yaron 1996). Its objective can
have multiple local optima, and Baum, Schaffer & Stillman (2007, p. 479)
warn that CUE estimation can be slow and problematic when the number of
parameters is substantial. `ivreg2r` uses BFGS optimization with
Nelder-Mead refinement, starting from the 2SLS estimates under i.i.d.
weighting and from the two-step GMM estimates under a robust, cluster,
or kernel variance, matching Stata’s documented behavior. On
ill-conditioned models the objective can possess economically degenerate
minima with negative R-squared, a documented property of CUE itself
(Hausman, Lewis, Menzel & Newey 2011); when it does, both Stata and
`ivreg2r` converge to the same degenerate solution, so inspect the J
statistic and R-squared before trusting a CUE fit.

### Reusing the weighting matrix

GMM works with every variance type. On the Griliches (1976) wage
equation — 758 young men from the U.S. National Longitudinal Survey, the
help file’s own heteroskedastic-GMM illustration following Hayashi
(2000, p. 255) — two-step GMM with heteroskedasticity-robust weighting
estimates the return to measured ability (`iq`) instrumented by test
scores and family background. Here `factor(year)` expands to a set of
year dummies, R’s spelling of Stata’s `i.year`:

``` r

gril_formula <- lw ~ s + expr + tenure + rns + smsa + factor(year) |
  iq | med + kww + age + mrt
fit_g_gmm <- ivreg2(gril_formula, data = griliches,
                    method = "gmm2s", vcov = "robust")
tidy(fit_g_gmm) |> filter(term == "iq")
#> # A tibble: 1 × 7
#>   term  estimate std.error statistic p.value conf.low conf.high
#>   <chr>    <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 iq    -0.00140   0.00411    -0.341   0.733 -0.00946   0.00666
```

Every fit stores the estimated moment-condition covariance as `fit$S`,
and every fit with a defined weighting matrix — all but LIML and
k-class, exactly as in Stata — stores it as `fit$W`, matching `e(S)` and
`e(W)`. This makes the help file’s weighting-matrix reuse workflow
direct. A robust 2SLS fit of the same equation computes the very moment
covariance the two-step estimator uses in its second step, so passing
that `S` back through `smatrix` reproduces the GMM fit exactly:

``` r

fit_g_2sls <- ivreg2(gril_formula, data = griliches, vcov = "robust")
refit <- ivreg2(gril_formula, data = griliches, method = "gmm2s",
                vcov = "robust", smatrix = fit_g_2sls$S)
all.equal(fit_g_2sls$S, fit_g_gmm$S)
#> [1] TRUE
all.equal(coef(refit), coef(fit_g_gmm))
#> [1] TRUE
```

Supplying `wmatrix` instead sets the *first-step* weighting matrix
rather than fixing the second-step `S`, so the two-step estimator still
recomputes `S` from the first-step residuals and the result need not
match. For the reuse workflow, `smatrix` is the argument to use. For
two-step GMM and CUE fits, `fit$W` is simply the inverse of the defining
`fit$S`.

### Centering the moment conditions

The `center` option subtracts the mean of the moment conditions before
forming the covariance matrix `S` that feeds the variance and the
diagnostics, which can improve the finite-sample behavior of the
standard errors and the J test. No canonical example exists — the option
postdates the help-file examples and both Baum-Schaffer-Stillman papers
— so we show its effect directly on the Griliches GMM fit, reading the
Hansen J from each fit’s
[`diagnostics()`](https://restatr.com/ivreg2r/reference/diagnostics.md)
table (the package’s one-row-per-test accessor, the diagnostic
counterpart to
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) for
coefficients):

``` r

fit_g_cen <- ivreg2(gril_formula, data = griliches,
                    method = "gmm2s", vcov = "robust", center = TRUE)
bind_rows(uncentered = tidy(fit_g_gmm), centered = tidy(fit_g_cen),
          .id = "fit") |>
  filter(term == "iq") |>
  select(fit, term, estimate, std.error)
#> # A tibble: 2 × 4
#>   fit        term  estimate std.error
#>   <chr>      <chr>    <dbl>     <dbl>
#> 1 uncentered iq    -0.00140   0.00411
#> 2 centered   iq    -0.00157   0.00411
overid_g_gmm <- diagnostics(fit_g_gmm) |> filter(test == "overid")
overid_g_cen <- diagnostics(fit_g_cen) |> filter(test == "overid")
overid_g_gmm
#> # A tibble: 1 × 7
#>   test   test_name statistic    df   df2  p_value tested_vars
#>   <chr>  <chr>         <dbl> <int> <int>    <dbl> <chr>      
#> 1 overid Hansen J       74.2     3    NA 5.47e-16 NA
overid_g_cen
#> # A tibble: 1 × 7
#>   test   test_name statistic    df   df2  p_value tested_vars
#>   <chr>  <chr>         <dbl> <int> <int>    <dbl> <chr>      
#> 1 overid Hansen J       82.2     3    NA 1.03e-17 NA
coef_shift <- max(abs(coef(fit_g_cen) - coef(fit_g_gmm)))
coef_shift
#> [1] 0.004038233
```

Because two-step GMM uses `S` as its second-step weighting matrix,
centering changes the coefficients as well as the inference — here by at
most 0.004 — and shifts the J statistic from 74.2 to 82.2, with a slight
change in the standard errors on the `iq` row shown above. For
estimators whose coefficients do not depend on `S`, such as 2SLS with a
centered robust variance, only the standard errors and the J statistic
move.

## Partialling out regressors

The `partial` argument removes specified exogenous regressors from every
variable — dependent, endogenous, and excluded instruments — by
Frisch-Waugh-Lovell projection before estimation, which absorbs
high-dimensional fixed effects without carrying them as reported
coefficients. The coefficients on partialled-out variables are not
reported, but the remaining coefficients are numerically identical to
including those variables directly, for 2SLS, LIML, and two-step GMM
(the identity does not extend to CUE or to iterated GMM beyond two
steps).

Partialling also has a substantive role that the Griliches wage equation
illustrates. Clustering that equation on `year` gives far fewer clusters
than moment conditions, so the cluster-robust moment covariance is
singular and the overidentification statistics cannot be computed. The
warnings are the point, not a malfunction: they report exactly the rank
deficiency that Baum, Schaffer & Stillman (2007, pp. 484–485) describe.
This equation instruments `iq` with `med`, `kww`, and `age`, without the
`mrt` term used in the GMM example above. Partialling out the year
dummies leaves the reported coefficients unchanged but, as the counts
after the fit show, does not by itself repair the deficiency:

``` r

gril_cl <- lw ~ s + expr + tenure + rns + smsa + factor(year) |
  iq | med + kww + age
fit_full    <- ivreg2(gril_cl, data = griliches, clusters = ~ year)
#> Warning: Hansen J statistic not computed; the moment covariance or Hessian
#> matrix is singular.
#> Warning: Stock-Wright: Omega is rank-deficient; S statistic not computed.
fit_partial <- ivreg2(gril_cl, data = griliches, clusters = ~ year,
                      partial = "factor(year)")
#> Warning: Hansen J statistic not computed; the moment covariance or Hessian
#> matrix is singular.
#> Warning: Stock-Wright: Omega is rank-deficient; S statistic not computed.
shared <- c("s", "expr", "tenure", "rns", "smsa", "iq")
all.equal(coef(fit_full)[shared], coef(fit_partial)[shared])
#> [1] TRUE
```

The two fits agree on the shared coefficients to floating-point
precision, confirming the projection identity even though both leave the
Hansen J uncomputed. The counts explain why partialling alone does not
help: the full fit has 15 moment conditions against 7 clusters, and the
projection — which absorbs the year dummies and the constant from both
the regressor and instrument sets — still leaves 8 moment conditions
against the same 7 clusters. Two-step GMM is stricter than 2SLS here: it
must invert the moment covariance, so with more moment conditions than
clusters it fails outright, reproducing Stata’s `r(506)` error rather
than proceeding with a singular weighting matrix:

``` r

ivreg2(gril_cl, data = griliches, clusters = ~ year,
       partial = "factor(year)", method = "gmm2s")
#> Error:
#> ! The estimated covariance matrix of the moment conditions is not of full rank, so the optimal GMM weighting matrix is not unique. Too few clusters or too many instruments is the usual cause; reduce the instrument set or change the VCE.
```

Partialling out one further regressor brings the instrument count down
to the number of clusters, and two-step GMM becomes feasible:

``` r

fit_g_feas <- ivreg2(gril_cl, data = griliches, clusters = ~ year,
                     partial = c("factor(year)", "rns"), method = "gmm2s")
tidy(fit_g_feas) |> filter(term == "iq")
#> # A tibble: 1 × 7
#>   term  estimate std.error statistic p.value conf.low conf.high
#>   <chr>    <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 iq     0.00183   0.00387     0.472   0.637 -0.00576   0.00942
diagnostics(fit_g_feas) |> filter(test == "overid")
#> # A tibble: 1 × 7
#>   test   test_name statistic    df   df2 p_value tested_vars
#>   <chr>  <chr>         <dbl> <int> <int>   <dbl> <chr>      
#> 1 overid Hansen J       5.35     2    NA  0.0689 NA
```

With the instrument count now matched to the 7 year clusters, the Hansen
J is computed and the `iq` coefficient is reported. The special tokens
`"_cons"` (partial out the constant, equivalent to demeaning every
variable) and `"_all"` (partial out all exogenous regressors) are also
accepted. Partialling adjusts `sdofminus` for the degrees of freedom it
consumes; `nopartialsmall = TRUE` suppresses that adjustment. When fixed
effects have been partialled out *manually*, outside
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md), the
`dofminus` and `sdofminus` arguments make the same adjustment, with a
self-verifying demonstration in the degrees-of-freedom section of
[`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md).

## PSD corrections

Kernels outside the Bartlett/Parzen/Quadratic-Spectral family do not
guarantee a positive semi-definite estimate of `S`, the covariance of
the orthogonality conditions that feeds both the variance and the
chi-squared diagnostics. When `S` has negative eigenvalues, `psd`
applies an eigenvalue correction: `psd = "psd0"` sets negative
eigenvalues to zero (Politis 2007), and `psd = "psda"` replaces them
with their absolute values (Stock & Watson 2008, Remark 8).
Driscoll-Kraay with its default Bartlett kernel is positive
semi-definite by construction, but the truncated kernel carries no such
guarantee. On `wagepan`, Driscoll-Kraay with a truncated kernel produces
an indefinite `S` — a genuinely binding case, and here `S` is also
singular, so the Hansen J and Stock-Wright (2000) S statistics cannot be
computed and are reported as `NA`, exactly as Stata leaves them blank:

``` r

psd_formula <- lwage ~ exper + expersq + married + union | hours | educ + black
fit_nopsd <- ivreg2(psd_formula, data = wagepan, dkraay = 2,
                    kernel = "truncated", tvar = "year", ivar = "nr")
#> Warning: Hansen J statistic not computed; the moment covariance or Hessian
#> matrix is singular.
#> Warning: Stock-Wright: Omega is rank-deficient; S statistic not computed.
min(eigen(fit_nopsd$S, symmetric = TRUE)$values)
#> [1] -212.1851
```

Refitting with `psd = "psda"` repairs `S`; a warning reports each
corrected matrix, where Stata performs the same correction silently:

``` r

fit_psda <- ivreg2(psd_formula, data = wagepan, dkraay = 2,
                   kernel = "truncated", tvar = "year", ivar = "nr",
                   psd = "psda")
#> Warning: The covariance matrix was not positive semidefinite; 2 negative
#> eigenvalues corrected via the 'psda' method.
#> Warning: The covariance matrix was not positive semidefinite; 2 negative
#> eigenvalues corrected via the 'psda' method.
#> Warning: Hansen J statistic not computed; the moment covariance or Hessian
#> matrix is singular.
#> Warning: The covariance matrix was not positive semidefinite; 2 negative
#> eigenvalues corrected via the 'psda' method.
#> Warning: The covariance matrix was not positive semidefinite; 1 negative
#> eigenvalue corrected via the 'psda' method.
min(eigen(fit_psda$S, symmetric = TRUE)$values)
#> [1] 9.457119e-14
```

The correction is applied to `S` before the coefficient covariance is
assembled from it, exactly as in Stata, so the standard errors change
while the coefficients do not:

``` r

bind_rows(uncorrected = tidy(fit_nopsd), psda = tidy(fit_psda),
          .id = "fit") |>
  select(fit, term, estimate, std.error)
#> # A tibble: 12 × 4
#>    fit         term        estimate std.error
#>    <chr>       <chr>          <dbl>     <dbl>
#>  1 uncorrected (Intercept) -9.20     1.73    
#>  2 uncorrected hours        0.00622  0.000932
#>  3 uncorrected exper       -0.664    0.0997  
#>  4 uncorrected expersq      0.0347   0.00590 
#>  5 uncorrected married     -0.860    0.109   
#>  6 uncorrected union        0.693    0.194   
#>  7 psda        (Intercept) -9.20     1.73    
#>  8 psda        hours        0.00622  0.000938
#>  9 psda        exper       -0.664    0.112   
#> 10 psda        expersq      0.0347   0.00773 
#> 11 psda        married     -0.860    0.441   
#> 12 psda        union        0.693    0.214
```

Two quantities are deliberately untouched, again matching Stata: the
conventional `vcov = "iid"` variance, in which `S` plays no role, and
the Kleibergen-Paap identification statistics, which never receive the
`psd` option.

## Suppressing identification statistics

For simulation loops where the identification diagnostics are not
needed, `noid = TRUE` suppresses the underidentification and
weak-identification statistics (and the redundancy test), saving the
rank computations they require:

``` r

fit_noid <- ivreg2(f_ph_iv, data = phillips, tvar = "year", noid = TRUE)
diagnostics(fit_noid)
#> # A tibble: 5 × 7
#>   test                test_name        statistic    df   df2 p_value tested_vars
#>   <chr>               <chr>                <dbl> <int> <int>   <dbl> <chr>      
#> 1 overid              Sargan              6.93       2    NA 0.0313  NA         
#> 2 anderson_rubin_f    Anderson-Rubin …    3.50       3    42 0.0236  NA         
#> 3 anderson_rubin_chi2 Anderson-Rubin …   11.5        3    NA 0.00932 NA         
#> 4 stock_wright        Stock-Wright LM…    9.20       3    NA 0.0268  NA         
#> 5 endogeneity         Endogeneity         0.0962     1    NA 0.756   unem
```

The `underid`, `weak_id`, and Stock-Yogo rows are absent from the
returned tibble — the whole identification block is dropped from
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) as well —
while the `overid` test and the other non-identification diagnostics are
still reported.

## Estimator equivalences

Several estimators coincide under specific conditions, which is useful
for verification.

| Equivalence | Conditions | Where shown |
|----|----|----|
| Driscoll-Kraay = cluster(time) + kernel | Same kernel and bandwidth | In code above (grunfeld) |
| Kiefer = truncated kernel | Panel data, AC, `bw = T - 1` (all within-panel lags weighted one) | In code above (abdata) |
| Cluster-robust = truncated kernel, robust | Panel data, `bw = T - 1`, `vcov = "robust"` | In code above (abdata) |
| LIML + COVIV = CUE (COVIV: `coviv = TRUE`, the 2SLS-bread covariance for LIML) | Conditional-homoskedasticity (i.i.d.) weighting; agreement up to optimizer tolerance | Demonstrated in [`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md) |
| LIML = 2SLS | Exactly identified | Stated in [`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md) |
| GMM2S = 2SLS | Exactly identified | Stated in the GMM section above |

## References

- Andrews, D.W.K. (1991). “Heteroskedasticity and Autocorrelation
  Consistent Covariance Matrix Estimation.” *Econometrica*, 59(3),
  817–858.
- Baum, C.F., Schaffer, M.E., and Stillman, S. (2003). “Instrumental
  Variables and GMM: Estimation and Testing.” *Stata Journal*, 3(1),
  1–31.
- Baum, C.F., Schaffer, M.E., and Stillman, S. (2007). “Enhanced
  Routines for Instrumental Variables/Generalized Method of Moments
  Estimation and Testing.” *Stata Journal*, 7(4), 465–506.
- Cameron, A.C., Gelbach, J.B., and Miller, D.L. (2011). “Robust
  Inference with Multiway Clustering.” *Journal of Business & Economic
  Statistics*, 29(2), 238–249.
- Driscoll, J.C. and Kraay, A.C. (1998). “Consistent Covariance Matrix
  Estimation with Spatially Dependent Panel Data.” *Review of Economics
  and Statistics*, 80(4), 549–560.
- Griliches, Z. (1976). “Wages of Very Young Men.” *Journal of Political
  Economy*, 84(4, Part 2), S69–S85.
- Hansen, L.P. (1982). “Large Sample Properties of Generalized Method of
  Moments Estimators.” *Econometrica*, 50(4), 1029–1054.
- Hansen, L.P., Heaton, J., and Yaron, A. (1996). “Finite-Sample
  Properties of Some Alternative GMM Estimators.” *Journal of Business &
  Economic Statistics*, 14(3), 262–280.
- Hausman, J., Lewis, R., Menzel, K., and Newey, W. (2011). “Properties
  of the CUE Estimator and a Modification with Moments.” *Journal of
  Econometrics*, 165(1), 45–57.
- Hayashi, F. (2000). *Econometrics*. Princeton University Press.
- Kiefer, N.M. (1980). “Estimation of Fixed Effect Models for Time
  Series of Cross-Sections with Arbitrary Intertemporal Covariance.”
  *Journal of Econometrics*, 14(2), 195–202.
- Kleibergen, F. and Paap, R. (2006). “Generalized Reduced Rank Tests
  Using the Singular Value Decomposition.” *Journal of Econometrics*,
  133(1), 97–126.
- Newey, W.K. and West, K.D. (1987). “A Simple, Positive Semi-Definite,
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix.”
  *Econometrica*, 55(3), 703–708.
- Newey, W.K. and West, K.D. (1994). “Automatic Lag Selection in
  Covariance Matrix Estimation.” *Review of Economic Studies*, 61(4),
  631–653.
- Politis, D.N. (2007). “Higher-Order Accurate, Positive Semi-Definite
  Estimation of Large-Sample Covariance and Spectral Density Matrices.”
  Working paper, UC San Diego; published (2011) in *Econometric Theory*,
  27(4), 703–744.
- Stock, J.H. and Watson, M.W. (2008). “Heteroskedasticity-Robust
  Standard Errors for Fixed Effects Panel Data Regression.”
  *Econometrica*, 76(1), 155–174.
- Stock, J.H. and Wright, J.H. (2000). “GMM with Weak Identification.”
  *Econometrica*, 68(5), 1055–1096.
- Stock, J.H. and Yogo, M. (2005). “Testing for Weak Instruments in
  Linear IV Regression.” In D.W.K. Andrews and J.H. Stock (Eds.),
  *Identification and Inference for Econometric Models: Essays in Honor
  of Thomas Rothenberg*. Cambridge University Press.
- Thompson, S.B. (2011). “Simple Formulas for Standard Errors that
  Cluster by Both Firm and Time.” *Journal of Financial Economics*,
  99(1), 1–10.
