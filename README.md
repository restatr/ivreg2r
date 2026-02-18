# ivreg2r

Extended instrumental variables estimation for R, inspired by Stata's
[`ivreg2`](https://ideas.repec.org/c/boc/bocode/s425401.html) (Baum,
Schaffer & Stillman). Automatically computes and reports diagnostics for weak
identification, underidentification, and overidentification alongside 2SLS
estimates.

## Installation

```r
# install.packages("devtools")
devtools::install_github("restatr/ivreg2r")
```

## Quick start

`ivreg2r` uses a three-part formula: `y ~ exog | endo | instruments`.

```r
library(ivreg2r)
data(card)

# IV regression: instrument education with college proximity
fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
              data = card)
summary(fit)
```

```
2SLS Estimation

Call:
ivreg2(formula = lwage ~ exper + expersq + black + south | educ |
    nearc4, data = card)

Observations: 3,010
VCV type:     Classical (iid)

Coefficients:
              Estimate Std. Error z value Pr(>|z|)
(Intercept)  2.3248059  0.7064713   3.291 0.000999 ***
exper        0.1439330  0.0186794   7.705 1.30e-14 ***
expersq     -0.0024018  0.0004016  -5.980 2.23e-09 ***
black       -0.0369386  0.0457925  -0.807 0.419866
south       -0.0888885  0.0256982  -3.459 0.000542 ***
educ         0.2213902  0.0408605   5.418 6.02e-08 ***
---
R-squared:      -0.1343
Adj. R-squared: -0.1362
Wald chi2(5):  83.1 (p < 2.2e-16)
Root MSE:       0.4726

Weak identification test:
  Cragg-Donald Wald F:           35.20
  Stock-Yogo critical values (IV size):
     10%  maximal IV size       16.38
     15%  maximal IV size       8.96

Underidentification test (Anderson canon. corr. LM statistic):
  Chi-sq(1) = 34.86 (p = 0.0000)

Overidentification test (Sargan):  excluded (exactly identified)

Endogeneity test:
  Chi-sq(1) = 19.17 (p = 0.0000)
  Tested: educ
```

## Formula syntax

Each variable appears exactly once. Exogenous regressors are automatically
included as instruments (you never list them twice).

| Part | Contents | Example |
|------|----------|---------|
| LHS | Dependent variable | `lwage` |
| 1st RHS | Exogenous regressors | `exper + expersq + black + south` |
| 2nd RHS | Endogenous regressors | `educ` |
| 3rd RHS | Excluded instruments | `nearc4` |

A one-part formula (`y ~ x1 + x2`) estimates OLS with the same diagnostics
infrastructure.

**Stata equivalent:**

| Stata | R |
|-------|---|
| `ivreg2 lwage exper expersq black south (educ = nearc4)` | `ivreg2(lwage ~ exper + expersq + black + south \| educ \| nearc4, data = card)` |

## Robust and cluster-robust standard errors

```r
# Heteroskedasticity-robust (equivalent to Stata's ", robust")
fit_hc0 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                   data = card, vcov = "HC0")

# With small-sample correction (equivalent to Stata's ", robust small")
fit_hc1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                   data = card, vcov = "HC1", small = TRUE)

# Cluster-robust (equivalent to Stata's ", cluster(smsa)")
fit_cl <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                  data = card, clusters = ~ smsa, small = TRUE)
```

The `vcov` argument accepts `"iid"` (default), `"HC0"`, or `"HC1"`. When
`clusters` is supplied, VCE is automatically cluster-robust. The `small`
argument controls whether t/F (instead of z/chi-squared) and finite-sample
corrections are used.

## Broom integration

`ivreg2` objects work with `tidy()`, `glance()`, and `augment()`:

```r
library(generics)

tidy(fit)
#> # A tibble: 6 x 7
#>   term        estimate std.error statistic  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)  2.32     0.706        3.29  9.99e- 4  0.940     3.71
#> 2 exper        0.144    0.0187       7.71  1.30e-14  0.107     0.181
#> 3 expersq     -0.00240  0.000402    -5.98  2.23e- 9 -0.00319  -0.00161
#> 4 black       -0.0369   0.0458      -0.807 4.20e- 1 -0.127     0.0528
#> 5 south       -0.0889   0.0257      -3.46  5.42e- 4 -0.139    -0.0385
#> 6 educ         0.221    0.0409       5.42  6.02e- 8  0.141     0.301

glance(fit)
#> # A tibble: 1 x 17
#>   r.squared adj.r.squared sigma statistic p.value    df df.residual  nobs ...
#>       <dbl>         <dbl> <dbl>     <dbl>   <dbl> <int>       <int> <int>
#> 1    -0.134        -0.136 0.473      83.1       0     5        3004  3010

augment(fit)
#> # A tibble: 3,010 x 4
#>   .fitted .resid .hat .cooksd
#>     <dbl>  <dbl> ...
```

## Diagnostics

All diagnostic tests are computed automatically at estimation time and
displayed by `summary()`. No post-estimation commands needed.

| Test | What it tells you | Statistic |
|------|-------------------|-----------|
| **Weak identification** | Are instruments strong enough? | Cragg-Donald Wald F (iid) / Kleibergen-Paap rk Wald F (robust) |
| **Stock-Yogo critical values** | Thresholds for weak identification | IV relative bias, IV size distortion |
| **Underidentification** | Is the equation identified? | Anderson LM (iid) / Kleibergen-Paap rk LM (robust) |
| **Overidentification** | Are instruments valid? | Sargan (iid) / Hansen J (robust) |
| **Endogeneity** | Are endogenous regressors actually endogenous? | C-statistic (difference-in-Sargan/J) |
| **Anderson-Rubin** | Joint significance robust to weak instruments | AR F and chi-sq |
| **First-stage F** | Excluded instrument relevance per endogenous variable | Partial F, Shea partial R-sq |
| **Sanderson-Windmeijer** | Conditional first-stage F (multiple endogenous) | SW F and chi-sq |
| **Angrist-Pischke** | Alternative conditional first-stage F | AP F and chi-sq |

See `?ivreg2` for full argument documentation.

## Stata parity

`ivreg2r` is designed as a Stata `ivreg2` translation. The Tier 1 release
matches Stata output within tight numerical tolerances (coefficients/SEs to
1e-6 relative, test statistics to 1e-4 relative) across all VCE types.

### Tier 1 (current release)

| Category | Outputs |
|----------|---------|
| **Estimation** | 2SLS coefficients, SEs, VCV, residuals, fitted values |
| **Model statistics** | R-squared, adjusted R-squared, RMSE, model F |
| **Robust VCE** | Classical (iid), HC0, HC1, one-way cluster-robust |
| **Small-sample** | t/F vs z/chi-sq, N-K denominator, cluster corrections |
| **Overidentification** | Sargan (iid), Hansen J (robust/cluster) |
| **Identification** | Anderson LM, Cragg-Donald F, Kleibergen-Paap rk LM/F |
| **Weak ID** | Stock-Yogo critical values (size/bias) |
| **Endogeneity** | C-statistic (all VCE types) |
| **Anderson-Rubin** | AR F and chi-sq (all VCE types) |
| **First stage** | F-stat, partial R-sq, Shea partial R-sq, SW F/chi-sq, AP F/chi-sq |
| **Weights** | Analytic weights (aweight) |
| **S3 methods** | coef, vcov, residuals, fitted, nobs, formula, confint, predict, summary, print |
| **Broom** | tidy, glance, augment |

### Future releases

| Feature | Tier |
|---------|------|
| LIML / Fuller / k-class | 2 |
| Two-way clustering | 2 |
| `orthog()` instrument orthogonality test | 2 |
| Frequency / probability weights | 2 |
| GMM / CUE | 3 |
| HAC / AC kernels | 3 |

## Intentional differences from Stata

| Area | Stata `ivreg2` | `ivreg2r` | Rationale |
|------|----------------|-----------|-----------|
| Formula | `ivreg2 y x (endo = z)` | `y ~ x \| endo \| z` | R convention; each variable listed once |
| Numerical method | Cross-products via Mata | QR decomposition | Better conditioning |
| Diagnostics | Some computed at post-estimation | All computed at estimation time | Always available in `glance()` |
| Sigma (weighted) | Weights normalized to sum to N | Same | Matches Stata; differs from `lm()` by `sqrt(N/sum(w))` |

## License

MIT
