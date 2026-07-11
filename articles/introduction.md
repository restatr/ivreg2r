# Instrumental Variables Estimation with ivreg2r

## Introduction

`ivreg2r` provides comprehensive instrumental variables and GMM
estimation with automatic diagnostic tests. It supports 2SLS, LIML,
Fuller, k-class, two-step efficient GMM, and continuously-updated (CUE)
estimators with classical, robust, cluster-robust, HAC, and
Driscoll-Kraay standard errors. It is an R implementation inspired by
Stata’s `ivreg2` (Baum, Schaffer & Stillman), designed for applied
researchers — especially those migrating from Stata who expect rich
diagnostics reported alongside their IV estimates.

This vignette covers the core workflow: OLS, 2SLS, diagnostics, robust
and cluster-robust inference, weights, and broom integration. For LIML,
Fuller, and k-class estimation, weak-instrument-robust inference, and
model-comparison workflows, see
[`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md).
For HAC/AC standard errors, GMM, CUE, partialling, and panel VCE
options, see
[`vignette("time-series-gmm")`](https://restatr.com/ivreg2r/articles/time-series-gmm.md).

We’ll use the Card (1995) dataset throughout this vignette. Card uses
geographic proximity to a four-year college as an instrument for years
of education in a log-wage equation, exploiting the idea that growing up
near a college reduces the cost of attending and therefore shifts
educational attainment. The exclusion restriction — that proximity
affects wages only through schooling — is an identifying *assumption*,
not something the data can fully verify: proximity could in principle be
correlated with unobserved determinants of wages, a concern Card (1995)
discusses at length.

``` r

library(ivreg2r)
library(dplyr)
data(card)
```

## OLS baseline

Give [`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md) an
ordinary regression formula with no `|` separators — the same
`y ~ x1 + x2` you would pass to
[`lm()`](https://rdrr.io/r/stats/lm.html) — and it fits by ordinary
least squares. In Stata this is the OLS special case of the same
command: `ivreg2 lwage educ exper ...` with no `(educ = ...)` block. We
use the control set from Card (1995, Table 2), which is also the
specification in Example 15.4 of Wooldridge (2020): years of schooling,
experience and its square, race, urban residence (`smsa`), current
region (`south`), and the 1966 residence indicators (`smsa66` together
with the census-region dummies). The Card data ship the census-region
indicators already expanded as separate dummies, `reg661`–`reg669` with
`reg661` the omitted base, so we list `reg662`–`reg669` directly:

``` r

fit_ols <- ivreg2(
  lwage ~ educ + exper + expersq + black + smsa + south + smsa66 +
    reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669,
  data = card
)
summary(fit_ols)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ educ + exper + expersq + black + smsa + 
#>     south + smsa66 + reg662 + reg663 + reg664 + reg665 + reg666 + 
#>     reg667 + reg668 + reg669, data = card)
#> 
#> Observations: 3,010 
#> VCV type:     Classical (iid) 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  4.6208069  0.0740352  62.414  < 2e-16 ***
#> educ         0.0746932  0.0034890  21.408  < 2e-16 ***
#> exper        0.0848320  0.0066066  12.841  < 2e-16 ***
#> expersq     -0.0022870  0.0003158  -7.242 4.41e-13 ***
#> black       -0.1990123  0.0181997 -10.935  < 2e-16 ***
#> smsa         0.1363846  0.0200470   6.803 1.02e-11 ***
#> south       -0.1479550  0.0259107  -5.710 1.13e-08 ***
#> smsa66       0.0262417  0.0193959   1.353  0.17607    
#> reg662       0.0963671  0.0358024   2.692  0.00711 ** 
#> reg663       0.1445400  0.0350309   4.126 3.69e-05 ***
#> reg664       0.0550755  0.0415465   1.326  0.18496    
#> reg665       0.1280248  0.0417281   3.068  0.00215 ** 
#> reg666       0.1405174  0.0451265   3.114  0.00185 ** 
#> reg667       0.1179810  0.0446832   2.640  0.00828 ** 
#> reg668      -0.0564361  0.0511215  -1.104  0.26961    
#> reg669       0.1185697  0.0387267   3.062  0.00220 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.2998 
#> Adj. R-squared: 0.2963 
#> Wald chi2(15):  85.5 (p < 2.2e-16)
#> Root MSE:       0.3713
```

The point estimates are identical to
[`lm()`](https://rdrr.io/r/stats/lm.html), and the model uses the same
estimation infrastructure as the IV models below, so switching between
OLS and 2SLS requires only a formula change. The default *reporting*
follows Stata’s conventions, however: z statistics rather than t
statistics, standard errors without the finite-sample N - K correction,
and Root MSE computed as sqrt(RSS/N). Setting `small = TRUE` reproduces
[`lm()`](https://rdrr.io/r/stats/lm.html)’s standard errors, t
statistics, and sigma exactly. (For *weighted* models, coefficients and
`small = TRUE` standard errors still match
[`lm()`](https://rdrr.io/r/stats/lm.html), but sigma is defined
differently — see Weighted estimation below.)

## Two-stage least squares

To estimate by two-stage least squares, mark the endogenous regressor
and its instruments with the `|` separator. Where Stata writes
`ivreg2 lwage <controls> (educ = nearc4)`, `ivreg2r` writes
`lwage ~ <controls> | educ | nearc4`: the exogenous regressors first,
then the endogenous regressor after the first `|`, then the excluded
instruments after the second. Each variable is listed once, and the
exogenous regressors are used as instruments automatically.

Here we instrument `educ` with `nearc4` (grew up near a four-year
college). The robust, first-stage, and weighted fits below all reuse
this same specification, so we write the formula out once and store it:

``` r

iv_formula <- lwage ~ exper + expersq + black + smsa + south + smsa66 +
  reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669 |
  educ | nearc4

fit_iv <- ivreg2(iv_formula, data = card)
summary(fit_iv)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = iv_formula, data = card)
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
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.2382 
#> Adj. R-squared: 0.2343 
#> Wald chi2(15):  51.0 (p < 2.2e-16)
#> Root MSE:       0.3873 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(1) = 13.27 (p = 0.0003)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           13.26 
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       16.38 
#>      15%  maximal IV size       8.96 
#>      20%  maximal IV size       6.66 
#>      25%  maximal IV size       5.53 
#> 
#> Overidentification test (Sargan):  (equation exactly identified)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(1,2994) = 5.42 (p = 0.0200)
#>   Anderson-Rubin Wald Chi-sq(1) = 5.44 (p = 0.0196)
#>   Stock-Wright LM S Chi-sq(1) = 5.43 (p = 0.0197)
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
```

All 3,010 observations are used, and the IV estimate of the return to
schooling (0.1315, versus 0.0747 by OLS) replicates the published
estimates for this specification (Card, 1995, Table 2; Wooldridge, 2020,
Example 15.4).

The formula has three groups, separated by the two `|` bars, mapping
onto the pieces of the Stata command:

| Group | Contents | In `ivreg2 lwage ... (educ = nearc4)` |
|----|----|----|
| Before the first `|` | Dependent variable, then exogenous regressors | `lwage`, then the controls |
| Between the two `|` | Endogenous regressor(s) | `educ` |
| After the second `|` | Excluded instruments | `nearc4` |

The full instrument set is the union of the exogenous regressors and the
excluded instruments; you never list the exogenous regressors twice.

## Understanding the diagnostics

The [`summary()`](https://rdrr.io/r/base/summary.html) output reports
several diagnostic tests automatically. This section walks through each
one, using the values from the fit above.

### Weak identification

The fit above reports a **Cragg-Donald Wald F** of 13.26. This statistic
tests whether the excluded instruments are *weak* — that is, only weakly
correlated with the endogenous regressors. Weak instruments lead to
biased 2SLS estimates and unreliable standard errors.

**Stock-Yogo critical values** provide formal thresholds for the F
statistic. Two types exist:

- **IV size**: the critical value such that the *maximal* rejection rate
  of a nominal 5% Wald test on the endogenous regressors’ coefficients
  is no more than 10%, 15%, 20%, or 25%. “Maximal” means the worst case
  over all possible true parameter values. For one endogenous variable
  and one instrument, the 10% maximal-size critical value is 16.38.
- **IV relative bias**: the critical value such that the *maximum*
  relative bias of 2SLS (relative to OLS bias) is no more than 5%, 10%,
  20%, or 30%. Stock and Yogo tabulate these only when the number of
  overidentifying restrictions is at least two: the measure relies on
  the finite-sample moments of the 2SLS estimator, which exist only up
  to the degree of overidentification (Stock & Yogo, 2005; Baum,
  Schaffer & Stillman, 2007, p. 490). That is why neither this
  just-identified fit nor the overidentified fit below displays them:
  both have fewer than two overidentifying restrictions.

Here the Cragg-Donald F of 13.26 falls short of the 10% maximal-size
critical value of 16.38 but clears the 15% value of 8.96 — a genuinely
borderline case. A 5% Wald test on `educ` could reject as often as 15%
of the time under the null, so the conventional t-statistic should be
read with some suspicion. This is exactly the situation the
weak-instrument-robust Anderson-Rubin test (below) is designed for.

One caution: the Stock-Yogo critical values are tabulated for the
Cragg-Donald statistic under i.i.d. errors. Under a robust, cluster, or
HAC VCE the reported statistic is instead the Kleibergen-Paap rk Wald F,
and Baum, Schaffer & Stillman (2007, p. 490) recommend applying the
i.i.d. values to it only with caution, or falling back on the Staiger &
Stock (1997) rule that F should exceed 10. See
[`?ivreg2`](https://restatr.com/ivreg2r/reference/ivreg2.md) and
[`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md)
for the fuller discussion.

### Underidentification

The fit above reports an Anderson LM statistic of 13.27 (p = 0.00027).
The underidentification test asks whether the excluded instruments have
*any* correlation with the endogenous regressors. Formally, the null
hypothesis is that the matrix of reduced-form coefficients on the
excluded instruments has rank K1 - 1, where K1 is the number of
endogenous regressors — i.e., that the equation is underidentified. With
a single endogenous regressor, as here, this reduces to the intuitive
statement that the instruments are irrelevant. Rejection is good news:
it means the model is identified.

With *several* endogenous regressors the test covers the whole system,
failing to reject if any one of them is unidentified, so the
per-regressor first-stage statistics below are the finer diagnostic.
Under i.i.d. errors this is the Anderson (1951) canonical-correlations
LM test; under heteroskedasticity or clustering it becomes the
Kleibergen-Paap (2006) rk LM statistic.

Underidentification is distinct from weak identification. The
underidentification test asks “is there *any* correlation at all?” while
the weak identification test asks “is there *enough* correlation for
reliable inference?” Here the underidentification test rejects
decisively even though the instrument is borderline-weak by the
Stock-Yogo standard — a common pattern, and a reminder that rejecting
underidentification is not by itself evidence of strong instruments.

### Overidentification

When there are more excluded instruments than endogenous regressors, the
model is overidentified and we can test whether the instruments are
valid (i.e., uncorrelated with the structural error). Adding proximity
to a *two*-year college (`nearc2`) as a second instrument gives the
specification of Computer Exercise C5, Chapter 15, in Wooldridge (2020):

``` r

fit_overid <- ivreg2(
  lwage ~ exper + expersq + black + smsa + south + smsa66 +
    reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + reg669 |
    educ | nearc2 + nearc4,
  data = card
)
summary(fit_overid)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq + black + smsa + south + 
#>     smsa66 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + 
#>     reg668 + reg669 | educ | nearc2 + nearc4, data = card)
#> 
#> Observations: 3,010 
#> VCV type:     Classical (iid) 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  3.2367114  0.8825567   3.667 0.000245 ***
#> educ         0.1570593  0.0524383   2.995 0.002743 ** 
#> exper        0.1188149  0.0227454   5.224 1.75e-07 ***
#> expersq     -0.0023565  0.0003466  -6.799 1.05e-11 ***
#> black       -0.1232778  0.0520112  -2.370 0.017778 *  
#> smsa         0.1007530  0.0314355   3.205 0.001350 ** 
#> south       -0.1431945  0.0283691  -5.048 4.48e-07 ***
#> smsa66       0.0150626  0.0222765   0.676 0.498937    
#> reg662       0.1027473  0.0391861   2.622 0.008741 ** 
#> reg663       0.1499316  0.0382896   3.916 9.01e-05 ***
#> reg664       0.0475676  0.0454799   1.046 0.295605    
#> reg665       0.1544801  0.0484336   3.190 0.001425 ** 
#> reg666       0.1729728  0.0532743   3.247 0.001167 ** 
#> reg667       0.1420355  0.0509858   2.786 0.005340 ** 
#> reg668      -0.0950611  0.0608178  -1.563 0.118041    
#> reg669       0.1029760  0.0433068   2.378 0.017415 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.1702 
#> Adj. R-squared: 0.1660 
#> Wald chi2(15):  47.1 (p < 2.2e-16)
#> Root MSE:       0.4042 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(2) = 15.79 (p = 0.0004)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           7.89 
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       19.93 
#>      15%  maximal IV size       11.59 
#>      20%  maximal IV size       8.75 
#>      25%  maximal IV size       7.25 
#> 
#> Overidentification test (Sargan):
#>   Chi-sq(1) = 1.25 (p = 0.2639)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(2,2993) = 5.24 (p = 0.0053)
#>   Anderson-Rubin Wald Chi-sq(2) = 10.55 (p = 0.0051)
#>   Stock-Wright LM S Chi-sq(2) = 10.51 (p = 0.0052)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 2.94 (p = 0.0864)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ                7.89    0.0004      0.0052    0.0052      7.89      7.89
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq, black, smsa, south, smsa66, reg662, reg663, reg664, reg665, reg666, reg667, reg668, reg669 
#> Excluded instruments:  nearc2, nearc4
```

The Sargan (1958) test (for iid models) or Hansen (1982) J test (for
robust/cluster models) tests the null hypothesis that all instruments
are valid. Here the Sargan statistic is 1.25 (p = 0.2639), so the single
overidentifying restriction is not rejected. A rejection would have
indicated that the instruments are not satisfying the required
orthogonality conditions — either because an instrument is not truly
exogenous, or because it is incorrectly excluded from the regression.
The test is a joint test of instrument validity *and* correct model
specification (Baum, Schaffer & Stillman, 2003, pp. 16–17).

The test exists only when the model is overidentified; with exact
identification there are no testable restrictions. It has low power when
instruments are weak, and no power at all against violations that leave
every instrument consistent for the same wrong parameter value, so a
failure to reject does not by itself establish validity (see
[`?ivreg2`](https://restatr.com/ivreg2r/reference/ivreg2.md) for the
fuller treatment).

To test a *subset* of instruments while maintaining the rest, use the
`orthog` argument — see
[`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md)
for details.

### Endogeneity test (C-statistic)

The endogeneity test (also called the C-statistic or
difference-in-Sargan test) asks whether the regressors treated as
endogenous are actually endogenous. The null hypothesis is that they are
exogenous. For the just-identified fit above, the statistic is 1.17 (p =
0.2786): we cannot reject that `educ` is exogenous. In that case OLS may
be preferred because it is more efficient — this fit is a live example
of the trade-off, since the IV point estimate is similar to OLS but far
less precise.

The test maintains instrument validity under both hypotheses, so it is
uninformative about endogeneity if the instruments are themselves
invalid, and it shares the low power under weak instruments noted above.
Using its result to choose between OLS and IV is a pre-test whose
selection step can badly distort the size of subsequent inference
(Guggenberger, 2010); see
[`?ivreg2`](https://restatr.com/ivreg2r/reference/ivreg2.md) and Hayashi
(2000, pp. 218–221, 232–234) for the underlying theory.

Stata’s `ivreg2` computes this test only when `endog()` is requested,
whereas `ivreg2r` reports it automatically for all endogenous
regressors; use the `endog` argument to test a specific subset (e.g.,
`endog = "educ"`).

### Anderson-Rubin test

The Anderson-Rubin (1949) test provides inference on the endogenous
regressors that remains valid even when instruments are weak. Unlike the
standard Wald or t-test on 2SLS coefficients — which can be severely
distorted by weak instruments — the AR test has correct size regardless
of instrument strength.

The null hypothesis is that the coefficients on the endogenous
regressors are jointly zero *and* the overidentifying restrictions hold.
When the model is overidentified, rejection could indicate either
nonzero coefficients or invalid overidentifying restrictions (or both).
In the just-identified case, the overidentifying restrictions are
trivially satisfied and the test is purely about the coefficients.

This serves as a useful cross-check: if the standard t-test says an
endogenous regressor is significant but the AR test does not reject,
weak instruments may be distorting inference. In the running example the
two tests agree but the gap is visible: the conventional z-test on
`educ` gives p = 0.0164, while the Anderson-Rubin chi-squared test gives
p = 0.0196 — still a rejection at the 5% level, but a weaker one,
consistent with the borderline first stage.

Robustness to weak instruments comes at a price: as instruments weaken,
the power of the AR test declines — Baum, Schaffer & Stillman (2007,
p. 491) describe this as an intuitively appealing feature, since weak
instruments *should* make confident inference harder. In overidentified
models the AR chi-squared statistic has L1 degrees of freedom (L1 =
number of excluded instruments), so power also falls as instruments are
added without strengthening identification. More powerful
weak-instrument-robust procedures of the conditional likelihood ratio
type exist for those cases;
[`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md)
discusses what this package does and does not offer and names software
that constructs weak-instrument-robust confidence sets.

The Stock-Wright (2000) S statistic tests the same joint null hypothesis
using a GMM-distance formulation. It is also automatically computed. See
[`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md)
for more on weak-instrument-robust inference.

## First-stage diagnostics

The [`summary()`](https://rdrr.io/r/base/summary.html) output reports a
table of first-stage diagnostics for each endogenous variable. Fitting
with `first_stage = TRUE` additionally stores the underlying first-stage
regressions, which
[`first_stage()`](https://restatr.com/ivreg2r/reference/first_stage.md)
returns as a named list — one `ivreg2_first_stage` object per endogenous
regressor, each supporting
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html), and
[`glance()`](https://generics.r-lib.org/reference/glance.html):

``` r

fit_fs <- ivreg2(iv_formula, data = card, first_stage = TRUE)
fs <- first_stage(fit_fs)
summary(fs$educ)
#> 
#> First-stage regression: educ ~ instruments
#> Observations: 3,010
#> VCV type:     iid
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) 16.6382529  0.2406297  69.145  < 2e-16 ***
#> exper       -0.4125334  0.0336996 -12.241  < 2e-16 ***
#> expersq      0.0008686  0.0016504   0.526 0.598728    
#> black       -0.9355287  0.0937348  -9.981  < 2e-16 ***
#> smsa         0.4021825  0.1048112   3.837 0.000127 ***
#> south       -0.0516126  0.1354284  -0.381 0.703152    
#> smsa66       0.0254805  0.1057692   0.241 0.809644    
#> reg662      -0.0786363  0.1871154  -0.420 0.674329    
#> reg663      -0.0279390  0.1833745  -0.152 0.878913    
#> reg664       0.1171820  0.2172531   0.539 0.589665    
#> reg665      -0.2726165  0.2184204  -1.248 0.212082    
#> reg666      -0.3028147  0.2370712  -1.277 0.201590    
#> reg667      -0.2168177  0.2343879  -0.925 0.355021    
#> reg668       0.5238914  0.2674749   1.959 0.050246 .  
#> reg669       0.2102710  0.2024568   1.039 0.299076    
#> nearc4       0.3198989  0.0878638   3.641 0.000276 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Root MSE:       1.9405 
#> F(1, 2994) = 13.26 (p = 0.0003)
#> Partial R2:     0.0044 
#> Shea PR2:       0.0044
```

The headline number is the partial F-statistic of the excluded
instruments, the standard measure of instrument relevance. With one
endogenous regressor and one instrument it equals the Cragg-Donald Wald
F reported above (13.26 here), so the borderline weak-identification
reading carries over directly. The
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) method
returns the same fit as a coefficient tibble for use in tables:

``` r

tidy(fs$educ)
#> # A tibble: 16 × 7
#>    term         estimate std.error statistic  p.value  conf.low conf.high
#>    <chr>           <dbl>     <dbl>     <dbl>    <dbl>     <dbl>     <dbl>
#>  1 (Intercept) 16.6        0.241      69.1   0        16.2       17.1    
#>  2 exper       -0.413      0.0337    -12.2   1.17e-33 -0.479     -0.346  
#>  3 expersq      0.000869   0.00165     0.526 5.99e- 1 -0.00237    0.00410
#>  4 black       -0.936      0.0937     -9.98  4.24e-23 -1.12      -0.752  
#>  5 smsa         0.402      0.105       3.84  1.27e- 4  0.197      0.608  
#>  6 south       -0.0516     0.135      -0.381 7.03e- 1 -0.317      0.214  
#>  7 smsa66       0.0255     0.106       0.241 8.10e- 1 -0.182      0.233  
#>  8 reg662      -0.0786     0.187      -0.420 6.74e- 1 -0.446      0.288  
#>  9 reg663      -0.0279     0.183      -0.152 8.79e- 1 -0.387      0.332  
#> 10 reg664       0.117      0.217       0.539 5.90e- 1 -0.309      0.543  
#> 11 reg665      -0.273      0.218      -1.25  2.12e- 1 -0.701      0.156  
#> 12 reg666      -0.303      0.237      -1.28  2.02e- 1 -0.768      0.162  
#> 13 reg667      -0.217      0.234      -0.925 3.55e- 1 -0.676      0.243  
#> 14 reg668       0.524      0.267       1.96  5.02e- 2 -0.000562   1.05   
#> 15 reg669       0.210      0.202       1.04  2.99e- 1 -0.187      0.607  
#> 16 nearc4       0.320      0.0879      3.64  2.76e- 4  0.148      0.492
```

The remaining columns are per-regressor diagnostics that become
informative only with several endogenous regressors: Shea’s (1997)
partial R-squared adjusts for intercorrelation among the instruments,
while the Sanderson-Windmeijer (2016) and Angrist-Pischke (2009)
conditional F statistics assess each endogenous regressor’s
identification given the others. With a single endogenous regressor the
two conditional F statistics collapse to the ordinary first-stage F, and
Shea’s partial R-squared to the ordinary partial R-squared. See
[`?ivreg2`](https://restatr.com/ivreg2r/reference/ivreg2.md) for their
definitions and
[`?first_stage`](https://restatr.com/ivreg2r/reference/first_stage.md)
for the stored per-regressor regressions.

## Robust and clustered inference

### Heteroskedasticity-robust standard errors

Use `vcov = "robust"` for White (1980) heteroskedasticity-robust
standard errors. Add `small = TRUE` for a finite-sample correction
(equivalent to Stata’s `robust small`):

``` r

fit_robust <- ivreg2(iv_formula, data = card, vcov = "robust", small = TRUE)
summary(fit_robust)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = iv_formula, data = card, vcov = "robust", small = TRUE)
#> 
#> Observations: 3,010 
#> VCV type:     Robust, small-sample corrected 
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  3.6661519  0.9109598   4.024 5.85e-05 ***
#> educ         0.1315038  0.0541436   2.429 0.015208 *  
#> exper        0.1082711  0.0234089   4.625 3.90e-06 ***
#> expersq     -0.0023349  0.0003488  -6.695 2.57e-11 ***
#> black       -0.1467758  0.0525019  -2.796 0.005213 ** 
#> smsa         0.1118084  0.0311448   3.590 0.000336 ***
#> south       -0.1446715  0.0291429  -4.964 7.28e-07 ***
#> smsa66       0.0185311  0.0205651   0.901 0.367610    
#> reg662       0.1007677  0.0365488   2.757 0.005867 ** 
#> reg663       0.1482588  0.0355971   4.165 3.20e-05 ***
#> reg664       0.0498971  0.0436162   1.144 0.252714    
#> reg665       0.1462719  0.0492259   2.971 0.002988 ** 
#> reg666       0.1629029  0.0517655   3.147 0.001666 ** 
#> reg667       0.1345722  0.0505567   2.662 0.007814 ** 
#> reg668      -0.0830770  0.0572432  -1.451 0.146801    
#> reg669       0.1078142  0.0410761   2.625 0.008716 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.2382 
#> Adj. R-squared: 0.2343 
#> F(15, 2994):     55.8 (p < 2.2e-16)
#> Root MSE:       0.3883 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(1) = 14.03 (p = 0.0002)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           13.26 
#>   Kleibergen-Paap rk Wald F:     14.14 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       16.38 
#>      15%  maximal IV size       8.96 
#>      20%  maximal IV size       6.66 
#>      25%  maximal IV size       5.53 
#> 
#> Overidentification test (Hansen J):  (equation exactly identified)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(1,2994) = 5.76 (p = 0.0164)
#>   Anderson-Rubin Wald Chi-sq(1) = 5.80 (p = 0.0161)
#>   Stock-Wright LM S Chi-sq(1) = 5.83 (p = 0.0157)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 1.22 (p = 0.2696)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ               14.14    0.0002      0.0044    0.0044     14.14     14.14
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq, black, smsa, south, smsa66, reg662, reg663, reg664, reg665, reg666, reg667, reg668, reg669 
#> Excluded instruments:  nearc4
```

Notice that the diagnostics automatically adapt. The underidentification
test is now the Kleibergen-Paap rk LM statistic instead of the Anderson
LM statistic. In the weak-identification block, the Kleibergen-Paap rk
Wald F (14.14 here) is reported *alongside* the Cragg-Donald Wald F
(13.26), which remains in the output for reference — only the rk
statistic is valid under non-i.i.d. errors. In overidentified models the
Hansen J replaces the Sargan statistic (this fit is exactly identified,
so the line reads “(equation exactly identified)”).

The `small = TRUE` option applies finite-sample corrections throughout:
t-statistics and F-tests replace z-statistics and chi-squared tests.

### Cluster-robust standard errors

When observations are correlated within groups (e.g., repeated
observations on the same individual in panel data), use the `clusters`
argument. Panel data is the natural home of clustering, so we switch
here to the NLS young women’s panel (`nlswork`) and regress log wages on
grade, age, total labor-market experience, and job tenure, clustering on
the individual (`idcode`). The example is ordinary least squares with no
instruments; the `clusters` argument behaves identically once the
formula gains endogenous and excluded-instrument parts, and the
clustered-IV demonstrations live in the other two vignettes.

``` r

data(nlswork)
fit_cluster <- ivreg2(
  ln_wage ~ grade + age + ttl_exp + tenure,
  data = nlswork, clusters = ~ idcode
)
summary(fit_cluster)
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
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.3170 
#> Adj. R-squared: 0.3169 
#> Wald chi2(4):  966.0 (p < 2.2e-16)
#> Root MSE:       0.3949
```

Supplying `clusters` selects the cluster-robust VCE automatically,
without setting `vcov`. To see how much clustering changes the standard
error on `grade`, we compare it against a classical (i.i.d.) fit on the
same data. broom’s
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns a
fit’s coefficient table as a tibble, standard errors and all; stacking
the tidied i.i.d. and clustered fits and keeping the `grade` row puts
the two standard errors side by side:

``` r

fit_iid <- ivreg2(ln_wage ~ grade + age + ttl_exp + tenure, data = nlswork)

grade_se <- bind_rows(iid = tidy(fit_iid), cluster = tidy(fit_cluster), .id = "vce") |>
  filter(term == "grade")
grade_se
#> # A tibble: 2 × 8
#>   vce     term  estimate std.error statistic   p.value conf.low conf.high
#>   <chr>   <chr>    <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
#> 1 iid     grade   0.0744   0.00104      71.5 0           0.0724    0.0765
#> 2 cluster grade   0.0744   0.00216      34.4 8.28e-260   0.0702    0.0787
```

Fitting on the full 28,099 observations, the cluster-robust standard
error on `grade` is 0.0022, versus 0.0010 under the i.i.d. VCE.
Clustering roughly doubles the standard error because wages are serially
correlated within person across the panel’s multiple survey years, a
dependence that classical and heteroskedasticity-robust standard errors
both ignore.

Cluster-robust inference is asymptotic in the number of clusters rather
than observations: with only a handful of clusters the standard errors
can be biased downward even when N is large, and `small = TRUE` only
partially compensates. Here the 4,697 clusters, one per individual, put
us well inside the asymptotic regime. Researchers clustering on a dozen
states or regions should consider wild-cluster bootstrap or other
few-cluster methods.

When `small = TRUE` is added, `ivreg2r` applies the standard cluster
finite-sample correction factor, `(N-1)/(N-K) * M/(M-1)`, where M is the
number of clusters, and reports t- and F-statistics instead of z- and
chi-squared statistics.

For two-way clustering (e.g., by firm and year in a panel), use
`clusters = ~ var1 + var2`. See
[`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md)
for details.

### VCE option summary

| `vcov`     | `clusters` | `small` | Equivalent Stata command         |
|------------|------------|---------|----------------------------------|
| `"iid"`    | —          | `FALSE` | `ivreg2 ...`                     |
| `"iid"`    | —          | `TRUE`  | `ivreg2 ..., small`              |
| `"robust"` | —          | `FALSE` | `ivreg2 ..., robust`             |
| `"robust"` | —          | `TRUE`  | `ivreg2 ..., robust small`       |
| —          | `~ var`    | `FALSE` | `ivreg2 ..., cluster(var)`       |
| —          | `~ var`    | `TRUE`  | `ivreg2 ..., cluster(var) small` |

When `clusters` is supplied, VCE is automatically cluster-robust
regardless of the `vcov` argument.

For HAC, AC, Kiefer, and Driscoll-Kraay VCE options, see
[`vignette("time-series-gmm")`](https://restatr.com/ivreg2r/articles/time-series-gmm.md).

## Weighted estimation

`ivreg2r` supports Stata’s four weight types through `weight_type`:
sampling weights (`"pweight"`), analytic weights (`"aweight"`, the
default), frequency weights (`"fweight"`), and importance weights
(`"iweight"`). All four are passed with the `weights` argument. We work
through two of them — sampling weights on the genuine survey weight the
Card data ship with, and analytic weights on a grouped-data example —
then summarize how the remaining two differ.

### Sampling weights (pweight)

The Card data include a genuine survey weight: `weight` is the NLS
sampling weight, reflecting each respondent’s probability of selection
into the survey. Probability/sampling weights (equivalent to Stata’s
`[pw=varname]`) are specified with `weights` plus
`weight_type = "pweight"`:

``` r

fit_pw <- ivreg2(iv_formula, data = card, weights = weight, weight_type = "pweight",
                 vcov = "robust")
summary(fit_pw)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = iv_formula, data = card, weights = weight, vcov = "robust", 
#>     weight_type = "pweight")
#> 
#> Observations: 3,010 
#> VCV type:     Robust 
#> Weight type:   pweight 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  3.1946899  0.9575119   3.336 0.000849 ***
#> educ         0.1578176  0.0576981   2.735 0.006234 ** 
#> exper        0.1210899  0.0211739   5.719 1.07e-08 ***
#> expersq     -0.0023122  0.0004505  -5.132 2.87e-07 ***
#> black       -0.1461836  0.0505733  -2.891 0.003846 ** 
#> smsa         0.0919528  0.0335047   2.744 0.006061 ** 
#> south       -0.1085573  0.0348286  -3.117 0.001828 ** 
#> smsa66       0.0309925  0.0223835   1.385 0.166171    
#> reg662       0.0972339  0.0373200   2.605 0.009176 ** 
#> reg663       0.1436910  0.0368092   3.904 9.47e-05 ***
#> reg664       0.0378907  0.0463862   0.817 0.414013    
#> reg665       0.1199726  0.0530193   2.263 0.023647 *  
#> reg666       0.1372865  0.0604205   2.272 0.023075 *  
#> reg667       0.1071151  0.0534122   2.005 0.044916 *  
#> reg668      -0.0902232  0.0610340  -1.478 0.139343    
#> reg669       0.0999564  0.0447433   2.234 0.025483 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.1251 
#> Adj. R-squared: 0.1208 
#> Wald chi2(15):  36.1 (p < 2.2e-16)
#> Root MSE:       0.4054 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(1) = 12.85 (p = 0.0003)
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           15.85 
#>   Kleibergen-Paap rk Wald F:     12.96 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       16.38 
#>      15%  maximal IV size       8.96 
#>      20%  maximal IV size       6.66 
#>      25%  maximal IV size       5.53 
#> 
#> Overidentification test (Hansen J):  (equation exactly identified)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(1,2994) = 7.59 (p = 0.0059)
#>   Anderson-Rubin Wald Chi-sq(1) = 7.63 (p = 0.0058)
#>   Stock-Wright LM S Chi-sq(1) = 7.82 (p = 0.0052)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 2.43 (p = 0.1189)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ               12.96    0.0003      0.0053    0.0053     12.96     12.96
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq, black, smsa, south, smsa66, reg662, reg663, reg664, reg665, reg666, reg667, reg668, reg669 
#> Excluded instruments:  nearc4
```

Because pweights describe the sampling design rather than the error
variance, valid standard errors must be robust: `ivreg2r` forces a
robust VCE whenever `weight_type = "pweight"`, regardless of the `vcov`
setting, exactly as Stata does for `[pw=]`. We pass `vcov = "robust"`
explicitly above; leaving `vcov` at its default `"iid"` would trigger
the same forced override, with a warning noting it. The weighted point
estimate answers a population-level question — what is the return to
schooling in the population the NLS sample was drawn from — whereas the
unweighted fit describes the sample.

### Analytic weights (aweight)

An analytic weight is inversely proportional to an observation’s
variance. The textbook case is grouped data: each row is a cell mean,
weighted by the number of individuals that went into the cell, and the
defining property is that a weighted regression on the cell means
reproduces the coefficients of the individual-level regression. We can
watch that identity hold. Collapse the individuals to the cells of three
binary covariates, taking each cell’s mean log wage and its size —
`group_by() |> summarize()` is Stata’s
`collapse (mean) lwage (count) n = lwage, by(black smsa south)`:

``` r

cells <- card |>
  group_by(black, smsa, south) |>
  summarize(lwage = mean(lwage), n = n(), .groups = "drop")

micro   <- ivreg2(lwage ~ black + smsa + south, data = card)
grouped <- ivreg2(lwage ~ black + smsa + south, data = cells,
                  weights = n, weight_type = "aweight")

all.equal(coef(micro), coef(grouped))
#> [1] TRUE
```

The cell means are computed in double precision, so nothing is lost in
the collapse, and the two coefficient vectors agree to numerical
tolerance: `aweight` on the cell means recovers exactly what the
individual-level regression would have given.

### The other weight types

The remaining two types are passed the same way — `weights =` plus
`weight_type =` — and differ only in what the weight *means*:

- **Frequency weights** (`weight_type = "fweight"`, Stata’s `[fw=]`)
  record how many identical observations a row stands for, so a row with
  weight 3 counts as three copies and the sample size becomes the sum of
  the weights:
  `ivreg2(y ~ x, data = d, weights = count, weight_type = "fweight")`.
  They obey an identity of the same kind as analytic weights: an
  fweighted fit reproduces the fit on the fully expanded rows, standard
  errors and sample size included. The [validation
  article](https://restatr.com/ivreg2r/articles/validation.html) checks
  that identity against Stata.
- **Importance weights** (`weight_type = "iweight"`, Stata’s `[iw=]`)
  attach an arbitrary multiplier to each observation without a
  probabilistic interpretation. Like frequency weights they are not
  normalized, so N becomes the sum of the weights, but they are
  restricted to the classical i.i.d. VCE and are incompatible with
  robust, cluster, HAC, and GMM estimation. They exist mainly for
  building weighted computations on top of `ivreg2r`.

See [`?ivreg2`](https://restatr.com/ivreg2r/reference/ivreg2.md) for the
precise rules.

For aweights and pweights the weights are normalized internally to sum
to N, following Stata, so coefficients, standard errors, and test
statistics do not depend on the scale of the weights. One consequence is
that for these two types sigma (the Root MSE) differs from
[`lm()`](https://rdrr.io/r/stats/lm.html)’s by a factor of
`sqrt(N / sum(w))`, because [`lm()`](https://rdrr.io/r/stats/lm.html)
uses the raw weights.

## Beyond 2SLS

This vignette focuses on 2SLS, but `ivreg2r` also supports LIML, Fuller,
and k-class estimation via the `method`, `fuller`, and `kclass`
arguments. LIML has smaller approximate bias than 2SLS when there are
many or weak instruments. Fuller’s (1977) modified LIML further improves
on LIML by having finite moments. These finite-sample advantages are
established under i.i.d. errors; neither estimator is robust to
violations of that assumption (Baum, Schaffer & Stillman, 2007,
pp. 478–479). See
[`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md)
for examples and comparisons.

## Tidyverse workflows

`ivreg2` objects integrate with the tidyverse through
[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html), and
[`augment()`](https://generics.r-lib.org/reference/augment.html).

### Coefficient table

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns a
tibble of coefficient estimates with standard errors, test statistics,
p-values, and confidence intervals:

``` r

tidy(fit_iv)
#> # A tibble: 16 × 7
#>    term        estimate std.error statistic  p.value conf.low conf.high
#>    <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#>  1 (Intercept)  3.67     0.922        3.97  7.05e- 5  1.86      5.47   
#>  2 educ         0.132    0.0548       2.40  1.64e- 2  0.0241    0.239  
#>  3 exper        0.108    0.0236       4.59  4.46e- 6  0.0620    0.155  
#>  4 expersq     -0.00233  0.000333    -7.02  2.22e-12 -0.00299  -0.00168
#>  5 black       -0.147    0.0538      -2.73  6.33e- 3 -0.252    -0.0414 
#>  6 smsa         0.112    0.0316       3.54  3.99e- 4  0.0499    0.174  
#>  7 south       -0.145    0.0272      -5.32  1.06e- 7 -0.198    -0.0913 
#>  8 smsa66       0.0185   0.0216       0.860 3.90e- 1 -0.0237    0.0608 
#>  9 reg662       0.101    0.0376       2.68  7.34e- 3  0.0271    0.174  
#> 10 reg663       0.148    0.0367       4.04  5.39e- 5  0.0763    0.220  
#> 11 reg664       0.0499   0.0436       1.14  2.53e- 1 -0.0356    0.135  
#> 12 reg665       0.146    0.0469       3.12  1.83e- 3  0.0543    0.238  
#> 13 reg666       0.163    0.0518       3.15  1.65e- 3  0.0614    0.264  
#> 14 reg667       0.135    0.0493       2.73  6.31e- 3  0.0380    0.231  
#> 15 reg668      -0.0831   0.0592      -1.40  1.60e- 1 -0.199     0.0329 
#> 16 reg669       0.108    0.0417       2.59  9.73e- 3  0.0261    0.190
```

For the intervals alone, `confint(fit_iv)` returns the same bounds as a
matrix; they are based on the t-distribution when `small = TRUE` was
used at estimation time and on the standard normal otherwise.

### Model summary

[`glance()`](https://generics.r-lib.org/reference/glance.html) returns a
single-row tibble with fit statistics and the headline diagnostic tests
(weak identification, underidentification, overidentification) — a
compact set chosen so that model-comparison tables render sensibly:

``` r

glance(fit_iv)
#> # A tibble: 1 × 15
#>   r.squared adj.r.squared sigma statistic   p.value    df df.residual  nobs
#>       <dbl>         <dbl> <dbl>     <dbl>     <dbl> <int>       <int> <int>
#> 1     0.238         0.234 0.387      51.0 8.58e-136    15        2994  3010
#> # ℹ 7 more variables: vcov_type <chr>, weak_id_stat <dbl>,
#> #   weak_id_robust_stat <dbl>, underid_stat <dbl>, underid_p <dbl>,
#> #   overid_stat <dbl>, overid_p <dbl>
```

### Diagnostic tests as data

[`diagnostics()`](https://restatr.com/ivreg2r/reference/diagnostics.md)
returns the full set of specification tests that
[`summary()`](https://rdrr.io/r/base/summary.html) prints — including
the Stock-Yogo critical values — as a tibble with one row per test:

``` r

diagnostics(fit_iv)
#> # A tibble: 11 × 8
#>    test               test_name statistic    df   df2  p_value tested_vars note 
#>    <chr>              <chr>         <dbl> <int> <int>    <dbl> <chr>       <chr>
#>  1 underid            Anderson…     13.3      1    NA  2.70e-4 NA          NA   
#>  2 weak_id            Cragg-Do…     13.3     NA    NA NA       NA          NA   
#>  3 sy_iv_size_10      Stock-Yo…     16.4     NA    NA NA       NA          NA   
#>  4 sy_iv_size_15      Stock-Yo…      8.96    NA    NA NA       NA          NA   
#>  5 sy_iv_size_20      Stock-Yo…      6.66    NA    NA NA       NA          NA   
#>  6 sy_iv_size_25      Stock-Yo…      5.53    NA    NA NA       NA          NA   
#>  7 overid             Sargan         0        0    NA NA       NA          NA   
#>  8 anderson_rubin_f   Anderson…      5.42     1  2994  2.00e-2 NA          NA   
#>  9 anderson_rubin_ch… Anderson…      5.44     1    NA  1.96e-2 NA          NA   
#> 10 stock_wright       Stock-Wr…      5.43     1    NA  1.97e-2 NA          NA   
#> 11 endogeneity        Endogene…      1.17     1    NA  2.79e-1 educ        NA
```

The `test` column holds stable keys, so pulling one test out for a table
or a threshold comparison is a one-line `filter(test == "weak_id")`.
Tests you request through `endog`, `orthog`, or `redundant` — Stata’s
`endog()`, `orthog()`, and `redundant()` options — appear as additional
rows.

### Augmented data

[`augment()`](https://generics.r-lib.org/reference/augment.html) returns
the original data with fitted values and residuals appended:

``` r

augment(fit_iv) |> head()
#> # A tibble: 6 × 19
#>   lwage exper expersq black  smsa south smsa66 reg662 reg663 reg664 reg665
#>   <dbl> <int>   <int> <int> <int> <int>  <int>  <int>  <int>  <int>  <int>
#> 1  6.31    16     256     1     1     0      1      0      0      0      0
#> 2  6.18     9      81     0     1     0      1      0      0      0      0
#> 3  6.58    16     256     0     1     0      1      0      0      0      0
#> 4  5.52    10     100     0     1     0      1      1      0      0      0
#> 5  6.59    16     256     0     1     0      1      1      0      0      0
#> 6  6.21     8      64     0     1     0      1      1      0      0      0
#> # ℹ 8 more variables: reg666 <int>, reg667 <int>, reg668 <int>, reg669 <int>,
#> #   educ <int>, nearc4 <int>, .fitted <dbl>, .resid <dbl>
```

## Stata migration quick reference

#### Estimation basics and accessors

| Task | Stata | R (`ivreg2r`) |
|----|----|----|
| Basic IV | `ivreg2 y x (endo = z)` | `ivreg2(y ~ x | endo | z, data = d)` |
| Extra instruments, no endogenous regressor | `ivreg2 y x (=z1 z2)` | `ivreg2(y ~ x | 0 | z1 + z2, data = d)` (the `0` marks an empty endogenous group) |
| Lag / difference | `L.x` / `D.x` | `l(x, 1)` / `d(x, 1)` (requires `tvar`, plus `ivar` for panels) |
| Coefficients | `e(b)` | `coef(fit)` |
| VCV matrix | `e(V)` | `vcov(fit)` |
| Residuals | `predict resid, resid` | `residuals(fit)` |
| Fitted values | `predict yhat` | `fitted(fit)` |
| Prediction SE | `predict s, stdp` | `predict(fit, se.fit = TRUE)` |
| Confidence intervals | — | `confint(fit)` |
| Tidy coef table | — | `tidy(fit)` |
| Summary statistics | `ereturn list` | `glance(fit)` |
| Diagnostics | Displayed after estimation | Displayed by `summary(fit)`; as a tibble via `diagnostics(fit)` |

#### Estimators

| Task | Stata | R (`ivreg2r`) |
|----|----|----|
| LIML | `ivreg2 ..., liml` | `ivreg2(..., method = "liml")` |
| Fuller(1) | `ivreg2 ..., fuller(1)` | `ivreg2(..., fuller = 1)` |
| k-class | `ivreg2 ..., kclass(0.5)` | `ivreg2(..., kclass = 0.5)` |
| COVIV | `ivreg2 ..., liml coviv` | `ivreg2(..., method = "liml", coviv = TRUE)` |
| GMM2S | `ivreg2 ..., gmm2s` | `ivreg2(..., method = "gmm2s")` |
| CUE | `ivreg2 ..., cue` | `ivreg2(..., method = "cue")` |
| Evaluate CUE at fixed b | `ivreg2 ..., b0(bmat)` | `ivreg2(..., b0 = c(...))` |
| User W matrix | `ivreg2 ..., wmatrix(name)` | `ivreg2(..., wmatrix = W)` |
| User S matrix | `ivreg2 ..., smatrix(name)` | `ivreg2(..., smatrix = S)` |
| Moment covariance / weight matrix | `e(S)` / `e(W)` | `fit$S` / `fit$W` |

#### Variance estimation and clustering

| Task | Stata | R (`ivreg2r`) |
|----|----|----|
| Robust SEs | `ivreg2 ..., robust` | `ivreg2(..., vcov = "robust")` |
| Robust + small | `ivreg2 ..., robust small` | `ivreg2(..., vcov = "robust", small = TRUE)` |
| Cluster SEs | `ivreg2 ..., cluster(g)` | `ivreg2(..., clusters = ~ g)` |
| Cluster + small | `ivreg2 ..., cluster(g) small` | `ivreg2(..., clusters = ~ g, small = TRUE)` |
| Two-way cluster | `ivreg2 ..., cluster(g1 g2)` | `ivreg2(..., clusters = ~ g1 + g2)` |
| AC kernel | `ivreg2 ..., kernel(bar) bw(3)` | `ivreg2(..., kernel = "bartlett", bw = 3, tvar = "year")` |
| HAC kernel | `ivreg2 ..., robust kernel(bar) bw(3)` | `ivreg2(..., vcov = "robust", kernel = "bartlett", bw = 3, tvar = "year")` |
| Auto bandwidth | `ivreg2 ..., kernel(bar) bw(auto)` | `ivreg2(..., kernel = "bartlett", bw = "auto", tvar = "year")` |
| Driscoll-Kraay | `ivreg2 ..., dkraay(3)` | `ivreg2(..., dkraay = 3, tvar = "year", ivar = "id")` |
| Kiefer VCE | `ivreg2 ..., kiefer` | `ivreg2(..., kiefer = TRUE, tvar = "year", ivar = "id")` |
| Stock-Watson VCE | `ivreg2 ..., sw` | `ivreg2(..., sw = TRUE, ivar = "id")` |
| Center moments | `ivreg2 ..., center` | `ivreg2(..., center = TRUE)` |

#### Weights, tests, and other options

| Task | Stata | R (`ivreg2r`) |
|----|----|----|
| Analytic weights | `ivreg2 ... [aw=w]` | `ivreg2(..., weights = w)` |
| Frequency weights | `ivreg2 ... [fw=w]` | `ivreg2(..., weights = w, weight_type = "fweight")` |
| Probability weights | `ivreg2 ... [pw=w]` | `ivreg2(..., weights = w, weight_type = "pweight")` |
| Importance weights | `ivreg2 ... [iw=w]` | `ivreg2(..., weights = w, weight_type = "iweight")` |
| Orthogonality test | `ivreg2 ..., orthog(z1)` | `ivreg2(..., orthog = c("z1"))` |
| Endogeneity test | `ivreg2 ..., endog(z1)` | `ivreg2(..., endog = c("z1"))` |
| Redundancy test | `ivreg2 ..., redundant(z1)` | `ivreg2(..., redundant = c("z1"))` |
| DoF adjustment | `ivreg2 ..., dofminus(5)` | `ivreg2(..., dofminus = 5)` |
| Small-sample DoF | `ivreg2 ..., sdofminus(3)` | `ivreg2(..., sdofminus = 3)` |
| Partialling | `ivreg2 ..., partial(x1)` | `ivreg2(..., partial = c("x1"))` |
| Partial constant / all | `ivreg2 ..., partial(_cons)` | `ivreg2(..., partial = "_cons")` or `"_all"` |
| No partial DoF adjustment | `ivreg2 ..., nopartialsmall` | `ivreg2(..., nopartialsmall = TRUE)` |
| Save first stage | `ivreg2 ..., savefirst` | `ivreg2(..., first_stage = TRUE)`, then `first_stage(fit)` |
| Save reduced form | `ivreg2 ..., saverf` | `ivreg2(..., reduced_form = "rf")` (`rf` displays, `saverf` stores) |
| Save first-stage system | `ivreg2 ..., savesfirst` | `ivreg2(..., reduced_form = "system")` (`savesfirst` is undocumented in Stata) |
| Suppress ID stats | `ivreg2 ..., noid` | `ivreg2(..., noid = TRUE)` |

See
[`vignette("advanced-iv")`](https://restatr.com/ivreg2r/articles/advanced-iv.md)
for LIML, k-class, two-way clustering, and model comparison examples.
See
[`vignette("time-series-gmm")`](https://restatr.com/ivreg2r/articles/time-series-gmm.md)
for HAC/AC, GMM, CUE, partialling, and panel VCE examples.

### Intentional differences from Stata

| Area | Stata `ivreg2` | `ivreg2r` | Rationale |
|----|----|----|----|
| Formula | `ivreg2 y x (endo = z)` | `y ~ x | endo | z` | R convention; each variable listed once |
| Numerical method | Cross-products via Mata | QR-based regression via `lm.fit` | Better conditioning for core estimation |
| Diagnostics | Some computed at post-estimation | All computed at estimation time | Returned as a tibble by [`diagnostics()`](https://restatr.com/ivreg2r/reference/diagnostics.md); headline tests also in [`glance()`](https://generics.r-lib.org/reference/glance.html) |
| Factor variables | `i.` prefix syntax | R [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) / contrasts | R convention |
| Time-series | `tsset` declares time/panel | `tvar=` / `ivar=` arguments | R has no `tsset`; explicit args |
| CUE optimizer | Mata [`optimize()`](https://rdrr.io/r/stats/optimize.html) | R [`optim()`](https://rdrr.io/r/stats/optim.html) (BFGS + Nelder-Mead) | Different optimizer, same objective |
| `redundant` with no endogenous regressors | Silently ignored | Warns that the test is skipped | An explicitly requested test should not disappear silently |
| Missing values in a cluster variable | Silently drops the affected observations from the estimation sample | Raises an error asking the user to drop or handle those rows explicitly | The estimation sample should never change silently |

## References

- Anderson, T.W. (1951). “Estimating Linear Restrictions on Regression
  Coefficients for Multivariate Normal Distributions.” *Annals of
  Mathematical Statistics*, 22(3), 327–351.
- Anderson, T.W. and Rubin, H. (1949). “Estimation of the Parameters of
  a Single Equation in a Complete System of Stochastic Equations.”
  *Annals of Mathematical Statistics*, 20(1), 46–63.
- Angrist, J.D. and Pischke, J.-S. (2009). *Mostly Harmless
  Econometrics*. Princeton University Press.
- Baum, C.F., Schaffer, M.E., and Stillman, S. (2003). “Instrumental
  Variables and GMM: Estimation and Testing.” *Stata Journal*, 3(1),
  1–31.
- Baum, C.F., Schaffer, M.E., and Stillman, S. (2007). “Enhanced
  Routines for Instrumental Variables/Generalized Method of Moments
  Estimation and Testing.” *Stata Journal*, 7(4), 465–506.
- Card, D. (1995). “Using Geographic Variation in College Proximity to
  Estimate the Return to Schooling.” In L.N. Christofides, E.K. Grant,
  and R. Swidinsky (Eds.), *Aspects of Labour Market Behaviour: Essays
  in Honour of John Vanderkamp*. University of Toronto Press.
- Cragg, J.G. and Donald, S.G. (1993). “Testing Identifiability and
  Specification in Instrumental Variable Models.” *Econometric Theory*,
  9(2), 222–240.
- Fuller, W.A. (1977). “Some Properties of a Modification of the Limited
  Information Estimator.” *Econometrica*, 45(4), 939–953.
- Guggenberger, P. (2010). “The Impact of a Hausman Pretest on the
  Asymptotic Size of a Hypothesis Test.” *Econometric Theory*, 26(2),
  369–382.
- Hansen, L.P. (1982). “Large Sample Properties of Generalized Method of
  Moments Estimators.” *Econometrica*, 50(4), 1029–1054.
- Hayashi, F. (2000). *Econometrics*. Princeton University Press.
- Kleibergen, F. and Paap, R. (2006). “Generalized Reduced Rank Tests
  Using the Singular Value Decomposition.” *Journal of Econometrics*,
  133(1), 97–126.
- Sanderson, E. and Windmeijer, F. (2016). “A Weak Instrument F-Test in
  Linear IV Models with Multiple Endogenous Variables.” *Journal of
  Econometrics*, 190(2), 212–221.
- Sargan, J.D. (1958). “The Estimation of Economic Relationships Using
  Instrumental Variables.” *Econometrica*, 26(3), 393–415.
- Shea, J. (1997). “Instrument Relevance in Multivariate Linear Models:
  A Simple Measure.” *Review of Economics and Statistics*, 79(2),
  348–352.
- Staiger, D. and Stock, J.H. (1997). “Instrumental Variables Regression
  with Weak Instruments.” *Econometrica*, 65(3), 557–586.
- Stock, J.H. and Wright, J.H. (2000). “GMM with Weak Identification.”
  *Econometrica*, 68(5), 1055–1096.
- Stock, J.H. and Yogo, M. (2005). “Testing for Weak Instruments in
  Linear IV Regression.” In D.W.K. Andrews and J.H. Stock (Eds.),
  *Identification and Inference for Econometric Models: Essays in Honor
  of Thomas Rothenberg*. Cambridge University Press.
- White, H. (1980). “A Heteroskedasticity-Consistent Covariance Matrix
  Estimator and a Direct Test for Heteroskedasticity.” *Econometrica*,
  48(4), 817–838.
- Wooldridge, J.M. (2020). *Introductory Econometrics: A Modern
  Approach*, 7th edition. Cengage.
