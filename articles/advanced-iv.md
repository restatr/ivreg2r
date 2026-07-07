# Advanced IV Estimation with ivreg2r

## Introduction

This vignette covers features of `ivreg2r` beyond the basic two-stage
least squares workflow of
[`vignette("introduction")`](https://restatr.com/ivreg2r/articles/introduction.md):
alternative estimators, the instrument-validity tests,
weak-instrument-robust inference, two-way clustering, and the reduced
form. It is organized around the canonical example groups — a Mroz
labor-supply equation that carries most of the help file’s testing arc,
the Griliches wage equation for the help file’s orthogonality example
and for Baum, Schaffer & Stillman’s (2007) weak-instrument and
redundancy illustrations, and Klein’s consumption function that
separates the estimators — and adds a two-way clustered IV example on a
cigarette-demand panel.

The topics are:

- LIML, Fuller, and k-class estimation, including the canonical Klein
  consumption function
- Instrument-validity tests: overidentification, the orthogonality
  C-statistic, endogeneity, and redundancy
- Weak-instrument-robust inference (the Anderson-Rubin and Stock-Wright
  statistics)
- Two-way clustering, including with instrumental variables
- Reduced-form regression and degrees-of-freedom adjustments

For HAC and AC standard errors, the
[`l()`](https://restatr.com/ivreg2r/reference/ts-operators.md)/[`d()`](https://restatr.com/ivreg2r/reference/ts-operators.md)
time-series operators, GMM and the continuously-updated estimator, the
panel variance estimators (Kiefer, Driscoll-Kraay, and
cluster-plus-kernel), partialling out regressors, and user-supplied
weighting matrices, see
[`vignette("time-series-gmm")`](https://restatr.com/ivreg2r/articles/time-series-gmm.md).
The four weight types are documented in
[`vignette("introduction")`](https://restatr.com/ivreg2r/articles/introduction.md)
and in [`?ivreg2`](https://restatr.com/ivreg2r/reference/ivreg2.md).

Five bundled datasets appear below: the Mroz (1987) female labor-supply
data (`mroz`), Klein’s (1950) Model I macroeconomic series (`klein`),
the Griliches (1976) young-men’s wage data (`griliches`), the NLS
young-women panel (`nlswork`, distributed with Stata’s `[XT]` manual),
and the Baltagi cigarette-demand panel (`cigar`).

``` r

library(ivreg2r)
data(mroz)
data(klein)
data(griliches)
data(nlswork)
```

## The Mroz baseline

The Mroz (1987) dataset contains 753 married women from the 1975 Panel
Study of Income Dynamics, of whom 428 were working (`inlf == 1`) and
have observed wages. Years of education (`educ`) is potentially
endogenous in a wage equation because of ability bias. Following the
illustrative specification used throughout Stata’s `ivreg2` help file
and in Baum, Schaffer & Stillman (2007, pp. 471–472), the excluded
instruments are `age` and the numbers of children at home under 6
(`kidslt6`) and 6–18 (`kidsge6`). These exclusion restrictions are
illustrative and contestable — children plausibly affect wage offers
through labor-supply channels — and the help file adopts them precisely
because they make the instrument-testing machinery interesting.

``` r

mroz_work <- subset(mroz, inlf == 1)
nrow(mroz_work)
#> [1] 428
```

The baseline regresses log wages on experience and its square, with
education instrumented by the three excluded instruments. Three excluded
instruments for one endogenous regressor leave two overidentifying
restrictions, which is what makes the model useful for comparing
estimators and testing instruments:

``` r

fit_2sls <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                   data = mroz_work)
summary(fit_2sls)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq | educ | age + kidslt6 + 
#>     kidsge6, data = mroz_work)
#> 
#> Observations: 428 
#> VCV type:     Classical (iid) 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)   
#> (Intercept) -0.3848718  1.0115512  -0.380  0.70359   
#> educ         0.0964002  0.0814278   1.184  0.23646   
#> exper        0.0421930  0.0138831   3.039  0.00237 **
#> expersq     -0.0008323  0.0004204  -1.980  0.04773 * 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.1556 
#> Adj. R-squared: 0.1496 
#> Wald chi2(3):  7.5 (p = 0.0001)
#> Root MSE:       0.6638 
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           4.34 
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
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(3) = 12.82 (p = 0.0051)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(3,422) = 0.61 (p = 0.6076)
#>   Anderson-Rubin Wald Chi-sq(3) = 1.86 (p = 0.6016)
#>   Stock-Wright LM S Chi-sq(3) = 1.85 (p = 0.6033)
#> 
#> Overidentification test (Sargan):
#>   Chi-sq(2) = 0.70 (p = 0.7042)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 0.02 (p = 0.8899)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ                4.34    0.0050      0.0299    0.0299      4.34      4.34
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq 
#> Excluded instruments:  age, kidslt6, kidsge6
```

## LIML, Fuller, and k-class estimation

### LIML

Limited-information maximum likelihood (LIML) is an alternative to 2SLS
that behaves better under weak instruments: it is approximately
median-unbiased and its Wald test suffers less size distortion (Stock &
Yogo, 2005). LIML has no finite moments of any order, whatever the
degree of overidentification (Mariano & Sawa, 1972), so individual LIML
estimates can be extreme in small samples; Fuller’s modification,
introduced below, restores finite moments. Request it with
`method = "liml"`:

``` r

fit_liml <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                   data = mroz_work, method = "liml")
summary(fit_liml)
#> 
#> LIML Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq | educ | age + kidslt6 + 
#>     kidsge6, data = mroz_work, method = "liml")
#> 
#> Observations: 428 
#> VCV type:     Classical (iid) 
#> lambda:       1.001642 
#> kclass:       1.001642 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)   
#> (Intercept) -0.3769294  1.0394246  -0.363  0.71688   
#> educ         0.0957581  0.0836906   1.144  0.25254   
#> exper        0.0422292  0.0139270   3.032  0.00243 **
#> expersq     -0.0008335  0.0004220  -1.975  0.04827 * 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.1555 
#> Adj. R-squared: 0.1495 
#> Wald chi2(3):  7.5 (p = 0.0001)
#> Root MSE:       0.6638 
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           4.34 
#>   Stock-Yogo critical values (LIML size):
#>      10%  maximal LIML size      6.46 
#>      15%  maximal LIML size      4.36 
#>      20%  maximal LIML size      3.69 
#>      25%  maximal LIML size      3.32 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(3) = 12.82 (p = 0.0051)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(3,422) = 0.61 (p = 0.6076)
#>   Anderson-Rubin Wald Chi-sq(3) = 1.86 (p = 0.6016)
#>   Stock-Wright LM S Chi-sq(3) = 1.85 (p = 0.6033)
#> 
#> Overidentification test (Sargan):
#>   Chi-sq(2) = 0.70 (p = 0.7042)
#> 
#> Anderson-Rubin overidentification:
#>   LR Chi-sq(2) = 0.702 (p = 0.7040)
#>   Linearized Chi-sq(2) = 0.703 (p = 0.7038)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 0.02 (p = 0.8899)
#>   Tested: educ
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   educ                4.34    0.0050      0.0299    0.0299      4.34      4.34
#> 
#> Instrumented:          educ 
#> Included instruments:  exper, expersq 
#> Excluded instruments:  age, kidslt6, kidsge6
```

The output reports `lambda`, the LIML eigenvalue. When `lambda` is close
to 1, LIML and 2SLS give similar results, since 2SLS is the k-class
member with `k = 1`. A `lambda` far above 1 is evidence against the
overidentifying restrictions, not of weak instruments: `N * log(lambda)`
is the Anderson-Rubin (1950) overidentification LR statistic reported in
the Anderson-Rubin section below, so judge how far above 1 is far by
that statistic rather than by eye. Under weak but valid instruments
`lambda` stays near 1, and weak identification is diagnosed by the
Cragg-Donald and Kleibergen-Paap statistics instead. Here `lambda` is
1.0016, so LIML barely moves the education coefficient away from 2SLS.

### Fuller modification

Fuller’s (1977) modified LIML uses `k = lambda - alpha / (N - L)`, where
`L` is the number of instruments including the exogenous regressors and
`alpha` is the Fuller parameter. The standard choice `fuller = 1` is
recommended by Fuller (1977) and by Davidson & MacKinnon (1993,
pp. 649–650), as cited in Baum, Schaffer & Stillman (2007, p. 479):

``` r

fit_fuller1 <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                      data = mroz_work, fuller = 1)
c(k = fit_fuller1$kclass_value, educ = unname(coef(fit_fuller1)["educ"]))
#>          k       educ 
#> 0.99927193 0.09666366
```

### k-class estimation

The k-class family nests OLS (`k = 0`), 2SLS (`k = 1`), and LIML
(`k = lambda`). The `kclass` argument accepts an arbitrary `k`; the
documented non-trivial choice is Nagar’s (1959) bias-adjusted estimator
with `k = 1 + (L - K) / N`, where `L - K` is the number of
overidentifying restrictions. Rather than hand-counting `L` and `K`,
derive them from the fitted baseline — `rankzz` is the rank of the full
instrument matrix and `rank` is the number of estimated coefficients:

``` r

N <- nobs(fit_2sls)
L <- fit_2sls$rankzz
K <- fit_2sls$rank
k_nagar <- 1 + (L - K) / N
fit_nagar <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                    data = mroz_work, kclass = k_nagar)
c(L = L, K = K, k = k_nagar, educ = unname(coef(fit_nagar)["educ"]))
#>          L          K          k       educ 
#> 6.00000000 4.00000000 1.00467290 0.09436094
```

Both the Fuller and Nagar estimators tend to have better finite-sample
behavior than 2SLS under weak instruments, though neither is robust to
violations of the i.i.d. assumption (Baum, Schaffer & Stillman, 2007,
pp. 478–479).

### COVIV

With LIML or k-class under robust or cluster-robust standard errors, the
sandwich bread uses the k-class projection matrix by default. Setting
`coviv = TRUE` computes the covariance from the 2SLS projection
(`k = 1`) instead, keeping the LIML residuals in the meat. Under i.i.d.
errors, LIML with COVIV coincides with the continuously-updated
estimator — the demonstration appears in the Klein suite below:

``` r

fit_coviv <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                    data = mroz_work, method = "liml", vcov = "robust",
                    small = TRUE, coviv = TRUE)
all.equal(coef(fit_coviv), coef(fit_liml))
#> [1] TRUE
subset(tidy(fit_coviv), term == "educ", c(term, estimate, std.error))
#> # A tibble: 1 × 3
#>   term  estimate std.error
#>   <chr>    <dbl>     <dbl>
#> 1 educ    0.0958    0.0869
```

The asserted identity makes the division of labor explicit: COVIV
changes only the covariance, so the coefficients are exactly the LIML
coefficients, while the displayed standard error comes from the
2SLS-bread sandwich.

### Choosing among the estimators

The case for LIML and its relatives is a many-instruments argument.
Two-stage least squares is biased toward OLS, and that bias grows with
the number of overidentifying restrictions relative to the strength of
the instruments; LIML removes its leading term, which is why LIML
tolerates weaker instruments before its Wald test becomes unreliable
(Stock & Yogo, 2005). When the model is just-identified the argument is
moot: with one excluded instrument per endogenous regressor, LIML,
Fuller, and 2SLS coincide numerically, and the advantage is entirely a
feature of overidentified models. Under genuinely weak instruments
LIML’s lack of finite moments cuts the other way, since individual
estimates can be far from the truth even while the estimator stays
median-unbiased. The ranking of LIML over 2SLS is itself established
within the i.i.d. Stock-Yogo framework (Baum, Schaffer & Stillman, 2007,
pp. 489–490), and switching estimator because a pretest flagged weak
instruments carries the usual consequences of pre-testing for the
inference that follows.

On the Mroz baseline these distinctions barely register, because
`lambda` is essentially 1 and all four estimators return nearly the same
education coefficient. To see the estimators genuinely diverge we need a
small, heavily overidentified problem, which is exactly what the help
file’s Klein example provides.

### The Klein consumption function

Klein’s (1950) Model I consumption function regresses consumption on
lagged profits and treats current profits and the total wage bill as
endogenous, with six excluded instruments: government spending, indirect
taxes plus net exports, a linear trend, the government wage bill, lagged
capital, and lagged total income. This is the canonical LIML example
from Stata’s `ivreg2` help file. Two of the terms are lags built on the
fly with the
[`l()`](https://restatr.com/ivreg2r/reference/ts-operators.md)
time-series operator against the year variable `yr`; the operator is
documented in
[`vignette("time-series-gmm")`](https://restatr.com/ivreg2r/articles/time-series-gmm.md).

``` r

klein_form <- consump ~ l(profits, 1) | profits + wagetot |
  govt + taxnetx + year + wagegovt + capital1 + l(totinc, 1)
fit_k_liml <- ivreg2(klein_form, data = klein, tvar = "yr", method = "liml")
summary(fit_k_liml)
#> 
#> LIML Estimation
#> 
#> Call:
#> ivreg2(formula = klein_form, data = klein, method = "liml", tvar = "yr")
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
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.9566 
#> Adj. R-squared: 0.9489 
#> Wald chi2(3):  118.4 (p = 0.0000)
#> Root MSE:       1.3953 
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           2.89 
#>   Stock-Yogo critical values (LIML size):
#>      10%  maximal LIML size      4.06 
#>      15%  maximal LIML size      2.95 
#>      20%  maximal LIML size      2.63 
#>      25%  maximal LIML size      2.46 
#> 
#> Underidentification test (Anderson canon. corr. LM statistic):
#>   Chi-sq(5) = 12.01 (p = 0.0347)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(6,13) = 18.01 (p = 0.0000)
#>   Anderson-Rubin Wald Chi-sq(6) = 174.61 (p < 2.2e-16)
#>   Stock-Wright LM S Chi-sq(6) = 18.75 (p = 0.0046)
#> 
#> Overidentification test (Sargan):
#>   Chi-sq(4) = 6.99 (p = 0.1365)
#> 
#> Anderson-Rubin overidentification:
#>   LR Chi-sq(4) = 8.497 (p = 0.0750)
#>   Linearized Chi-sq(4) = 10.474 (p = 0.0332)
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
```

With two endogenous regressors and six excluded instruments the model
has four overidentifying restrictions on only 21 observations, and here
`lambda` is 1.4987, well above the near-1 value on Mroz. That is the
regime in which the choice of estimator matters. Fitting 2SLS, LIML,
Fuller(1), and Nagar’s k-class on the same equation shows the two
endogenous coefficients moving with the estimator:

``` r

fit_k_2sls   <- ivreg2(klein_form, data = klein, tvar = "yr")
fit_k_fuller <- ivreg2(klein_form, data = klein, tvar = "yr",
                       method = "liml", fuller = 1)
# Nagar's k = 1 + (L - K)/N = 1 + 4/21, which the help file rounds to 1.19
fit_k_nagar  <- ivreg2(klein_form, data = klein, tvar = "yr",
                       method = "kclass", kclass = 1.19)
klein_models <- list("2SLS" = fit_k_2sls, LIML = fit_k_liml,
                     "Fuller(1)" = fit_k_fuller, Nagar = fit_k_nagar)
round(sapply(klein_models, \(m) coef(m)[c("profits", "wagetot")]), 4)
#>           2SLS    LIML Fuller(1)   Nagar
#> profits 0.0173 -0.2225   -0.1686 -0.0497
#> wagetot 0.8102  0.8226    0.8201  0.8140
```

The estimators genuinely diverge here: the `profits` coefficient changes
sign between 2SLS and LIML, and the wage-bill coefficient moves in the
second decimal place — differences that are invisible on the strongly
identified Mroz equation, and the reason the help file uses Klein to
introduce these estimators.

### LIML with COVIV equals CUE

The Klein example also demonstrates an equivalence. Under
conditional-homoskedasticity (i.i.d.) weighting, LIML with COVIV and the
continuously-updated estimator (CUE) are mathematically the same
estimator, so they should return the same coefficients. CUE reaches its
answer by numerical optimization, so the two agree up to optimizer
tolerance rather than exactly:

``` r

fit_k_coviv <- ivreg2(klein_form, data = klein, tvar = "yr",
                      method = "liml", coviv = TRUE)
fit_k_cue   <- ivreg2(klein_form, data = klein, tvar = "yr", method = "cue")
max(abs(coef(fit_k_coviv) - coef(fit_k_cue)))
#> [1] 1.077372e-06
```

On this tiny (N = 21) problem the cross-platform agreement is on the
order of 1e-5; the residual gap is the optimizer’s endpoint, not a
difference in the estimators. CUE generalizes LIML with COVIV to the
non-i.i.d. case, where a closed-form eigenvalue solution no longer
exists.

## Instrument-validity tests

The tests below follow the escalating order of the help file’s testing
examples: a joint overidentification test, then a test of a subset of
instruments (the help file’s Griliches example), then a test of a single
regressor’s endogeneity, then the empty-endogenous form for an included
regressor, and finally the LIML overidentification statistics. Every
automatic test lives in `fit$diagnostics`, a named list whose elements
[`glance()`](https://generics.r-lib.org/reference/glance.html) exposes
as columns; the chunks below use whichever access route reads better.

### Overidentification

With three excluded instruments and one endogenous regressor, the
baseline has two overidentifying restrictions. The Sargan (1958) test is
a joint test of those restrictions: a rejection signals that some
instrument violates the orthogonality conditions but cannot attribute
the failure to a specific instrument.

``` r

glance(fit_2sls)[, c("overid_stat", "overid_p")]
#> # A tibble: 1 × 2
#>   overid_stat overid_p
#>         <dbl>    <dbl>
#> 1       0.702    0.704
```

### Orthogonality (the C-statistic)

The `orthog` argument tests whether a *subset* of instruments is
exogenous, conditional on the rest being valid — the
difference-in-Sargan C-statistic (Hayashi, 2000). Stata’s help file
demonstrates it on the Griliches (1976) wage equation, challenging the
exogeneity of `age` and `mrt` (marital status) while maintaining `med`
(mother’s education) and `kww` (a test score), estimated by two-step
GMM:

``` r

fit_orthog <- ivreg2(
  lw ~ s + expr + tenure + rns + smsa + factor(year) | iq |
    med + kww + age + mrt,
  data = griliches, method = "gmm2s", orthog = c("age", "mrt")
)
glance(fit_orthog)[, c("orthog_stat", "orthog_p")]
#> # A tibble: 1 × 2
#>   orthog_stat orthog_p
#>         <dbl>    <dbl>
#> 1        86.6 1.55e-19
```

The null is that `age` and `mrt` satisfy the exclusion restrictions
given that the remaining instruments do. The C-statistic rejects
emphatically, so the data are inconsistent with treating `age` and `mrt`
as valid exclusions in this equation. The same pair reappears below as
Baum, Schaffer & Stillman’s weak-instruments illustration. The equation
must remain identified after the challenged instruments are dropped, and
under robust or cluster-robust standard errors the test is a difference
in Hansen J statistics rather than in Sargan statistics.

### Endogeneity of a regressor

The `endog` argument tests whether a nominated endogenous regressor is
in fact exogenous, comparing the model that instruments it against the
model that treats it as exogenous. On the baseline this reproduces
Stata’s `endog(educ)` example:

``` r

fit_endog <- ivreg2(lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
                    data = mroz_work, endog = "educ")
glance(fit_endog)[, c("endogeneity_stat", "endogeneity_p")]
#> # A tibble: 1 × 2
#>   endogeneity_stat endogeneity_p
#>              <dbl>         <dbl>
#> 1           0.0191         0.890
```

With a single endogenous regressor this coincides numerically with the
endogeneity test that
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md) already
reports for the baseline, so `endog` earns its keep only with two or
more endogenous regressors, where it can test a subset:

``` r

c(automatic = glance(fit_2sls)$endogeneity_stat,
  endog_arg = glance(fit_endog)$endogeneity_stat)
#>  automatic  endog_arg 
#> 0.01914714 0.01914714
```

### Testing an included regressor: the empty-endogenous form

The C-statistic can also challenge an *included* regressor’s exogeneity.
Write the model with no endogenous regressors — the empty-endogenous
form `y ~ exog | 0 | instruments`, Stata’s `(= z1 z2)` — so the equation
is estimated by OLS while the excluded instruments supply surplus
orthogonality conditions, then test the regressor with `orthog`. This is
the help file’s C-test of `educ` in the Mroz equation:

``` r

fit_hols <- ivreg2(
  lwage ~ exper + expersq + educ | 0 | age + kidslt6 + kidsge6,
  data = mroz, orthog = "educ"
)
summary(fit_hols)
#> 
#> OLS Estimation
#> 
#> Call:
#> ivreg2(formula = lwage ~ exper + expersq + educ | 0 | age + kidslt6 + 
#>     kidsge6, data = mroz, orthog = "educ")
#> 
#> Observations: 428 
#> VCV type:     Classical (iid) 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.5220407  0.1977017  -2.641  0.00828 ** 
#> exper        0.0415665  0.0131135   3.170  0.00153 ** 
#> expersq     -0.0008112  0.0003914  -2.073  0.03822 *  
#> educ         0.1074896  0.0140802   7.634 2.27e-14 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.1568 
#> Adj. R-squared: 0.1509 
#> Wald chi2(3):  26.3 (p = 0.0000)
#> Root MSE:       0.6633 
#> 
#> 
#> Overidentification test (Sargan):
#>   Chi-sq(3) = 0.72 (p = 0.8681)
#> 
#> Orthogonality test (C-stat):
#>   Chi-sq(1) = 0.02 (p = 0.8899)
#>   Tested: educ
#> 
#> Included instruments:  exper, expersq, educ 
#> Excluded instruments:  age, kidslt6, kidsge6
```

Estimation is exactly OLS — the instruments play no role in the
coefficients — but the Sargan statistic now tests all surplus
conditions, and the C-statistic isolates `educ` by comparing this model
against one in which `educ` is instrumented by `age`, `kidslt6`, and
`kidsge6`. A small C-statistic means the data do not reject treating
`educ` as exogenous, conditional on the three excluded instruments being
valid. The same form with `method = "gmm2s"` and `vcov = "robust"` gives
Cragg’s (1983) heteroskedastic OLS (HOLS) estimator, which exploits the
surplus orthogonality conditions for efficiency under
heteroskedasticity. Under the default i.i.d. weighting there is no
efficiency gain to exploit: the two-step weighting matrix is
proportional to `(Z'Z)^{-1}` and the estimates are identical to OLS, in
Stata just as here. A cluster or HAC weighting generalizes the
efficiency gain to those error structures.

### Anderson-Rubin LIML overidentification

Under LIML with i.i.d. errors the summary reports two further
overidentification statistics. The Anderson-Rubin (1950) LR statistic is
`N * log(lambda)`, and its linearized form is `N * (lambda - 1)`; both
are chi-squared with `L - K` degrees of freedom under the null that the
overidentifying restrictions hold. When `lambda` is close to 1 the three
overidentification tests — Sargan, AR LR, and AR linearized — agree
closely:

``` r

glance(fit_liml)[, c("ar_overid_lr_stat", "ar_overid_lr_p",
                      "ar_overid_lin_stat", "ar_overid_lin_p")]
#> # A tibble: 1 × 4
#>   ar_overid_lr_stat ar_overid_lr_p ar_overid_lin_stat ar_overid_lin_p
#>               <dbl>          <dbl>              <dbl>           <dbl>
#> 1             0.702          0.704              0.703           0.704
```

## Weak-instrument-robust inference

### The Stock-Wright S statistic

The Stock-Wright (2000) S statistic is a weak-instrument-robust test of
the joint null that the coefficients on the endogenous regressors are
zero and the overidentifying restrictions hold. The Anderson-Rubin
statistic tests the same null, and under it both are chi-squared with
`L1` degrees of freedom, where `L1` is the number of excluded
instruments: Anderson-Rubin is the Wald form of the test and S is the LM
/ GMM-distance form (Baum, Schaffer & Stillman, 2007, pp. 491–492).
Unlike Wald tests on 2SLS coefficients, both keep correct size even when
instruments are weak. The S statistic is the value of the CUE objective
with the exogenous regressors partialled out, and it is computed for
every IV model:

``` r

glance(fit_2sls)[, c("stock_wright_stat", "stock_wright_p")]
#> # A tibble: 1 × 2
#>   stock_wright_stat stock_wright_p
#>               <dbl>          <dbl>
#> 1              1.85          0.603
```

### When instruments are genuinely weak

Baum, Schaffer & Stillman (2007, pp. 492–493) illustrate the
weak-instrument problem on the Griliches equation, instrumenting `iq` (a
measure of ability) with `age` and `mrt` alone, under year dummies and a
heteroskedasticity-robust variance:

``` r

weak_form <- lw ~ s + expr + tenure + rns + smsa + factor(year) | iq | age + mrt
fit_weak <- ivreg2(weak_form, data = griliches, vcov = "robust")
sy10 <- with(fit_weak$diagnostics$weak_id_sy,
             critical_value[threshold == "10%"])
summary(fit_weak)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = weak_form, data = griliches, vcov = "robust")
#> 
#> Observations: 758 
#> VCV type:     Robust 
#> 
#> Coefficients:
#>                 Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)    10.550965   2.781762   3.793 0.000149 ***
#> iq             -0.094890   0.041890  -2.265 0.023500 *  
#> s               0.339712   0.118327   2.871 0.004092 ** 
#> expr           -0.006604   0.029255  -0.226 0.821406    
#> tenure          0.084885   0.030668   2.768 0.005643 ** 
#> rns            -0.376939   0.155997  -2.416 0.015678 *  
#> smsa            0.218119   0.103112   2.115 0.034399 *  
#> factor(year)67  0.007775   0.166325   0.047 0.962717    
#> factor(year)68  0.037799   0.152359   0.248 0.804061    
#> factor(year)69  0.334703   0.163799   2.043 0.041016 *  
#> factor(year)70  0.628643   0.246846   2.547 0.010875 *  
#> factor(year)71  0.444610   0.186188   2.388 0.016942 *  
#> factor(year)73  0.439027   0.166866   2.631 0.008513 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      -6.4195 
#> Adj. R-squared: -6.5390 
#> Wald chi2(12):  4.4 (p = 0.0000)
#> Root MSE:       1.1676 
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           2.72 
#>   Kleibergen-Paap rk Wald F:     2.93 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       19.93 
#>      15%  maximal IV size       11.59 
#>      20%  maximal IV size       8.75 
#>      25%  maximal IV size       7.25 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(2) = 5.90 (p = 0.0524)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(2,744) = 46.95 (p < 2.2e-16)
#>   Anderson-Rubin Wald Chi-sq(2) = 95.66 (p < 2.2e-16)
#>   Stock-Wright LM S Chi-sq(2) = 71.33 (p = 0.0000)
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(1) = 1.56 (p = 0.2111)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 56.06 (p = 0.0000)
#>   Tested: iq
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   iq                  2.93    0.0539      0.0073    0.0073      2.93      2.93
#> 
#> Instrumented:          iq 
#> Included instruments:  s, expr, tenure, rns, smsa, factor(year)67, factor(year)68, factor(year)69, factor(year)70, factor(year)71, factor(year)73 
#> Excluded instruments:  age, mrt
```

Because the variance is robust, the relevant weak-identification
statistic is the Kleibergen-Paap rk Wald F (2.93), far below the
Stock-Yogo 10% maximal-size critical value of 19.93. The Stock-Yogo
critical values are tabulated for the Cragg-Donald statistic under
i.i.d. errors, so under a robust variance they apply to the rk statistic
only with caution, or give way to the Staiger-Stock rule of thumb that
the first-stage F should exceed 10 (Baum, Schaffer & Stillman, 2007,
p. 490). With instruments this weak, the Stock-Wright S and
Anderson-Rubin tests provide size-correct inference where 2SLS Wald
tests do not. These robust tests keep correct size but lose power as
instruments weaken (Baum, Schaffer & Stillman, 2007, p. 491), and
because the null is joint, a rejection can reflect either nonzero
coefficients on the endogenous regressors or failure of the
overidentifying restrictions — the test cannot say which.

### Redundant instruments

Redundancy is a separate question from validity: it asks whether an
instrument adds any identifying information once the others are present.
An instrument is redundant when dropping it costs no asymptotic
efficiency, and on a weak equation it is worth knowing whether the few
available instruments are all pulling their weight. Baum, Schaffer &
Stillman (2007, p. 492) test `mrt` for redundancy in exactly this weak
specification with `redundant`:

``` r

fit_redund <- ivreg2(weak_form, data = griliches, vcov = "robust",
                     redundant = "mrt")
glance(fit_redund)[, c("redundancy_stat", "redundancy_p")]
#> # A tibble: 1 × 2
#>   redundancy_stat redundancy_p
#>             <dbl>        <dbl>
#> 1         0.00176        0.967
```

The LM statistic is essentially zero (p = 0.97), the published result:
`mrt` provides no useful information for identifying the coefficient on
`iq`, so identification here rests on `age` alone. Redundancy is about
efficiency rather than validity, and it is a conditional concept — an
instrument redundant given one set of instruments may cease to be so
given another.

### Building the S statistic by hand with `b0`

The `b0` argument evaluates the CUE objective at a fixed parameter
vector rather than optimizing it, and doing so reconstructs the
Stock-Wright S statistic directly. The help file’s construction sets the
endogenous coefficient (`educ`) to zero and the exogenous coefficients
to their OLS values; fixing the exogenous coefficients at OLS is
equivalent to partialling them out, so the objective value at that point
is the S statistic for the joint null:

``` r

ols_fit <- lm(lwage ~ exper + expersq, data = mroz_work)
b0 <- c(educ = 0, coef(ols_fit))
fit_b0 <- ivreg2(
  lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
  data = mroz_work, vcov = "robust", b0 = b0
)
fit_b0$diagnostics$overid$stat  # the objective value, Stata's e(j)
#> [1] 1.661957
```

The identity is easy to confirm against the S statistic that the
ordinary robust fit reports automatically:

``` r

fit_reg <- ivreg2(
  lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
  data = mroz_work, vcov = "robust"
)
stopifnot(isTRUE(all.equal(
  fit_b0$diagnostics$overid$stat,
  fit_reg$diagnostics$stock_wright$stat
)))
c(b0_objective   = fit_b0$diagnostics$overid$stat,
  stock_wright_S = fit_reg$diagnostics$stock_wright$stat)
#>   b0_objective stock_wright_S 
#>       1.661957       1.661957
```

One caution: the p-value stored beside the raw objective in
`fit_b0$diagnostics$overid$p` is calibrated against `L - K`
overidentification degrees of freedom, which is the wrong reference
distribution for this use. Display the raw objective, as Stata’s
`di e(j)` does, and take the S-statistic p-value (0.6454 here, on L1 = 3
degrees of freedom) from the regular fit.

### What the package does and does not offer

The Anderson-Rubin and Stock-Wright statistics test one particular null
— that the endogenous coefficients are zero and the instruments valid.
`ivreg2r` reports these joint tests, but it does not test non-zero nulls
by re-optimizing the objective at a candidate value (the `b0` mechanism
evaluates the objective at a point but is not wired into a search), and
it does not invert the Anderson-Rubin test into a confidence region. A
reader whose Wald confidence interval is untrustworthy under weak
instruments therefore has a size-correct test of the zero null but not,
from this package, a ready weak-instrument-robust confidence set.
Constructing one means inverting a weak-robust test — the Anderson-Rubin
statistic, Kleibergen’s LM statistic, or the conditional
likelihood-ratio test of Moreira (2003), whose properties are analyzed
in Andrews, Moreira & Stock (2006). Stata’s `weakiv` command (Finlay,
Magnusson & Schaffer) builds these confidence sets for linear IV models
with any number of endogenous regressors, and in R the `ivmodel` package
computes Anderson-Rubin and conditional-likelihood-ratio confidence
intervals for the single-endogenous-regressor case; either accepts the
kind of specification used here.

## Two-way clustering

### The nlswork panel

The NLS young-women panel is Stata’s `[XT]`-manual extract
(`webuse nlswork`): person-year observations on women who were 14–24 in
1968 and were interviewed annually through 1988. Panel data of this
shape can be correlated within a person over time and within a year
across persons, and two-way clustering accounts for both.

``` r

nrow(nlswork)
#> [1] 28534
length(unique(nlswork$idcode))  # persons
#> [1] 4711
length(unique(nlswork$year))    # interview years
#> [1] 15
```

### One-way versus two-way clustering

Compare one-way clustering by person with two-way clustering by person
and interview year, the help file’s own paired example. Fitting on the
full data frame restricts to the complete cases the formula needs, so no
manual subsetting is required:

``` r

fit_1way <- ivreg2(
  ln_wage ~ grade + age + ttl_exp + tenure,
  data = nlswork, clusters = ~ idcode
)

fit_2way <- ivreg2(
  ln_wage ~ grade + age + ttl_exp + tenure,
  data = nlswork, clusters = ~ idcode + year
)
```

Of the 28,534 person-year rows, 28,099 have complete data and enter both
fits. Compare the standard errors:

``` r

se_comparison <- data.frame(
  term = names(coef(fit_1way)),
  se_1way = tidy(fit_1way)$std.error,
  se_2way = tidy(fit_2way)$std.error
)
se_comparison$ratio <- se_comparison$se_2way / se_comparison$se_1way
se_comparison
#>          term      se_1way     se_2way    ratio
#> 1 (Intercept) 0.0337636175 0.043107155 1.276734
#> 2       grade 0.0021612566 0.002654424 1.228185
#> 3         age 0.0009445326 0.001644370 1.740936
#> 4     ttl_exp 0.0018369910 0.002677937 1.457785
#> 5      tenure 0.0016316256 0.003051750 1.870374
```

Two-way clustering uses the Cameron, Gelbach & Miller (2011)
decomposition `V_twoway = V_cluster1 + V_cluster2 - V_intersection`.
This estimator needs both cluster dimensions to grow for consistency, so
the effective number of clusters governing the asymptotics is the
smaller of the two, `min(M1, M2)`. Here there are 4,697 person clusters
and only 15 interview-year clusters, so `min(M1, M2)` is 15: the year
dimension, not the much larger person dimension, limits how far the
two-way asymptotics can be trusted.

Two-way clustering helps only when the two variables define distinct,
crossed categories, as here, where every person appears in several years
and every year holds many persons. That panel-by-time pairing is the
typical two-way application. When one category is nested within the
other — counties within states, say — the two-way variance reduces to
the single-cluster variance on the outer category, and the second
dimension adds nothing.

### Two-way clustering with instrumental variables

The help file’s two-way example is estimated by ordinary least squares,
but two-way clustering composes with every estimator, so we pair it with
2SLS on a panel demand equation. The Baltagi cigarette-demand data
(`data(cigar)`) are an annual panel of 46 U.S. states over 1963–1992,
assembled by Baltagi and Levin (1992) and revisited by Baltagi, Griffin
and Xiong (2000).

The specification regresses log per-capita cigarette sales on the log
real price of cigarettes, with log real per-capita disposable income as
an exogenous control. Price is plausibly endogenous, because state
demand shocks feed back into prices, so it is instrumented by the
minimum real price in adjoining states — a spatial cost-shifter in the
Baltagi-Levin tradition — together with the one-year lag of the state’s
own real price. Clustering on both state and year allows arbitrary
within-state serial correlation and arbitrary within-year correlation
across states.

``` r

data(cigar)
cigar <- transform(cigar,
  lsales  = log(sales),
  lrprice = log(price / cpi),
  lrndi   = log(ndi / cpi),
  lrpimin = log(pimin / cpi)
)

fit_iv_2way <- ivreg2(
  lsales ~ lrndi | lrprice | lrpimin + l(lrprice, 1),
  data = cigar, tvar = "year", ivar = "state",
  clusters = ~ state + year
)
summary(fit_iv_2way)
#> 
#> 2SLS Estimation
#> 
#> Call:
#> ivreg2(formula = lsales ~ lrndi | lrprice | lrpimin + l(lrprice, 
#>     1), data = cigar, clusters = ~state + year, tvar = "year", 
#>     ivar = "state")
#> 
#> Observations: 1,334 
#> VCV type:     Cluster-robust 
#> Clusters:     46 (state), 29 (year)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  3.42880    0.34460   9.950  < 2e-16 ***
#> lrprice     -0.84886    0.11442  -7.419 1.18e-13 ***
#> lrndi        0.27957    0.07535   3.710 0.000207 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> R-squared:      0.3271 
#> Adj. R-squared: 0.3261 
#> Wald chi2(2):  27.7 (p = 0.0000)
#> Root MSE:       0.1840 
#> 
#> Weak identification test:
#>   Cragg-Donald Wald F:           4766.55 
#>   Kleibergen-Paap rk Wald F:     275.63 
#>   (Stock-Yogo critical values are for iid errors)
#>   Stock-Yogo critical values (IV size):
#>      10%  maximal IV size       19.93 
#>      15%  maximal IV size       11.59 
#>      20%  maximal IV size       8.75 
#>      25%  maximal IV size       7.25 
#> 
#> Underidentification test (Kleibergen-Paap rk LM statistic):
#>   Chi-sq(2) = 12.14 (p = 0.0023)
#> 
#> Weak-instrument-robust inference:
#>   H0: B1=0 and orthogonality conditions are valid
#>   Anderson-Rubin Wald F(2,28) = 27.10 (p = 0.0000)
#>   Anderson-Rubin Wald Chi-sq(2) = 56.26 (p = 0.0000)
#>   Stock-Wright LM S Chi-sq(2) = 12.43 (p = 0.0020)
#> 
#> Overidentification test (Hansen J):
#>   Chi-sq(1) = 0.94 (p = 0.3325)
#> 
#> Endogeneity test:
#>   Chi-sq(1) = 0.09 (p = 0.7636)
#>   Tested: lrprice
#> 
#> First-stage diagnostics:
#>   Endogenous        F-stat   p-value  Partial R2  Shea PR2      SW F      AP F
#>   lrprice           275.63    0.0000      0.8776    0.8776    275.63    275.63
#> 
#> Instrumented:          lrprice 
#> Included instruments:  lrndi 
#> Excluded instruments:  lrpimin, l(lrprice, 1)
```

The `l(lrprice, 1)` term uses the panel lag operator, which needs the
panel structure declared through `tvar` and `ivar`, and which is
documented in
[`vignette("time-series-gmm")`](https://restatr.com/ivreg2r/articles/time-series-gmm.md).
The summary reports both cluster dimensions — 46 states and 29 years,
the first year lost to the lag — and the effective cluster count
governing the asymptotics is again the smaller dimension, the 29 years.
The automatic diagnostics carry through the two-way variance unchanged:
the Kleibergen-Paap underidentification and weak-identification
statistics and the Hansen J of the single overidentifying restriction
are all computed from the two-way cluster-robust covariance. This is a
methods illustration rather than a substantive demand study, and the
price elasticity should be read in that light.

## Reduced-form regression

The reduced form regresses the outcome directly on the full instrument
set — exogenous regressors plus excluded instruments — bypassing the
first stage.

### Single-equation reduced form

Use `reduced_form = "rf"` to store the outcome-equation reduced form
alongside the structural estimates:

``` r

# Same robust fit as fit_reg above, re-estimated with the reduced form stored
fit_rf <- update(fit_reg, reduced_form = "rf")
fit_rf$reduced_form$coefficients
#>   (Intercept)         exper       expersq           age       kidslt6 
#>  0.9940565569  0.0452258622 -0.0009716342 -0.0027608640  0.0175153014 
#>       kidsge6 
#> -0.0385237413
```

The stored reduced form is not just bookkeeping: its F-statistic on the
excluded instruments *is* the Anderson-Rubin F. The AR test is precisely
this reduced-form F-test, which is why it stays valid under weak
instruments, and the identity holds exactly:

``` r

stopifnot(isTRUE(all.equal(fit_rf$reduced_form$f_stat,
                           fit_rf$diagnostics$anderson_rubin$f_stat)))
c(rf_F = fit_rf$reduced_form$f_stat,
  ar_F = fit_rf$diagnostics$anderson_rubin$f_stat)
#>      rf_F      ar_F 
#> 0.5795624 0.5795624
```

### System reduced form

Use `reduced_form = "system"` to estimate the full system — one equation
for each endogenous variable plus the outcome — with a cross-equation
covariance matrix. The cross-equation covariance is what joint
hypotheses spanning the first-stage and outcome reduced forms require,
and storing the system reproduces Stata’s undocumented `savesfirst`
behavior:

``` r

fit_system <- ivreg2(
  lwage ~ exper + expersq | educ | age + kidslt6 + kidsge6,
  data = mroz_work, reduced_form = "system"
)

# Equation names
fit_system$reduced_form$depvar
#> [1] "lwage" "educ"

# Cross-equation VCV dimensions
dim(fit_system$reduced_form$vcov)
#> [1] 12 12
```

## Degrees-of-freedom adjustments

When fixed effects or other controls have been partialled out before
estimation, `dofminus` and `sdofminus` adjust the degrees of freedom
without re-estimating:

- `dofminus` is subtracted from `N` in the large-sample variance
  formulas, so `sigma^2 = RSS / (N - dofminus)`. It is the large-sample
  adjustment (for example, the number of fixed effects), and it is how
  `xtivreg2, fe` reports correct statistics after within-demeaning
  consumes one degree of freedom per panel unit.
- `sdofminus` is subtracted from the residual degrees of freedom
  alongside `K` in the small-sample corrections. It is the small-sample
  adjustment (for example, the number of partialled-out regressors), and
  it is what Stata’s `partial()` sets internally.

Neither option appears in the help file’s examples, so both
demonstrations below are self-verifying identities on the Griliches
data. First `sdofminus`, via the Frisch-Waugh-Lovell theorem.
Residualize every variable on the year dummies (7 columns including the
intercept), refit without an intercept, and tell
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md) that 7
degrees of freedom were consumed. The shared coefficients, standard
errors, and sigma then match the full regression to floating-point
precision:

``` r

fit_full <- ivreg2(
  lw ~ s + expr + tenure + rns + smsa + factor(year) | iq | med + kww + age,
  data = griliches, small = TRUE
)

# Partial the year dummies out of every variable by hand
yd <- model.matrix(~ factor(year), data = griliches)
resid_on_years <- function(v) lm.fit(yd, v)$residuals
gril_fwl <- data.frame(lapply(
  griliches[c("lw", "s", "expr", "tenure", "rns", "smsa",
              "iq", "med", "kww", "age")],
  resid_on_years
))

fit_fwl <- ivreg2(
  lw ~ 0 + s + expr + tenure + rns + smsa | iq | med + kww + age,
  data = gril_fwl, small = TRUE, sdofminus = 7
)

shared <- names(coef(fit_fwl))
max(abs(coef(fit_full)[shared] - coef(fit_fwl)[shared]))
#> [1] 3.330669e-16
max(abs(sqrt(diag(vcov(fit_full)))[shared] - sqrt(diag(vcov(fit_fwl)))[shared]))
#> [1] 3.469447e-18
c(sigma_full = fit_full$sigma, sigma_fwl = fit_fwl$sigma)
#> sigma_full  sigma_fwl 
#>  0.3266951  0.3266951
```

Without `sdofminus = 7` the coefficients would still match but the
standard errors and sigma would not, because the refit would not know
that 7 parameters were already spent on the year dummies. The built-in
`partial = "factor(year)"` argument performs exactly this
residualize-and-adjust sequence in one step; see
[`vignette("time-series-gmm")`](https://restatr.com/ivreg2r/articles/time-series-gmm.md)
for the partialling workflow.

The `dofminus` twin is the fixed-effects recipe that `xtivreg2, fe`
automates. Demean each variable by its group, adding back the grand mean
so the intercept survives, and declare that the G - 1 absorbed group
effects consumed degrees of freedom. Treating the seven Griliches survey
years as the groups, the within fit reproduces the dummy-variable fit
exactly — shared coefficients, standard errors, sigma, and residual
degrees of freedom:

``` r

demean_gm <- function(v, g) v - ave(v, g) + mean(v)
gril_within <- data.frame(lapply(
  griliches[c("lw", "s", "expr", "tenure")],
  demean_gm, g = griliches$year
))
G <- length(unique(griliches$year))

fit_lsdv <- ivreg2(lw ~ s + expr + tenure + factor(year),
                   data = griliches, small = TRUE)
fit_fe <- ivreg2(lw ~ s + expr + tenure, data = gril_within,
                 small = TRUE, dofminus = G - 1)

slopes <- c("s", "expr", "tenure")
max(abs(coef(fit_lsdv)[slopes] - coef(fit_fe)[slopes]))
#> [1] 2.081668e-17
max(abs(sqrt(diag(vcov(fit_lsdv)))[slopes] - sqrt(diag(vcov(fit_fe)))[slopes]))
#> [1] 2.602085e-18
c(sigma_lsdv = fit_lsdv$sigma, sigma_fe = fit_fe$sigma,
  df_lsdv = fit_lsdv$df.residual, df_fe = fit_fe$df.residual)
#> sigma_lsdv   sigma_fe    df_lsdv      df_fe 
#>   0.337452   0.337452 748.000000 748.000000
```

For a migrating `xtivreg2, fe` user the practical rule is: pass the
number of absorbed effects beyond the constant, and use the identity
above as the check that you passed the right one.

## Estimator equivalences

Two equivalences involving LIML are worth stating here; the
authoritative table of variance-estimator equivalences (Driscoll-Kraay,
Kiefer, cluster-robust) lives in
[`vignette("time-series-gmm")`](https://restatr.com/ivreg2r/articles/time-series-gmm.md).

- **LIML with COVIV = CUE** under conditional-homoskedasticity
  weighting: demonstrated on the Klein data above, with agreement up to
  optimizer tolerance.
- **LIML = 2SLS when exactly identified**, as noted in the
  estimator-choice discussion above.

## Publication tables with modelsummary

The [`modelsummary`](https://modelsummary.com/) package produces
publication-quality regression tables from any model with
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) methods,
and it is where the estimator-comparison workflow belongs. `ivreg2r`
objects work with it directly.

``` r

library(modelsummary)
```

### A single model

``` r

modelsummary(fit_2sls, output = "markdown")
```

[TABLE]

### Comparing estimators side by side

The natural side-by-side comparison is the Klein suite, where the
estimator choice visibly matters — reusing the `klein_models` list built
above:

``` r

modelsummary(klein_models, output = "markdown",
             gof_map = c("nobs", "r.squared", "adj.r.squared", "sigma"))
```

|               | 2SLS    | LIML    | Fuller(1) | Nagar   |
|---------------|---------|---------|-----------|---------|
| (Intercept)   | 16.555  | 17.148  | 17.008    | 16.711  |
|               | (1.321) | (1.840) | (1.702)   | (1.437) |
| profits       | 0.017   | -0.223  | -0.169    | -0.050  |
|               | (0.118) | (0.202) | (0.180)   | (0.137) |
| wagetot       | 0.810   | 0.823   | 0.820     | 0.814   |
|               | (0.040) | (0.055) | (0.051)   | (0.044) |
| l(profits, 1) | 0.216   | 0.396   | 0.355     | 0.266   |
|               | (0.107) | (0.174) | (0.156)   | (0.122) |
| Num.Obs.      | 21      | 21      | 21        | 21      |
| R2            | 0.977   | 0.957   | 0.963     | 0.973   |
| R2 Adj.       | 0.973   | 0.949   | 0.956     | 0.968   |
| Sigma         | 1.022   | 1.395   | 1.296     | 1.105   |

### Exponentiated coefficients

For log-linear models, `exponentiate = TRUE` reports multiplicative
effects (`exp(beta)`). `modelsummary` transforms the standard errors to
the exponentiated scale by the delta method, so the displayed standard
error is `exp(beta) * se(beta)` and the table is internally consistent
on the multiplicative scale. This transformation is `modelsummary`’s
own: the package’s `tidy(fit, exponentiate = TRUE)` follows the broom
convention of exponentiating estimates and confidence bounds while
leaving `std.error` on the log scale.

``` r

modelsummary(fit_2sls, output = "markdown", exponentiate = TRUE)
```

[TABLE]

### Customizing output

Rename coefficients and select goodness-of-fit statistics:

``` r

modelsummary(
  list("2SLS" = fit_2sls, "LIML" = fit_liml),
  output = "markdown",
  coef_rename = c(educ = "Education", exper = "Experience",
                  expersq = "Experience²"),
  gof_map = c("nobs", "r.squared", "sigma", "weak_id_stat")
)
```

|              | 2SLS    | LIML    |
|--------------|---------|---------|
| (Intercept)  | -0.385  | -0.377  |
|              | (1.012) | (1.039) |
| Education    | 0.096   | 0.096   |
|              | (0.081) | (0.084) |
| Experience   | 0.042   | 0.042   |
|              | (0.014) | (0.014) |
| Experience²  | -0.001  | -0.001  |
|              | (0.000) | (0.000) |
| Num.Obs.     | 428     | 428     |
| R2           | 0.156   | 0.155   |
| Sigma        | 0.664   | 0.664   |
| weak_id_stat | 4       | 4       |

`modelsummary` supports LaTeX, HTML, Word, and gt output; see
[`?modelsummary::modelsummary`](https://modelsummary.com/man/modelsummary.html)
for the full list.

## Stata migration

For a Stata-to-R translation table covering the most common options, see
the migration section in
[`vignette("introduction")`](https://restatr.com/ivreg2r/articles/introduction.md);
for the full argument list, see
[`?ivreg2`](https://restatr.com/ivreg2r/reference/ivreg2.md). Factor
variables are handled through R’s formula interface: `factor(year)` in
the Griliches specifications above plays the role of Stata’s `xi i.year`
dummies, expanding to one dummy per non-reference level under treatment
contrasts.

## References

- Anderson, T.W. and Rubin, H. (1949). “Estimation of the Parameters of
  a Single Equation in a Complete System of Stochastic Equations.”
  *Annals of Mathematical Statistics*, 20(1), 46–63.
- Anderson, T.W. and Rubin, H. (1950). “The Asymptotic Properties of
  Estimates of the Parameters of a Single Equation in a Complete System
  of Stochastic Equations.” *Annals of Mathematical Statistics*, 21(4),
  570–582.
- Andrews, D.W.K., Moreira, M.J., and Stock, J.H. (2006). “Optimal
  Two-Sided Invariant Similar Tests for Instrumental Variables
  Regression.” *Econometrica*, 74(3), 715–752.
- Baltagi, B.H. and Levin, D. (1992). “Cigarette Taxation: Raising
  Revenues and Reducing Consumption.” *Structural Change and Economic
  Dynamics*, 3(2), 321–335.
- Baltagi, B.H., Griffin, J.M., and Xiong, W. (2000). “To Pool or Not to
  Pool: Homogeneous Versus Heterogeneous Estimators Applied to Cigarette
  Demand.” *Review of Economics and Statistics*, 82(1), 117–126.
- Baum, C.F., Schaffer, M.E., and Stillman, S. (2007). “Enhanced
  Routines for Instrumental Variables/Generalized Method of Moments
  Estimation and Testing.” *Stata Journal*, 7(4), 465–506.
- Cameron, A.C., Gelbach, J.B., and Miller, D.L. (2011). “Robust
  Inference with Multiway Clustering.” *Journal of Business & Economic
  Statistics*, 29(2), 238–249.
- Cragg, J.G. (1983). “More Efficient Estimation in the Presence of
  Heteroscedasticity of Unknown Form.” *Econometrica*, 51(3), 751–763.
- Davidson, R. and MacKinnon, J.G. (1993). *Estimation and Inference in
  Econometrics*. Oxford University Press.
- Fuller, W.A. (1977). “Some Properties of a Modification of the Limited
  Information Estimator.” *Econometrica*, 45(4), 939–953.
- Griliches, Z. (1976). “Wages of Very Young Men.” *Journal of Political
  Economy*, 84(4, Part 2), S69–S85.
- Hayashi, F. (2000). *Econometrics*. Princeton University Press.
- Klein, L.R. (1950). *Economic Fluctuations in the United States,
  1921–1941*. Cowles Commission Monograph No. 11. Wiley.
- Mariano, R.S. and Sawa, T. (1972). “The Exact Finite-Sample
  Distribution of the Limited-Information Maximum Likelihood Estimator
  in the Case of Two Included Endogenous Variables.” *Journal of the
  American Statistical Association*, 67(337), 159–163.
- Moreira, M.J. (2003). “A Conditional Likelihood Ratio Test for
  Structural Models.” *Econometrica*, 71(4), 1027–1048.
- Mroz, T.A. (1987). “The Sensitivity of an Empirical Model of Married
  Women’s Hours of Work to Economic and Statistical Assumptions.”
  *Econometrica*, 55(4), 765–799.
- Nagar, A.L. (1959). “The Bias and Moment Matrix of the General k-Class
  Estimators of the Parameters in Simultaneous Equations.”
  *Econometrica*, 27(4), 575–595.
- Sargan, J.D. (1958). “The Estimation of Economic Relationships Using
  Instrumental Variables.” *Econometrica*, 26(3), 393–415.
- Stock, J.H. and Wright, J.H. (2000). “GMM with Weak Identification.”
  *Econometrica*, 68(5), 1055–1096.
- Stock, J.H. and Yogo, M. (2005). “Testing for Weak Instruments in
  Linear IV Regression.” In D.W.K. Andrews and J.H. Stock (Eds.),
  *Identification and Inference for Econometric Models: Essays in Honor
  of Thomas Rothenberg*. Cambridge University Press.
