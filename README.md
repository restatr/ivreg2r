# ivreg2r

Stata's [`ivreg2`](https://ideas.repec.org/c/boc/bocode/s425401.html) for R.
Comprehensive IV/GMM estimation with automatic diagnostics --- same estimators,
same tests, matching numbers. Numerical parity verified against Stata across
900+ fixture files.

## Installation

```r
# install.packages("devtools")
devtools::install_github("restatr/ivreg2r")
```

## A worked example

Card (1995) estimates returns to schooling, instrumenting education with
proximity to two- and four-year colleges. With two excluded instruments and one
endogenous variable, the model is overidentified --- so we get diagnostics for
both instrument strength and validity.

```r
library(ivreg2r)
data(card)

fit_2sls <- ivreg2(
  lwage ~ exper + expersq + black + south + smsa + married |
    educ | nearc2 + nearc4,
  data = card, vcov = "robust", small = TRUE
)
summary(fit_2sls)
```

The diagnostics tell a clear story:

- **Kleibergen-Paap rk F ≈ 9.7** --- below the Stock-Yogo 15% maximal-size
  critical value of 11.59 (tabulated for the Cragg-Donald statistic under iid
  and used here as an informal benchmark), suggesting the 2SLS Wald test may
  suffer from size distortion.
- **Hansen J p ≈ 0.10** --- does not reject the overidentifying restrictions
  at 5%. The J test has limited power, particularly with weak instruments, and
  only tests the *overidentifying* restrictions --- it cannot detect problems
  in exactly identified directions.
- **Endogeneity C-stat p ≈ 0.05** --- borderline evidence that education is
  endogenous (the null is exogeneity, conditional on the instruments being
  valid in the full model).

Since the KP F falls below the 2SLS Stock-Yogo thresholds, LIML is the
natural next step. LIML is more robust to weak instruments than 2SLS when
overidentified --- its approximate bias does not grow with the number of
instruments, and the Stock-Yogo size-distortion critical values for LIML are
lower than for 2SLS:

```r
fit_liml <- ivreg2(
  lwage ~ exper + expersq + black + south + smsa + married |
    educ | nearc2 + nearc4,
  data = card, method = "liml", vcov = "robust", small = TRUE
)
```

Programmatic access to all diagnostics makes model comparison easy --- this is
where R shines over Stata:

```r
comparison <- rbind(glance(fit_2sls), glance(fit_liml))
comparison[, c("method", "weak_id_robust_stat", "overid_stat", "overid_p",
               "endogeneity_p")]
```

## Formula syntax

Three-part formula: `y ~ exogenous | endogenous | excluded_instruments`. Each
variable appears once; exogenous regressors are automatically included as
instruments.

| Part | Contents | Example |
|------|----------|---------|
| LHS | Dependent variable | `lwage` |
| 1st RHS | Exogenous regressors | `exper + expersq + black + south` |
| 2nd RHS | Endogenous regressors | `educ` |
| 3rd RHS | Excluded instruments | `nearc2 + nearc4` |

A one-part formula (`y ~ x1 + x2`) estimates OLS with the same diagnostics.

**Stata equivalent:**

| Stata | R |
|-------|---|
| `ivreg2 lwage exper ... (educ = nearc2 nearc4)` | `ivreg2(lwage ~ exper + ... \| educ \| nearc2 + nearc4, data = d)` |

## Stata → R quick reference

| Stata | R (`ivreg2r`) |
|-------|---------------|
| `ivreg2 ..., robust` | `ivreg2(..., vcov = "robust")` |
| `ivreg2 ..., robust small` | `ivreg2(..., vcov = "robust", small = TRUE)` |
| `ivreg2 ..., cluster(g)` | `ivreg2(..., clusters = ~ g, small = TRUE)` |
| `ivreg2 ..., cluster(g1 g2)` | `ivreg2(..., clusters = ~ g1 + g2)` |
| `ivreg2 ..., liml` | `ivreg2(..., method = "liml")` |
| `ivreg2 ..., fuller(1)` | `ivreg2(..., fuller = 1)` |
| `ivreg2 ..., gmm2s` | `ivreg2(..., method = "gmm2s")` |
| `ivreg2 ..., cue` | `ivreg2(..., method = "cue")` |
| `ivreg2 ..., kernel(bar) bw(3)` | `ivreg2(..., kernel = "bartlett", bw = 3, tvar = "year")` |
| `ivreg2 ..., orthog(z1)` | `ivreg2(..., orthog = c("z1"))` |
| `ivreg2 ..., endog(x1)` | `ivreg2(..., endog = c("x1"))` |
| `ivreg2 ..., redundant(z1)` | `ivreg2(..., redundant = c("z1"))` |
| `ivreg2 ..., partial(x1)` | `ivreg2(..., partial = c("x1"))` |
| `ivreg2 ... [aw=w]` | `ivreg2(..., weights = w)` |
| `ivreg2 ... [fw=w]` | `ivreg2(..., weights = w, weight_type = "fweight")` |
| `ivreg2 ... [pw=w]` | `ivreg2(..., weights = w, weight_type = "pweight")` |
| `e(b)` / `e(V)` | `coef(fit)` / `vcov(fit)` |
| `ereturn list` | `glance(fit)` |

## Stata parity

All outputs match Stata within tight numerical tolerances (coefficients/SEs to
1e-6 relative, test statistics to 1e-4 relative). Verified against 900+ Stata
fixture files and 8,600+ automated tests. Stata statistical conventions are
followed by default (e.g., sigma normalization for weighted models).

## Learn more

- `vignette("introduction")` --- Getting started: OLS, 2SLS, diagnostics,
  robust/cluster SEs, weights, broom integration, Stata migration reference.
- `vignette("advanced-iv")` --- LIML, Fuller, k-class, COVIV, orthogonality
  tests, Stock-Wright, two-way clustering, reduced-form, model comparison.
- `vignette("time-series-gmm")` --- HAC/AC standard errors, Kiefer,
  Driscoll-Kraay, GMM, CUE, partialling, user matrices, panel VCE.

## License

MIT
