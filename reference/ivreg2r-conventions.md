# Statistical conventions in ivreg2r

`ivreg2r` is a translation of Stata's `ivreg2`, and it follows Stata's
statistical conventions rather than those of R's
[`lm()`](https://rdrr.io/r/stats/lm.html)/[`vcov()`](https://rdrr.io/r/stats/vcov.html)
ecosystem. This topic collects those conventions in one place. The API
itself is idiomatic R (formula interface, tidyverse-friendly output);
only the *statistical* choices below lean on Stata. See
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md) for the
estimation function and
[ivreg2r-glossary](https://restatr.com/ivreg2r/reference/ivreg2r-glossary.md)
for the diagnostic acronyms printed by
[`summary()`](https://rdrr.io/r/base/summary.html).

## Default VCE and small-sample corrections

`vcov` defaults to `"iid"` (classical) and `small` defaults to `FALSE`,
matching Stata's `ivreg2` defaults. With `small = FALSE`, inference uses
z and chi-squared statistics and the large-sample variance
normalization. Set `small = TRUE` for t and F statistics and the
finite-sample `N - K` correction, matching Stata's `small` option.
`small` is the single, uniform control for finite-sample corrections
across every VCE type (`"iid"`, `"robust"`, `"HAC"`, `"AC"`); it is not
a VCE selector. There are no `"HC0"`/`"HC1"` values — use
`vcov = "robust"` and toggle `small`.

Three statistics sit outside that uniform control, all matching Stata:
the first-stage F statistics
([`first_stage()`](https://restatr.com/ivreg2r/reference/first_stage.md))
always use small-sample degrees of freedom regardless of `small`; the
overidentification (Sargan/Hansen J) and endogeneity (C-statistic) tests
are never small-sample corrected, even when `small = TRUE`; and the
stored reduced-form variance-covariance matrix (`reduced_form`) always
applies OLS-style small-sample scaling regardless of the main model's
`small`.

## Sigma (RMSE) normalization

Under `small = FALSE`, \\\sigma^2 = RSS / (N - dofminus)\\. Under
`small = TRUE`, \\\sigma^2 = RSS / (N - K - dofminus - sdofminus)\\.
This differs from [`lm()`](https://rdrr.io/r/stats/lm.html), which
always divides the residual sum of squares by its residual degrees of
freedom. Because ivreg2r reports the large-sample sigma by default, its
RMSE and any quantity built from sigma differ from
[`lm()`](https://rdrr.io/r/stats/lm.html) unless `small = TRUE`.

## Standard errors relative to lm()

Coefficients match [`lm()`](https://rdrr.io/r/stats/lm.html) for OLS.
Standard errors and the full variance-covariance matrix match
[`lm()`](https://rdrr.io/r/stats/lm.html) only when `small = TRUE`; the
default `small = FALSE` reports Stata-convention standard errors without
the `N - K` finite-sample correction. Under analytic (`"aweight"`) or
probability (`"pweight"`) weights, sigma additionally differs from
`lm(..., weights = w)` by a factor of `sqrt(N / sum(w))`, because
[`lm()`](https://rdrr.io/r/stats/lm.html) uses raw (unnormalized)
weights while ivreg2r normalizes them (see below).

## Weight normalization

Analytic and probability weights are normalized internally to sum to
`N`, following Stata. This makes sigma scale-invariant for those types:
multiplying every weight by a constant leaves all coefficients, standard
errors, and test statistics unchanged. Frequency (`"fweight"`) and
importance (`"iweight"`) weights are *not* normalized; they redefine `N`
as `sum(weights)`, so their scale is meaningful by construction.
Probability weights force a robust VCE. See the `weight_type` argument
of [`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md).

## Degrees-of-freedom adjustments

`dofminus` is a large-sample adjustment subtracted from `N` in
large-sample variance formulas (for example when fixed effects have been
partialled out upstream). `sdofminus` is a small-sample adjustment
subtracted from the residual degrees of freedom alongside `K` (for
example the count of regressors removed by `partial`). Both mirror
Stata's `dofminus()` and `sdofminus()` options. Stata's own source
labels `dofminus` a large-sample adjustment (e.g. number of fixed
effects) and `sdofminus` a small-sample one (e.g. number of
partialled-out regressors).

## Confidence intervals and test statistics

[`confint()`](https://rdrr.io/r/stats/confint.html) and the model Wald
test use the standard normal and chi-squared distributions when
`small = FALSE`, and the t and F distributions when `small = TRUE`. The
reported model test is therefore a chi-squared Wald statistic by default
and an F statistic under `small = TRUE`.

## See also

[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md) for the
estimation function and the full argument list;
[ivreg2r-glossary](https://restatr.com/ivreg2r/reference/ivreg2r-glossary.md)
for the diagnostic-output acronyms.

Other ivreg2r reference:
[`ivreg2r-glossary`](https://restatr.com/ivreg2r/reference/ivreg2r-glossary.md)
