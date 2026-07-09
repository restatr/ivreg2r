# Glossary of diagnostic-output acronyms

[`summary()`](https://rdrr.io/r/base/summary.html) and
[`print()`](https://rdrr.io/r/base/print.html) report a battery of
identification, weak-instrument, and overidentification diagnostics,
many named by acronym. This topic glosses every acronym that reaches the
console. The underlying statistics are stored on the fitted object under
`fit$diagnostics`; the headline ones are also in
[`glance()`](https://generics.r-lib.org/reference/glance.html).

## Identification and weak instruments

- Underidentification test:

  Tests whether the excluded instruments have rank high enough to
  identify the endogenous regressors. Reported as the Anderson (1951)
  canonical-correlations LM statistic under `vcov = "iid"`, and the
  Kleibergen-Paap rk LM statistic under robust, cluster, or HAC errors.

- Weak identification test:

  Gauges instrument strength against the Stock-Yogo (2005) critical
  values. Reported as the Cragg-Donald (1993) Wald F under
  `vcov = "iid"`, and the Kleibergen-Paap rk Wald F under robust,
  cluster, or HAC errors.

- Cragg-Donald (CD):

  The iid weak-identification Wald F (and the underlying eigenvalue
  statistic); see "Weak identification test".

- Kleibergen-Paap (KP) rk:

  Rank-based LM and Wald statistics that replace the Anderson and
  Cragg-Donald statistics when errors are not iid.

- Anderson-Rubin (AR):

  A weak-instrument-robust test of the joint significance of the
  endogenous regressors: valid regardless of instrument strength.

- Stock-Wright (S):

  A weak-instrument-robust statistic for the same hypothesis, the score
  (LM) counterpart of the Anderson-Rubin Wald test. It is distinct from
  the Hansen J overidentification statistic: S constrains the endogenous
  coefficients to zero and has degrees of freedom equal to the number of
  excluded instruments.

## First-stage diagnostics

- Shea partial R-squared:

  A partial R-squared for the first stage that accounts for correlation
  among instruments; the relevant measure with more than one endogenous
  regressor.

- SW F (Sanderson-Windmeijer):

  A conditional first-stage F statistic for the identification of an
  individual endogenous regressor, given the others.

- AP (Angrist-Pischke):

  A conditional first-stage F (or chi-squared, depending on the VCE) for
  an individual endogenous regressor.

## Overidentification and endogeneity

- Sargan / Hansen J:

  The overidentification test: Sargan (1958) under `vcov = "iid"`, and
  the Hansen J statistic under robust, cluster, or HAC errors. Tests the
  joint validity of the overidentifying restrictions.

- C-statistic (C-stat):

  A difference-of-J statistic used for the endogeneity test (`endog`),
  the orthogonality test (`orthog`), and the redundancy test
  (`redundant`).

## VCE and estimator labels

- Kiefer:

  The Kiefer (1980) VCE: autocorrelation-consistent with a Truncated
  kernel at bandwidth equal to the full time span.

- Driscoll-Kraay (DK):

  The Driscoll-Kraay (1998) panel VCE, robust to cross-sectional
  dependence.

- dofminus / sdofminus:

  Large- and small-sample degrees-of-freedom adjustments; see
  [ivreg2r-conventions](https://restatr.com/ivreg2r/reference/ivreg2r-conventions.md).

- psd0 / psda:

  Positive-semidefinite corrections for the moment covariance matrix,
  zeroing (`psd0`) or taking the absolute value of (`psda`) any negative
  eigenvalues.

- COVIV:

  The "covariance at the IV estimates" variance matrix for LIML and
  k-class estimators, robust to misspecification of the LIML model.

- HOLS:

  Cragg's (1983) heteroskedastic OLS estimator, obtained from a
  no-endogenous-regressor model under `method = "gmm2s"` with a robust
  VCE.

- k-class:

  The family of IV estimators indexed by a scalar k (2SLS at k = 1, OLS
  at k = 0, LIML at k equal to the LIML eigenvalue).

## See also

[`summary.ivreg2()`](https://restatr.com/ivreg2r/reference/summary.ivreg2.md)
for the console output;
[ivreg2r-conventions](https://restatr.com/ivreg2r/reference/ivreg2r-conventions.md)
for the statistical conventions;
[`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md) for the
arguments that request each test.

Other ivreg2r reference:
[`ivreg2r-conventions`](https://restatr.com/ivreg2r/reference/ivreg2r-conventions.md)
