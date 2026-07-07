# Model F-test (all slopes = 0)

Computes the Wald F-statistic for the joint null that all slope
coefficients equal zero. The intercept (if present) is excluded from the
test.

## Usage

``` r
.compute_model_f(
  coefficients,
  vcov,
  N,
  K,
  has_intercept,
  vcov_type,
  small,
  M = NULL,
  dofminus = 0L,
  sdofminus = 0L
)
```

## Arguments

- coefficients:

  Named numeric vector of coefficient estimates.

- vcov:

  Variance-covariance matrix of coefficients.

- N:

  Integer: number of observations.

- K:

  Integer: number of regressors (including intercept if present).

- has_intercept:

  Logical: does the model include an intercept?

- vcov_type:

  Character: `"iid"`, `"robust"`, `"HAC"`, `"AC"`, or `"CL"`.

- small:

  Logical: whether small-sample corrections were applied.

- M:

  Integer or NULL: number of clusters (only for `vcov_type = "CL"`).

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

## Value

Named list with `model_f`, `model_f_p`, `model_f_df1`, `model_f_df2`
(all `NA` when df1 = 0 or VCV is not positive definite).

## Note

The model F-statistic is always reported as an F-distributed statistic.
For non-IID VCE types, the underlying Wald chi-squared is converted to F
using VCE-branch-specific degrees of freedom adjustments. This matches
Stata's `ivreg2` behavior.
