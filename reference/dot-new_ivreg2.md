# Construct an ivreg2 object

Construct an ivreg2 object

## Usage

``` r
.new_ivreg2(
  coefficients,
  residuals,
  fitted.values,
  vcov,
  sigma,
  df.residual,
  rank,
  r.squared,
  adj.r.squared,
  rss,
  r2u,
  r2c,
  mss,
  model_f = NULL,
  model_f_p = NULL,
  model_f_df1 = NULL,
  model_f_df2 = NULL,
  diagnostics = NULL,
  first_stage = NULL,
  first_stage_models = NULL,
  reduced_form = NULL,
  call,
  formula,
  terms,
  nobs,
  vcov_type,
  small,
  dofminus = 0L,
  sdofminus = 0L,
  cluster_var = NULL,
  n_clusters = NULL,
  n_clusters1 = NULL,
  n_clusters2 = NULL,
  na.action = NULL,
  weights = NULL,
  weight_type = "aweight",
  n_physical = NULL,
  endogenous = character(0),
  endo_colnames = character(0),
  instruments = character(0),
  dropped_regressors = character(0),
  dropped_instruments = character(0),
  reclassified_endogenous = character(0),
  method = "ols",
  lambda = NA_real_,
  kclass_value = NA_real_,
  fuller_parameter = 0,
  coviv = FALSE,
  S = NULL,
  W = NULL,
  wmatrix = NULL,
  smatrix = NULL,
  b0 = NULL,
  noid = FALSE,
  cue_convergence = NULL,
  cue_message = NULL,
  kernel = NULL,
  bw = NULL,
  tvar = NULL,
  kiefer = FALSE,
  dkraay = NULL,
  ivar = NULL,
  has_ts_ops = FALSE,
  center = FALSE,
  psd = NULL,
  sw = FALSE,
  partial_ct = 0L,
  partial_names = character(0),
  partialcons = FALSE,
  yy = NULL,
  yyc = NULL,
  rankzz = NULL,
  condxx = NULL,
  condzz = NULL,
  ll = NULL,
  contrasts = NULL,
  xlevels = NULL,
  model = NULL,
  x = NULL,
  y = NULL
)
```

## Arguments

- coefficients:

  Named numeric vector of coefficient estimates.

- residuals:

  Numeric vector of residuals.

- fitted.values:

  Numeric vector of fitted values.

- vcov:

  Variance-covariance matrix of coefficients.

- sigma:

  Residual standard error (scalar).

- df.residual:

  Residual degrees of freedom (integer).

- rank:

  Rank of the model matrix (integer).

- r.squared:

  R-squared.

- adj.r.squared:

  Adjusted R-squared.

- rss:

  Residual sum of squares (scalar).

- r2u:

  Uncentered R-squared (always `1 - rss / sum(w * y^2)`).

- r2c:

  Centered R-squared (always `1 - rss / sum(w * (y - wmean)^2)`).

- mss:

  Model sum of squares (`tss - rss`).

- model_f:

  Model F-statistic.

- model_f_p:

  p-value for model F-test.

- model_f_df1:

  Numerator df for model F-test.

- model_f_df2:

  Denominator df for model F-test.

- diagnostics:

  List of diagnostic test results (NULL for OLS).

- first_stage:

  List of first-stage results (NULL for OLS).

- reduced_form:

  List of reduced-form regression results (NULL if not requested or for
  OLS). See
  [`ivreg2()`](https://restatr.com/ivreg2r/reference/ivreg2.md)
  `reduced_form` argument.

- call:

  The original function call.

- formula:

  The parsed Formula object.

- terms:

  List of terms objects.

- nobs:

  Number of observations (integer).

- vcov_type:

  Character: `"iid"`, `"robust"`, `"HAC"`, `"AC"`, or `"CL"`.

- small:

  Logical: whether small-sample corrections were applied.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

- cluster_var:

  Name of cluster variable(s) (character scalar for one-way, character
  vector of length 2 for two-way, or NULL).

- n_clusters:

  Number of clusters: effective count used for small-sample corrections
  (min(M1, M2) for two-way), or NULL.

- n_clusters1:

  Number of clusters in first dimension (two-way only), or NULL.

- n_clusters2:

  Number of clusters in second dimension (two-way only), or NULL.

- na.action:

  Information about removed observations.

- weights:

  Weights used (or NULL).

- endogenous:

  Character vector of endogenous variable names.

- instruments:

  Character vector of excluded instrument names.

- dropped_regressors:

  Character vector of regressor names dropped due to collinearity (does
  not include reclassified endogenous variables).

- dropped_instruments:

  Character vector of instrument names dropped due to collinearity.

- reclassified_endogenous:

  Character vector of endogenous variable names reclassified as
  exogenous because they were collinear with the instruments.

- method:

  Character: estimation method (`"ols"`, `"2sls"`, `"liml"`, or
  `"kclass"`).

- lambda:

  Numeric: LIML eigenvalue (NA for OLS/2SLS/kclass).

- kclass_value:

  Numeric: the k value actually used (NA for OLS/2SLS).

- fuller_parameter:

  Numeric: Fuller modification parameter (0 if none).

- coviv:

  Logical: whether COVIV (2SLS bread) was used for VCV.

- S:

  Estimated moment-condition covariance matrix (Stata `e(S)`), with
  instrument-named dimnames.

- W:

  GMM weighting matrix (Stata `e(W)`), or NULL for LIML/k-class.

- wmatrix:

  User-supplied weighting matrix (or NULL).

- smatrix:

  User-supplied moment covariance matrix (or NULL).

- b0:

  User-supplied parameter vector for CUE evaluation (or NULL).

- cue_convergence:

  Integer: CUE optimizer convergence code (or NULL).

- cue_message:

  Character: CUE optimizer message (or NULL).

- kernel:

  Canonical kernel name (character) or NULL if no HAC/AC.

- bw:

  Numeric bandwidth or NULL if no HAC/AC.

- tvar:

  Character: name of time variable, or NULL.

- kiefer:

  Logical: whether Kiefer VCE was used.

- dkraay:

  Numeric Driscoll-Kraay bandwidth, or NULL.

- ivar:

  Character: name of panel variable, or NULL.

- contrasts:

  List of contrasts used for factor variables (or NULL).

- xlevels:

  List of factor levels (or NULL).

- model:

  Model frame (or NULL if `model = FALSE`).

- x:

  List with `X` and `Z` matrices (or NULL if `x = FALSE`).

- y:

  Response vector (or NULL if `y = FALSE`).

## Value

An object of class `"ivreg2"`.
