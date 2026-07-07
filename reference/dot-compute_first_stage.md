# First-stage regression diagnostics

For each endogenous variable, fits the first-stage OLS regression on the
full instrument set Z, then computes: F-test of excluded instruments,
partial R², Shea partial R², and RMSE.

## Usage

``` r
.compute_first_stage(
  X,
  Z,
  weights,
  cluster_vec,
  vcov_type,
  endo_names,
  excluded_names,
  N,
  K,
  L,
  K1,
  L1,
  M,
  bread_2sls,
  dofminus = 0L,
  sdofminus = 0L,
  weight_type = "aweight",
  kernel = NULL,
  bw = NULL,
  time_index = NULL,
  center = FALSE,
  sw = FALSE,
  ivar_vec = NULL,
  psd = NULL
)
```

## Arguments

- X:

  N x K regressor matrix.

- Z:

  N x L instrument matrix.

- weights:

  Normalized weights (sum to N), or NULL.

- cluster_vec:

  Cluster membership vector, or NULL.

- vcov_type:

  Character: `"iid"`, `"robust"`, `"HAC"`, `"AC"`, or `"CL"`.

- endo_names:

  Character vector of endogenous variable names.

- excluded_names:

  Character vector of excluded instrument names.

- N, K, L, K1, L1:

  Integer dimensions.

- M:

  Number of clusters (or NULL).

- bread_2sls:

  K x K bread from 2SLS: \\(X\_{hat}'WX\_{hat})^{-1}\\.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

## Value

Named list keyed by endogenous variable names, each element containing:
`f_stat`, `f_p`, `f_df1`, `f_df2`, `partial_r2`, `shea_partial_r2`,
`rmse`, `coefficients`, `residuals`, `fitted.values`, plus SW/AP
diagnostics: `sw_f`, `sw_f_p`, `sw_f_df1`, `sw_f_df2`, `sw_chi2`,
`sw_chi2_p`, `sw_partial_r2`, `ap_f`, `ap_f_p`, `ap_f_df1`, `ap_f_df2`,
`ap_chi2`, `ap_chi2_p`, `ap_partial_r2`.

## Details

First-stage F-statistics always use small-sample degrees of freedom
(F(L1, N-L) or F(L1, M-1) for clustered models), regardless of the main
model's `small` option. This matches Stata's `ivreg2` behavior.
