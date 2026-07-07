# Compute heteroskedasticity-consistent (robust) VCV

Sandwich estimator using the score-based approach: meat = X' diag(e^2)
X, where X is `X_hat` for IV or `X` for OLS. The base normalization is
`N/(N-dofminus)`. When `small = TRUE`, applies the additional
finite-sample correction `(N-dofminus)/(N-K-dofminus-sdofminus)`, giving
a net scaling of `N/(N-K-dofminus-sdofminus)`.

## Usage

``` r
.compute_hc_vcov(
  bread,
  X_hat,
  resid,
  N,
  K,
  small = FALSE,
  dofminus = 0L,
  sdofminus = 0L,
  weights = NULL,
  weight_type = "aweight",
  center = FALSE
)
```

## Arguments

- bread:

  K x K bread matrix: \\(X'P_Z X)^{-1}\\ for IV, \\(X'X)^{-1}\\ for OLS.

- X_hat:

  N x K matrix: projected regressors for IV, original X for OLS.

- resid:

  Length-N residual vector.

- N:

  Integer: number of observations.

- K:

  Integer: number of regressors.

- small:

  Logical: if `TRUE`, apply the finite-sample correction
  `(N-dofminus)/(N-K-dofminus-sdofminus)`.

- dofminus:

  Integer: large-sample DoF adjustment (default 0).

- sdofminus:

  Integer: small-sample DoF adjustment (default 0).

## Value

K x K variance-covariance matrix.
