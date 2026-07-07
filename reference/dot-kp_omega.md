# Kronecker score covariance for Kleibergen-Paap test

Computes the `(K1*L1) x (K1*L1)` score covariance matrix shat0. For HC:
shat0 = crossprod(scores) / N where row i of scores is kron(V_hat_i,
Z1_perp_i) for each observation i. For cluster: same but with rowsum
over clusters first.

## Usage

``` r
.kp_omega(
  Z1_perp,
  V_hat,
  weights,
  cluster_vec,
  N,
  K1,
  L1,
  weight_type = "aweight",
  kernel = NULL,
  bw = NULL,
  time_index = NULL,
  center = FALSE
)
```

## Arguments

- Z1_perp:

  N x L1 partialled excluded instruments.

- V_hat:

  N x K1 matrix (X1_perp for LM, residuals for Wald).

- weights:

  Normalized weights or NULL.

- cluster_vec:

  Cluster vector or NULL.

- N:

  Number of observations.

- K1:

  Number of endogenous regressors.

- L1:

  Number of excluded instruments.

## Value

`(K1*L1) x (K1*L1)` symmetric matrix.
