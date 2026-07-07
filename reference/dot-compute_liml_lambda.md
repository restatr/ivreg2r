# Compute the LIML eigenvalue from concentration matrices

Constructs Y = (y, X1) (response + endogenous regressors), computes the
residual concentration matrix QWW = Y'M_Z Y and the restricted
concentration matrix QWW1 (where Z2 = included instruments only), then
returns the minimum eigenvalue of the symmetrized generalized
eigenproblem.

## Usage

``` r
.compute_liml_lambda(y, X, Z, parsed, w)
```

## Arguments

- y:

  Response vector.

- X:

  Full regressor matrix.

- Z:

  Full instrument matrix.

- parsed:

  Parsed formula object (for endo_names, excluded_names).

- w:

  Weight vector (NULL for unweighted).

## Value

Numeric scalar: LIML eigenvalue lambda.
