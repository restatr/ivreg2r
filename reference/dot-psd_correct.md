# Apply PSD correction to a symmetric matrix

Performs eigenvalue-based correction to ensure a matrix is positive
semi-definite. Two modes are available:

- `"psd0"`: set negative eigenvalues to zero (Politis 2007).

- `"psda"`: replace negative eigenvalues with their absolute values
  (Stock & Watson 2008, Remark 8).

## Usage

``` r
.psd_correct(mat, psd)
```

## Arguments

- mat:

  Symmetric matrix to correct.

- psd:

  Character: `"psd0"` or `"psda"`, or `NULL` (no correction).

## Value

Corrected symmetric matrix (or `mat` unchanged if `psd` is `NULL` or no
negative eigenvalues are found).

## Details

When correction is applied, a warning is emitted indicating how many
eigenvalues were corrected and which mode was used.
