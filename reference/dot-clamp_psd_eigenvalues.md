# Clamp negative eigenvalues of a symmetric matrix to zero

Shared helper for the eigendecomposition-based routines that project a
matrix onto the nearest positive-semidefinite one. Negative eigenvalues
that are only floating-point round-off are clamped to zero silently; the
`non_psd` flag is set only when a negative eigenvalue is too large to be
noise, so the caller can decide whether to warn.

## Usage

``` r
.clamp_psd_eigenvalues(d)
```

## Arguments

- d:

  Numeric vector of eigenvalues (from `eigen(..., symmetric = TRUE)`).

## Value

List with `values` (eigenvalues with negatives clamped to zero),
`non_psd` (logical: a non-noise negative eigenvalue was present), and
`most_negative` (the smallest eigenvalue before clamping, or 0 if none
were negative) for use in a caller's warning message.

## Details

The tolerance is relative to the largest eigenvalue magnitude but
carries an absolute floor at machine epsilon, so a matrix that is
numerically zero throughout (all eigenvalues near the round-off floor)
is treated as noise rather than as a genuinely non-PSD input.
