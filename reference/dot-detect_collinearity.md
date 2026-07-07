# Detect and drop collinear columns via QR decomposition

Detect and drop collinear columns via QR decomposition

## Usage

``` r
.detect_collinearity(mat, label = "column")
```

## Arguments

- mat:

  A numeric matrix.

- label:

  `"regressor"` or `"instrument"` (for messaging).

## Value

A list with `matrix` (cleaned), `dropped` (character names), `rank`.
