# Map term labels to column names using the assign attribute

Expands user-supplied term labels (e.g., `"race"`) to the corresponding
model matrix column names (e.g., `c("race2", "race3")`) using the assign
attribute that maps each column to its term index.

## Usage

``` r
.expand_terms_to_colnames(term_labels, all_term_labels, all_colnames, assign)
```

## Arguments

- term_labels:

  Character vector of term labels to expand.

- all_term_labels:

  Full vector of term labels (from `endo_names` or `excluded_names`).

- all_colnames:

  Full vector of column names (from `endo_colnames` or
  `excluded_colnames`).

- assign:

  Integer vector mapping each column to its term index in
  `all_term_labels`.

## Value

Character vector of column names corresponding to the given terms.
