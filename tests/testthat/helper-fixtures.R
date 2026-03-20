# Shared test helpers for reading Stata benchmark fixtures
# Loaded automatically by testthat before test files run.

# Construct path to a Stata benchmark fixture file
fixture_path <- function(...) {
  file.path(testthat::test_path(), "..", "stata-benchmarks", "fixtures", ...)
}

# Read VCV fixture CSV with Stata->R name translation.
# Auto-detects format: if CSV has a `term` column (Group B format with vcov_
# prefixed columns), translates names. Otherwise returns plain matrix.
read_vcov_fixture <- function(path) {
  fixture <- read.csv(path)
  if ("term" %in% names(fixture)) {
    stata_names <- fixture$term
    r_names <- ifelse(stata_names == "_cons", "(Intercept)", stata_names)
    vcov_cols <- grep("^vcov_", names(fixture), value = TRUE)
    V <- as.matrix(fixture[, vcov_cols])
    rownames(V) <- r_names
    col_stata <- sub("^vcov_", "", vcov_cols)
    colnames(V) <- ifelse(col_stata == "_cons", "(Intercept)", col_stata)
    V
  } else {
    as.matrix(fixture)
  }
}

# Read diagnostics fixture CSV
read_diagnostics <- function(path) {
  read.csv(path, check.names = FALSE)
}

# Compare VCV matrices by name (for fixtures with term column)
expect_vcov_equal <- function(V_r, V_stata, tol = stata_tol$vcov) {
  shared <- intersect(rownames(V_r), rownames(V_stata))
  expect_true(length(shared) > 0, info = "No shared VCV terms found")
  for (rn in shared) {
    for (cn in shared) {
      expect_equal(V_r[rn, cn], V_stata[rn, cn],
                   tolerance = tol, info = paste("VCV mismatch:", rn, cn))
    }
  }
}

# Compare VCV matrices by position (for HAC/timeseries fixtures)
expect_vcov_match <- function(V_r, V_stata, tol = stata_tol$vcov, label = "") {
  V_r <- unname(as.matrix(V_r))
  V_stata <- unname(as.matrix(V_stata))
  P <- nrow(V_stata)
  for (i in seq_len(P)) {
    for (j in seq_len(P)) {
      expect_equal(V_r[i, j], V_stata[i, j],
                   tolerance = tol, info = paste0(label, " VCV[", i, ",", j, "]"))
    }
  }
}
