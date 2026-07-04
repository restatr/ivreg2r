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

# Read a scalar fixture CSV (quantity,value rows) into a named list.
read_scalar_fixture <- function(path) {
  fx <- read.csv(path, strip.white = TRUE)
  out <- as.list(fx$value)
  names(out) <- fx$quantity
  out
}

# Read a coefficient fixture CSV (term,estimate,std_error rows) with
# Stata->R name translation (_cons -> (Intercept)).
read_coef_fixture <- function(path) {
  d <- read.csv(path)
  nms <- ifelse(d$term == "_cons", "(Intercept)", d$term)
  list(
    estimate  = setNames(d$estimate, nms),
    std_error = setNames(d$std_error, nms)
  )
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


# Read a Stata matrix fixture (v1..vk CSV from save_matrix) together with its
# companion column-names CSV (from save_matrix_names), mapping Stata's
# "_cons" to R's "(Intercept)". Returns a named square matrix; align to an R
# fit with M[colnames(fit$S), colnames(fit$S)].
read_stata_matrix <- function(mat_path, names_path) {
  M <- as.matrix(read.csv(mat_path))
  nm <- read.csv(names_path)$name
  nm[nm == "_cons"] <- "(Intercept)"
  dimnames(M) <- list(nm, nm)
  M
}


# Shared Stata-to-R ts-operator name translator ("D2.w" -> "d(w, 2)", "L3.unem" -> "l(unem, 3)"), promoted from its two per-file copies at the M-25 review.
translate_stata_ts_names <- function(x) {
  out <- character(length(x))
  for (i in seq_along(x)) {
    s <- x[i]
    if (s %in% c("_cons", "(Intercept)")) {
      out[i] <- "(Intercept)"
      next
    }
    m <- regexec("^([LDld])([0-9]*)\\.(.+)$", s)
    parts <- regmatches(s, m)[[1L]]
    if (length(parts) == 0L) {
      out[i] <- s
      next
    }
    op <- if (toupper(parts[2L]) == "L") "l" else "d"
    k <- if (parts[3L] == "") 1L else as.integer(parts[3L])
    out[i] <- paste0(op, "(", parts[4L], ", ", k, ")")
  }
  out
}

# Shared Stata `xi i.year`-style name translator ("_Iyear_67" -> "factor(year)67", "_cons" -> "(Intercept)"), promoted from its two per-file copies at the M-25 review.
translate_stata_xi_names <- function(x) {
  out <- x
  out[x == "_cons"] <- "(Intercept)"
  m <- regmatches(x, regexec("^_I(.+)_([^_]+)$", x))
  for (i in seq_along(x)) {
    parts <- m[[i]]
    if (length(parts) == 3L) {
      out[i] <- paste0("factor(", parts[2L], ")", parts[3L])
    }
  }
  out
}

# Read a first-stage fixture CSV; identical to read_diagnostics() (read.csv with check.names = FALSE) but kept as a separate domain name for first-stage call sites, promoted from its two per-file copies at the M-25 review.
read_firststage <- function(path) {
  read.csv(path, check.names = FALSE)
}

# Extract a single first-stage statistic for one endogenous regressor from a long-format first-stage fixture (columns: statistic, <endo1>, <endo2>, ...), promoted from its two per-file copies at the M-25 review.
get_fs_value <- function(fixture, stat, endo_name) {
  as.numeric(fixture[fixture$statistic == stat, endo_name])
}

# Deterministic synthetic analytic weight (M-12/M-23 precedent): mod(age,5)+1 matches the .do file's `gen awt = mod(age,5)+1`; shared data object so the formula has one source of truth across the first-stage, redundancy, and edge-case test files.
griliches_awt <- transform(ivreg2r::griliches, awt = age %% 5 + 1)

# Deterministic synthetic analytic weight (M-12/M-23 precedent): mod(year,3)+1 matches the .do file's `gen abwt = mod(year,3)+1`; shared data object so the formula has one source of truth across the first-stage and diagnostics test files.
abdata_awt <- transform(ivreg2r::abdata, abwt = year %% 3 + 1)

# Deterministic card derivations (M-11): region = sum of k*reg66k (the 9 Census-region dummies form a partition; M=9 cluster replacing the retired binary smsa/smsa66 M=2 anti-pattern) and fwt = mod(age,5)+1 (M-12/M-23 precedent); matches the .do file's gen lines.
card_wt <- transform(
  ivreg2r::card,
  region = reg661 * 1 + reg662 * 2 + reg663 * 3 + reg664 * 4 + reg665 * 5 +
    reg666 * 6 + reg667 * 7 + reg668 * 8 + reg669 * 9,
  fwt = age %% 5 + 1
)

# Run a table of Stata fixture cells through the shared small-invariance harness used by the re-based diagnostic families (M-22 orthog, M-23 redundancy, ...). Each cell is a list(name, fixture, fit_args); the harness makes one test_that per cell that skips if the fixture CSV is missing, reads it once, fits the model with small = FALSE and small = TRUE via do.call, runs `compare(fit, fixture)` on both fits, and then asserts the two fits' diagnostics[[slot]] stat/p/df are directly equal at testthat's default ~1.5e-8 tolerance — the Stata statistics behind these families are small-invariant (verified byte-identical in each family's retired fixtures), so the fixtures do not vary `small` and the direct equality keeps the invariant pinned tightly rather than only transitively through the looser Stata tolerances.
test_stata_fixture_cells <- function(cells, compare, slot, label_prefix) {
  for (cell in cells) {
    test_that(paste0(label_prefix, ": ", cell$name), {
      fixture_file <- fixture_path(cell$fixture)
      skip_if(!file.exists(fixture_file), "fixture not found")
      fixture <- read_diagnostics(fixture_file)

      fit <- do.call(ivreg2, c(cell$fit_args, list(small = FALSE)))
      fit_small <- do.call(ivreg2, c(cell$fit_args, list(small = TRUE)))
      compare(fit, fixture)
      compare(fit_small, fixture)

      d <- fit$diagnostics[[slot]]
      d_small <- fit_small$diagnostics[[slot]]
      expect_equal(d_small$stat, d$stat)
      expect_equal(d_small$p, d$p)
      expect_identical(d_small$df, d$df)
    })
  }
}
