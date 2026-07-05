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

# Canonical Stata->R term translation: ts-operator names then xi factor names (the two regexes are disjoint, so composition is a no-op where inapplicable). Composed at the M-16 review to replace per-call-site chaining.
translate_stata_names <- function(x) translate_stata_xi_names(translate_stata_ts_names(x))

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

# Deterministic synthetic analytic weight for the M-17 CUE weighted cell (M-12 precedent): matches generate-cue-fixtures.do's `gen double swwt = mod(date, 4) + 1`; Stata mod() and R %% agree for the negative early-sample date values, so the two constructions are byte-identical.
stockwatson_swwt <- transform(ivreg2r::stockwatson, swwt = date %% 4 + 1)

# Deterministic card derivations (M-11): region = sum of k*reg66k (the 9 Census-region dummies form a partition; M=9 cluster replacing the retired binary smsa/smsa66 M=2 anti-pattern) and fwt = mod(age,5)+1 (M-12/M-23 precedent); matches the .do file's gen lines.
card_wt <- transform(
  ivreg2r::card,
  region = reg661 * 1 + reg662 * 2 + reg663 * 3 + reg664 * 4 + reg665 * 5 +
    reg666 * 6 + reg667 * 7 + reg668 * 8 + reg669 * 9,
  fwt = age %% 5 + 1
)

# Shared card_fw VCE cell table (M-11): the six diagnostics cells consumed by the cross-family weighted loops in test-model-f.R, test-summary-stats.R, test-stock-wright.R, test-diagnostics-endogeneity.R, and test-diagnostics-anderson-rubin.R; fit card_wt with weights = fwt, weight_type = "fweight".
card_fw_vce_cells <- list(
  list(vcov = "iid",    small = FALSE, clusters = NULL,     suffix = "iid"),
  list(vcov = "iid",    small = TRUE,  clusters = NULL,     suffix = "iid_small"),
  list(vcov = "robust", small = FALSE, clusters = NULL,     suffix = "hc1"),
  list(vcov = "robust", small = TRUE,  clusters = NULL,     suffix = "hc1_small"),
  list(vcov = "iid",    small = FALSE, clusters = ~ region, suffix = "cl"),
  list(vcov = "iid",    small = TRUE,  clusters = ~ region, suffix = "cl_small")
)

# Shared core Stata-parity assertion block (M-11 review): coefficients, SEs, VCV, sigma, RSS, and model F, promoted from its three per-file copies at the M-11 review (test_fweight_config and test_pw_style_config in test-weight-types.R, test_aw_cells_config in test-weights.R).
expect_stata_parity_core <- function(fit, stata_coef, stata_vcov, stata_diag) {
  r_coef <- coef(fit)[names(stata_coef$estimate)]
  expect_equal(r_coef, stata_coef$estimate, tolerance = stata_tol$coef)
  r_se <- sqrt(diag(vcov(fit)))[names(stata_coef$std_error)]
  expect_equal(r_se, stata_coef$std_error, tolerance = stata_tol$se)
  r_vcov <- vcov(fit)[rownames(stata_vcov), colnames(stata_vcov)]
  expect_equal(r_vcov, stata_vcov, tolerance = stata_tol$vcov, ignore_attr = TRUE)
  expect_equal(fit$sigma, stata_diag$rmse, tolerance = stata_tol$coef)
  expect_equal(fit$rss, stata_diag$rss, tolerance = stata_tol$stat)
  expect_equal(fit$model_f, stata_diag$F_stat, tolerance = stata_tol$stat)
}

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

# Compare coefficients and SEs against a Stata coef fixture, matching by translated term name (strict: every fit term must have a matching fixture term and vice versa); translates via translate_stata_names(), promoted from the two per-file copies (test-ts-operators.R and test-liml.R) at the M-15 review, keeping the stricter term-setequal semantics.
expect_coef_fixture <- function(fit, coef_file, tol_coef = stata_tol$coef,
                                tol_se = stata_tol$se) {
  fx <- read.csv(fixture_path(coef_file))
  fx$term_r <- translate_stata_names(fx$term)
  b <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  expect_setequal(fx$term_r, names(b))
  for (i in seq_len(nrow(fx))) {
    nm <- fx$term_r[i]
    expect_equal(unname(b[nm]), fx$estimate[i], tolerance = tol_coef,
                 info = paste(coef_file, "coef", nm))
    expect_equal(unname(se[nm]), fx$std_error[i], tolerance = tol_se,
                 info = paste(coef_file, "se", nm))
  }
}

# Compare the full VCV against a Stata vcov fixture, matching by translated term name (strict: row names must equal the fixture's term set); translates via translate_stata_names(), promoted from the two per-file copies (test-ts-operators.R and test-liml.R) at the M-15 review, keeping the stricter term-setequal semantics.
expect_vcov_fixture <- function(fit, vcov_file, tol = stata_tol$vcov) {
  fx <- read.csv(fixture_path(vcov_file))
  stata_terms <- translate_stata_names(fx$term)
  V_s <- as.matrix(fx[, grep("^vcov_", names(fx)), drop = FALSE])
  dimnames(V_s) <- list(stata_terms, stata_terms)
  V_r <- vcov(fit)
  expect_setequal(rownames(V_r), stata_terms)
  perm <- match(rownames(V_r), stata_terms)
  expect_vcov_match(V_r, V_s[perm, perm], tol = tol, label = vcov_file)
}

# Shared diagnostics-comparison helper for Stata-parity full cells: asserts overid, underid, Cragg-Donald and Kleibergen-Paap weak-ID stats, model F, sigma, rss, and N against a diagnostics fixture, each guarded on the fixture column being present (NA where inapplicable); promoted from test-center.R's local expect_center_diagnostics at the M-18 review, with the dead endog_stat block dropped and the CD assertion restored. The tol_stat/tol_pval/tol_coef arguments (added M-17) default to the standard Stata tolerances so existing call sites are unchanged; the CUE cells pass wider coef-class tolerances (sigma/rss use tol_coef) through them.
expect_diagnostics_fixture <- function(fit, diag_file,
                                       tol_stat = stata_tol$stat,
                                       tol_pval = stata_tol$pval,
                                       tol_coef = stata_tol$coef) {
  dx <- read_diagnostics(fixture_path(diag_file))
  d <- fit$diagnostics

  if (!is.na(dx$overid_stat)) {
    expect_equal(d$overid$stat, dx$overid_stat,
                 tolerance = tol_stat, info = "overid stat")
    expect_equal(d$overid$p, dx$overid_p,
                 tolerance = tol_pval, info = "overid p")
    expect_identical(d$overid$df, as.integer(dx$overid_df))
  }
  if (!is.na(dx$underid_stat)) {
    expect_equal(d$underid$stat, dx$underid_stat,
                 tolerance = tol_stat, info = "underid stat")
    expect_equal(d$underid$p, dx$underid_p,
                 tolerance = tol_pval, info = "underid p")
  }
  if (!is.na(dx$weak_id_cd_f)) {
    expect_false(is.null(d$weak_id))
    expect_equal(d$weak_id$stat, dx$weak_id_cd_f,
                 tolerance = tol_stat, info = "CD weak-ID F")
  }
  if (!is.na(dx$weak_id_kp_f)) {
    expect_false(is.null(d$weak_id_robust))
    expect_equal(d$weak_id_robust$stat, dx$weak_id_kp_f,
                 tolerance = tol_stat, info = "KP widstat")
  }
  if (!is.na(dx$model_f)) {
    expect_equal(fit$model_f, dx$model_f,
                 tolerance = tol_stat, info = "model F")
  }
  if (!is.na(dx$sigma)) {
    expect_equal(fit$sigma, dx$sigma,
                 tolerance = tol_coef, info = "sigma")
  }
  if (!is.na(dx$rss)) {
    expect_equal(fit$rss, dx$rss,
                 tolerance = tol_coef, info = "rss")
  }
  expect_identical(nobs(fit), as.integer(dx$N))
}

# Compare orthog() diagnostics against a Stata cstat/cstatp/cstatdf fixture; promoted from test-diagnostics-orthog.R's local compare_orthog() at the M-18 review because test-center.R needs the same comparator for its orthog-under-center cell.
compare_orthog_fixture <- function(fit, fixture) {
  orth <- fit$diagnostics$orthog
  expect_false(is.null(orth))
  if (fixture$cstat == 0) {
    # Underidentified restricted model: Stata reports cstat=0 with missing
    # p/df; R returns stat=0 with df=0
    expect_equal(orth$stat, 0)
    expect_identical(orth$df, 0L)
  } else {
    expect_equal(orth$stat, fixture$cstat, tolerance = stata_tol$stat)
    expect_equal(orth$p, fixture$cstatp, tolerance = stata_tol$pval)
    expect_identical(orth$df, as.integer(fixture$cstatdf))
  }
}

# Shared klein/abdata fixture formulas (one source of truth, M-25 promotion precedent): klein_formula anchors help-file case H72, ab_formula anchors H88.
klein_formula <- consump ~ l(profits, 1) | profits + wagetot |
  govt + taxnetx + year + wagegovt + capital1 + l(totinc, 1)
ab_formula <- n ~ 1 | w + k + ys |
  d(w, 1) + d(k, 1) + d(ys, 1) + d(w, 2) + d(k, 2) + d(ys, 2)

# Shared griliches/phillips fixture formulas (one source of truth, M-15/M-25 promotion precedent): gril_formula anchors help-file case H06 (help.txt:1154), phil_formula anchors the phillips H83-H85 base.
gril_formula <- lw ~ s + expr + tenure + rns + smsa + factor(year) | iq |
  med + kww + age + mrt
phil_formula <- cinf ~ 1 | unem | l(unem, 1:3)

# Shared stockwatson CUE fixture formula (M-17 re-base): anchors the Baum, Schaffer & Stillman (2007) Stata Journal pp. 476/480 worked CUE arc, ivreg2 dinf (UR = ggdp_2 TBILL_1 ER_1 TBON_1).
sw_formula <- dinf ~ 1 | UR | ggdp_2 + TBILL_1 + ER_1 + TBON_1

# Shared weighted-klein data object mirroring the griliches_awt/abdata_awt pattern; matches generate-liml-fixtures.do's `gen double awt = mod(yr, 5) + 1`.
klein_awt <- transform(ivreg2r::klein, awt = yr %% 5 + 1)
