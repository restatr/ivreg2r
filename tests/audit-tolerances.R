# Tolerance Audit Script
# ---------------------
# Measures actual worst-case R-vs-Stata discrepancies across fixture files.
#
# Run from the repo root:
#   Rscript pkg/tests/audit-tolerances.R
#
# Output: a table of worst-case relative errors by category, flagging
# any values within 5x of their tolerance threshold.

suppressPackageStartupMessages({
  devtools::load_all("pkg", quiet = TRUE)
})

fixture_dir <- file.path("pkg", "tests", "stata-benchmarks", "fixtures")
if (!dir.exists(fixture_dir)) {
  fixture_dir <- file.path("tests", "stata-benchmarks", "fixtures")
}
stopifnot(dir.exists(fixture_dir))

# --- Standard tolerances (from helper-tolerances.R) ---
STANDARD_TOL <- list(coef = 1e-6, se = 1e-6, vcov = 1e-6, stat = 1e-4, pval = 1e-4)

# --- Helpers ---

# ts-operator name translator ("L.profits" -> "l(profits, 1)", "D2.w" -> "d(w, 2)"). This standalone script cannot source the testthat helpers, so this is a copy of translate_stata_ts_names() from helper-fixtures.R (same precedent as the card_fwt comment below) -- keep the two in sync.
.translate_ts_names <- function(x) {
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

# Shared Stata `xi i.year`-style name translator ("_Iyear_67" -> "factor(year)67", "_cons" -> "(Intercept)"). This standalone script cannot source the testthat helpers, so this is a copy of translate_stata_xi_names() from helper-fixtures.R -- keep the two in sync.
.translate_xi_names <- function(x) {
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

.rel_err <- function(r_val, stata_val) {
  if (is.na(r_val) || is.na(stata_val)) return(NA_real_)
  if (abs(stata_val) < 1e-15) return(abs(r_val))
  abs((r_val - stata_val) / stata_val)
}

.read_coef_fixture <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$term <- .translate_xi_names(.translate_ts_names(d$term))
  d
}

.read_vcov_fixture <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  terms <- .translate_xi_names(.translate_ts_names(d$term))
  vcov_cols <- grep("^vcov_", names(d), value = TRUE)
  V <- as.matrix(d[, vcov_cols])
  rownames(V) <- terms
  # Column order matches row order (both derive from e(b)'s column order in
  # the .do generator), so reuse the row-derived terms for column names too --
  # the vcov_ prefixed header names have dots replaced by underscores and
  # can't be round-tripped through the ts-operator regex directly.
  colnames(V) <- terms
  V
}

.read_diag_fixture <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, na.strings = ".")
}

# --- Results accumulator ---
results <- data.frame(
  config = character(),
  category = character(),
  quantity = character(),
  stata_val = numeric(),
  r_val = numeric(),
  rel_error = numeric(),
  stringsAsFactors = FALSE
)

.add <- function(config, category, quantity, stata_val, r_val) {
  re <- .rel_err(r_val, stata_val)
  if (is.na(re)) return(invisible())
  results[nrow(results) + 1L, ] <<- list(
    config = config,
    category = category,
    quantity = quantity,
    stata_val = stata_val,
    r_val = r_val,
    rel_error = re
  )
}

.compare_coefs <- function(config, fit, fixture_path) {
  if (!file.exists(fixture_path)) return(invisible())
  fx <- .read_coef_fixture(fixture_path)
  r_coefs <- coef(fit)
  r_ses <- sqrt(diag(vcov(fit)))
  for (i in seq_len(nrow(fx))) {
    nm <- fx$term[i]
    if (nm %in% names(r_coefs)) {
      .add(config, "coef", nm, fx$estimate[i], unname(r_coefs[nm]))
      .add(config, "se", nm, fx$std_error[i], unname(r_ses[nm]))
    }
  }
}

.compare_vcov <- function(config, fit, fixture_path) {
  if (!file.exists(fixture_path)) return(invisible())
  V_stata <- .read_vcov_fixture(fixture_path)
  V_r <- vcov(fit)
  shared <- intersect(rownames(V_r), rownames(V_stata))
  for (rn in shared) {
    for (cn in shared) {
      .add(config, "vcov", paste0(rn, ":", cn), V_stata[rn, cn], V_r[rn, cn])
    }
  }
}

.compare_diagnostics <- function(config, fit, fixture_path) {
  if (!file.exists(fixture_path)) return(invisible())
  fx <- .read_diag_fixture(fixture_path)
  d <- fit$diagnostics

  mapping <- list(
    N        = list(val = fit$nobs, cat = "exact"),
    rss      = list(val = fit$rss, cat = "stat"),
    r2       = list(val = fit$r.squared, cat = "stat"),
    r2_a     = list(val = fit$adj.r.squared, cat = "stat"),
    rmse     = list(val = fit$sigma, cat = "stat"),
    F_stat   = list(val = fit$model_f, cat = "stat"),
    F_p      = list(val = fit$model_f_p, cat = "pval"),
    lambda   = list(val = fit$lambda, cat = "stat"),
    kclass   = list(val = fit$kclass_value, cat = "stat")
  )
  if (!is.null(d$overid)) {
    mapping$sargan   <- list(val = d$overid$stat, cat = "stat")
    mapping$sarganp  <- list(val = d$overid$p, cat = "pval")
    mapping$j        <- list(val = d$overid$stat, cat = "stat")
    mapping$jp       <- list(val = d$overid$p, cat = "pval")
  }
  if (!is.null(d$underid)) {
    mapping$idstat  <- list(val = d$underid$stat, cat = "stat")
    mapping$idp     <- list(val = d$underid$p, cat = "pval")
  }
  if (!is.null(d$weak_id)) {
    mapping$widstat <- list(val = d$weak_id$stat, cat = "stat")
  }
  if (!is.null(d$weak_id_robust)) {
    mapping$widstat <- list(val = d$weak_id_robust$stat, cat = "stat")
  }
  if (!is.null(d$anderson_rubin)) {
    mapping$arf     <- list(val = d$anderson_rubin$f_stat, cat = "stat")
    mapping$arfp    <- list(val = d$anderson_rubin$f_p, cat = "pval")
    mapping$archi2  <- list(val = d$anderson_rubin$chi2_stat, cat = "stat")
    mapping$archi2p <- list(val = d$anderson_rubin$chi2_p, cat = "pval")
  }
  if (!is.null(d$endogeneity)) {
    mapping$estat  <- list(val = d$endogeneity$stat, cat = "stat")
    mapping$estatp <- list(val = d$endogeneity$p, cat = "pval")
  }
  if (!is.null(d$stock_wright)) {
    mapping$sstat  <- list(val = d$stock_wright$stat, cat = "stat")
    mapping$sstatp <- list(val = d$stock_wright$p, cat = "pval")
  }

  for (qty in names(mapping)) {
    if (qty %in% names(fx) && !is.null(mapping[[qty]]$val)) {
      stata_val <- fx[[qty]]
      if (!is.na(stata_val)) {
        .add(config, mapping[[qty]]$cat, qty, stata_val, mapping[[qty]]$val)
      }
    }
  }
}

# Run a model and compare against all three fixture types
.audit_model <- function(label, expr, prefix, vce_suffix = "iid") {
  coef_path <- file.path(fixture_dir, paste0(prefix, "_coef_", vce_suffix, ".csv"))
  vcov_path <- file.path(fixture_dir, paste0(prefix, "_vcov_", vce_suffix, ".csv"))
  diag_path <- file.path(fixture_dir, paste0(prefix, "_diagnostics_", vce_suffix, ".csv"))

  # Skip if no fixtures exist for this config
  if (!file.exists(coef_path) && !file.exists(vcov_path) && !file.exists(diag_path)) {
    cat(sprintf("  %-40s no fixtures found\n", label))
    return(invisible())
  }

  fit <- tryCatch(
    suppressWarnings(expr),
    error = function(e) {
      cat(sprintf("  %-40s ERROR: %s\n", label, e$message))
      NULL
    }
  )
  if (is.null(fit)) return(invisible())

  .compare_coefs(label, fit, coef_path)
  .compare_vcov(label, fit, vcov_path)
  .compare_diagnostics(label, fit, diag_path)

  n <- sum(results$config == label)
  cat(sprintf("  %-40s %d quantities\n", label, n))
}

# =====================================================================
#  MODEL CONFIGURATIONS
# =====================================================================
data(card)
data(klein)
data(abdata)
data(griliches)
data(phillips)

# Fit the stored lwage, not log(wage): card.dta stores lwage as float32, and
# every Stata fixture was generated from that stored variable. Recomputing
# log(wage) in double precision perturbs the depvar at ~1e-8 relative per
# observation, which this audit then reports as spurious parity gaps of up to
# ~5e-6 relative on near-zero coefficients (measured 2026-07-05: the black
# coef gap is 5.5e-6 with log(wage) vs 1.6e-11 with lwage on the identical
# fit). The audit must measure implementation parity on the same data both
# sides saw, exactly as the test files do.
f_justid <- lwage ~ exper + expersq + black + south | educ | nearc4
f_overid <- lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2

# klein/abdata LIML fixture family (M-15 re-base): klein carries the native K1=2 multi-endogenous + panel-lag LIML/Fuller cells; abdata carries the K1=3 multi-endogenous cluster cell. See helper-fixtures.R's card_wt/fwt comment for the same "keep in sync" precedent on the awt formula below.
klein$awt <- klein$yr %% 5 + 1
klein_f <- consump ~ l(profits, 1) | profits + wagetot |
  govt + taxnetx + year + wagegovt + capital1 + l(totinc, 1)
ab_f <- n ~ 1 | w + k + ys | d(w, 1) + d(k, 1) + d(ys, 1) +
  d(w, 2) + d(k, 2) + d(ys, 2)

cat("Running model configurations...\n\n")

# --- 2SLS just-identified ---
.audit_model("2sls_justid_iid",
  ivreg2(f_justid, data = card),
  "card_just_id")
.audit_model("2sls_justid_hc1",
  ivreg2(f_justid, data = card, vcov = "robust"),
  "card_just_id", "hc1")
.audit_model("2sls_justid_hc1_small",
  ivreg2(f_justid, data = card, vcov = "robust", small = TRUE),
  "card_just_id", "hc1_small")

# --- 2SLS overidentified ---
.audit_model("2sls_overid_iid",
  ivreg2(f_overid, data = card),
  "card_overid")
.audit_model("2sls_overid_hc1",
  ivreg2(f_overid, data = card, vcov = "robust"),
  "card_overid", "hc1")
.audit_model("2sls_overid_hc1_small",
  ivreg2(f_overid, data = card, vcov = "robust", small = TRUE),
  "card_overid", "hc1_small")

# --- LIML (klein/abdata, M-15 re-base) ---
.audit_model("liml_klein_iid",
  ivreg2(klein_f, data = klein, tvar = "yr", method = "liml"),
  "tsop_klein", "liml")
.audit_model("liml_klein_hc1",
  ivreg2(klein_f, data = klein, tvar = "yr", method = "liml", vcov = "robust"),
  "klein_liml", "hc1")
.audit_model("liml_ab_cl_small",
  ivreg2(ab_f, data = abdata, tvar = "year", ivar = "id", method = "liml",
         clusters = ~id, small = TRUE),
  "ab_liml", "cl_small")

# --- Fuller (klein) ---
.audit_model("fuller1_klein_iid",
  ivreg2(klein_f, data = klein, tvar = "yr", fuller = 1),
  "tsop_klein", "fuller1")

# --- Weighted LIML (klein) ---
.audit_model("liml_klein_aw_iid",
  ivreg2(klein_f, data = klein, tvar = "yr", weights = awt, method = "liml"),
  "klein_liml_aw", "iid")

# --- Weighted (M-11 re-based cells: fweight fwt = mod(age,5)+1; pweight = genuine NLS sampling weight) ---
# This standalone script cannot source the testthat helpers, so the fwt
# formula is re-derived here; helper-fixtures.R's `card_wt` is the source of
# truth — keep the two in sync.
card_fwt <- transform(card, fwt = age %% 5 + 1)
.audit_model("2sls_overid_fw_iid",
  ivreg2(f_overid, data = card_fwt, weights = fwt, weight_type = "fweight"),
  "card_fw")
.audit_model("2sls_overid_fw_hc1",
  ivreg2(f_overid, data = card_fwt, weights = fwt, weight_type = "fweight",
         vcov = "robust"),
  "card_fw", "hc1")
.audit_model("2sls_overid_pw_hc0",
  ivreg2(f_overid, data = card, weights = weight, weight_type = "pweight",
         vcov = "robust"),
  "card_pw", "hc0")

# --- CUE ---
.audit_model("cue_overid_iid",
  ivreg2(f_overid, data = card, method = "cue"),
  "card_overid_cue", "iid")
.audit_model("cue_overid_robust",
  ivreg2(f_overid, data = card, method = "cue", vcov = "robust"),
  "card_overid_cue", "robust")
.audit_model("cue_overid_cluster",
  ivreg2(f_overid, data = card, method = "cue", clusters = ~smsa66),
  "card_overid_cue", "cluster")
.audit_model("cue_overid_weighted_aw_robust",
  ivreg2(f_overid, data = card, method = "cue", weights = age, vcov = "robust"),
  "card_overid_cue_weighted", "aw_robust")

# --- Center ---
.audit_model("center_hc0",
  ivreg2(f_overid, data = card, vcov = "robust", center = TRUE),
  "card_overid", "center_hc0")
.audit_model("center_hc1_small",
  ivreg2(f_overid, data = card, vcov = "robust", small = TRUE, center = TRUE),
  "card_overid", "center_hc1_small")
.audit_model("center_cue_hc0",
  ivreg2(f_overid, data = card, method = "cue", vcov = "robust", center = TRUE),
  "card_overid", "center_cue_hc0")
.audit_model("center_cue_cl",
  ivreg2(f_overid, data = card, method = "cue", clusters = ~smsa66, center = TRUE),
  "card_overid", "center_cue_cl")

# Note: dofminus fixtures use lwage (pre-computed) and sdofminus=1, which
# can't be easily reproduced from the bundled card dataset. Those configs
# are covered by the test suite but excluded from this audit script.

# --- GMM2S (M-16 re-base): the retired iid config is identity-covered
# (gmm2s+iid is algebraically identical to 2SLS, asserted in test-gmm2s.R),
# and the retired cluster config's ERROR (rank-deficient S under the M=2
# card clusters=~smsa66 anti-pattern) dies together with the retired cell.
# griliches H06 model formula (help.txt:1154); mirrors the inline formula in
# test-helpfile-examples.R's H06 test -- keep the two in sync.
gril_f <- lw ~ s + expr + tenure + rns + smsa + factor(year) | iq |
  med + kww + age + mrt
phil_f <- cinf ~ 1 | unem | l(unem, 1:3)

.audit_model("gmm2s_gril_robust_small",
  ivreg2(gril_f, data = griliches, method = "gmm2s", vcov = "robust",
         small = TRUE),
  "gril_gmm2s", "robust_small")
.audit_model("gmm2s_ab_cl_small",
  ivreg2(ab_f, data = abdata, tvar = "year", ivar = "id", method = "gmm2s",
         clusters = ~id, small = TRUE),
  "ab_gmm2s", "cl_small")
.audit_model("gmm2s_phil_ac",
  ivreg2(phil_f, data = phillips, tvar = "year", method = "gmm2s",
         kernel = "bartlett", bw = 3),
  "phil_gmm2s", "ac_bw3")

# =====================================================================
#  SUMMARY
# =====================================================================
cat("\n")
cat("============================================================\n")
cat("  TOLERANCE AUDIT RESULTS\n")
cat("============================================================\n\n")

if (nrow(results) == 0L) {
  cat("No comparisons recorded.\n")
  quit("no")
}

# Assign tolerance thresholds based on standard values
results$tolerance <- ifelse(
  results$category == "coef", STANDARD_TOL$coef,
  ifelse(results$category == "se", STANDARD_TOL$se,
  ifelse(results$category == "vcov", STANDARD_TOL$vcov,
  ifelse(results$category == "pval", STANDARD_TOL$pval,
  STANDARD_TOL$stat)))
)
results$margin <- results$tolerance / pmax(results$rel_error, .Machine$double.eps)

# Worst case per category
cat("Worst-case relative errors by category (vs STANDARD tolerances):\n\n")
cat(sprintf("  %-10s  %12s  %12s  %8s  %s\n",
            "Category", "Worst RE", "Tolerance", "Margin", "Config / Quantity"))
cat(sprintf("  %-10s  %12s  %12s  %8s  %s\n",
            "--------", "--------", "---------", "------", "----------------"))

for (cat_name in c("coef", "se", "vcov", "stat", "pval")) {
  subset <- results[results$category == cat_name, ]
  if (nrow(subset) == 0L) next
  worst_idx <- which.max(subset$rel_error)
  worst <- subset[worst_idx, ]
  flag <- if (worst$margin < 1) " *** EXCEEDS ***" else ""
  cat(sprintf("  %-10s  %12.3e  %12.3e  %7.1fx  %s / %s%s\n",
              cat_name,
              worst$rel_error,
              worst$tolerance,
              worst$margin,
              worst$config,
              worst$quantity,
              flag))
}

# Flag anything within 5x of tolerance (or exceeding it)
cat("\n")
close_calls <- results[results$margin < 5 & results$category != "exact", ]
if (nrow(close_calls) > 0L) {
  close_calls <- close_calls[order(close_calls$margin), ]
  n_exceed <- sum(close_calls$margin < 1)
  cat(sprintf("Comparisons within 5x of standard tolerance: %d (%d EXCEED)\n\n",
              nrow(close_calls), n_exceed))
  cat(sprintf("  %-35s  %-6s  %-15s  %12s  %12s  %8s\n",
              "Config", "Cat", "Quantity", "Rel Error", "Tolerance", "Margin"))
  cat(sprintf("  %-35s  %-6s  %-15s  %12s  %12s  %8s\n",
              "------", "---", "--------", "---------", "---------", "------"))
  for (i in seq_len(min(nrow(close_calls), 30L))) {
    r <- close_calls[i, ]
    flag <- if (r$margin < 1) " ***" else ""
    cat(sprintf("  %-35s  %-6s  %-15s  %12.3e  %12.3e  %7.1fx%s\n",
                r$config, r$category, r$quantity, r$rel_error, r$tolerance,
                r$margin, flag))
  }
  if (nrow(close_calls) > 30L) {
    cat(sprintf("  ... and %d more\n", nrow(close_calls) - 30L))
  }
} else {
  cat("All comparisons have > 5x margin.\n")
}

cat("\n")
cat(sprintf("Total comparisons: %d\n", nrow(results)))
cat(sprintf("Total configurations: %d\n", length(unique(results$config))))
cat("\n")
