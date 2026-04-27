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
.stata_to_r <- function(x) ifelse(x == "_cons", "(Intercept)", x)

.rel_err <- function(r_val, stata_val) {
  if (is.na(r_val) || is.na(stata_val)) return(NA_real_)
  if (abs(stata_val) < 1e-15) return(abs(r_val))
  abs((r_val - stata_val) / stata_val)
}

.read_coef_fixture <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$term <- .stata_to_r(d$term)
  d
}

.read_vcov_fixture <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  terms <- .stata_to_r(d$term)
  vcov_cols <- grep("^vcov_", names(d), value = TRUE)
  V <- as.matrix(d[, vcov_cols])
  rownames(V) <- terms
  col_names <- sub("^vcov_", "", vcov_cols)
  colnames(V) <- .stata_to_r(col_names)
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
    F_p      = list(val = fit$model_f_p, cat = "pval")
  )
  if (!is.null(d$overid)) {
    mapping$sargan   <- list(val = d$overid$stat, cat = "stat")
    mapping$sarganp  <- list(val = d$overid$p, cat = "pval")
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

f_justid <- log(wage) ~ exper + expersq + black + south | educ | nearc4
f_overid <- log(wage) ~ exper + expersq + black + south | educ | nearc4 + nearc2

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

# --- LIML ---
.audit_model("liml_overid_iid",
  ivreg2(f_overid, data = card, method = "liml"),
  "card_liml_overid")
.audit_model("liml_overid_hc1",
  ivreg2(f_overid, data = card, method = "liml", vcov = "robust"),
  "card_liml_overid", "hc1")
.audit_model("liml_overid_cl_small",
  ivreg2(f_overid, data = card, method = "liml", clusters = ~smsa66, small = TRUE),
  "card_liml_overid", "cl_small")

# --- Fuller ---
.audit_model("fuller1_overid_iid",
  ivreg2(f_overid, data = card, method = "liml", fuller = 1),
  "card_fuller1_overid")

# --- Weighted ---
.audit_model("2sls_justid_weighted_iid",
  ivreg2(f_justid, data = card, weights = weight),
  "card_just_id_weighted")
.audit_model("2sls_justid_weighted_hc1",
  ivreg2(f_justid, data = card, weights = weight, vcov = "robust"),
  "card_just_id_weighted", "hc1")

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

# --- GMM2S ---
.audit_model("gmm2s_overid_iid",
  ivreg2(f_overid, data = card, method = "gmm2s"),
  "card_overid_gmm2s")
.audit_model("gmm2s_overid_robust",
  ivreg2(f_overid, data = card, method = "gmm2s", vcov = "robust"),
  "card_overid_gmm2s", "robust")
.audit_model("gmm2s_overid_cluster",
  ivreg2(f_overid, data = card, method = "gmm2s", clusters = ~smsa66),
  "card_overid_gmm2s", "cluster")

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
