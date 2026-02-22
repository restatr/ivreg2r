#' @importFrom stats coef vcov residuals fitted nobs formula
#'   qt qnorm pt pnorm delete.response na.pass printCoefmat
NULL

# --------------------------------------------------------------------------
# .new_ivreg2
# --------------------------------------------------------------------------
#' Construct an ivreg2 object
#'
#' @param coefficients Named numeric vector of coefficient estimates.
#' @param residuals Numeric vector of residuals.
#' @param fitted.values Numeric vector of fitted values.
#' @param vcov Variance-covariance matrix of coefficients.
#' @param sigma Residual standard error (scalar).
#' @param df.residual Residual degrees of freedom (integer).
#' @param rank Rank of the model matrix (integer).
#' @param r.squared R-squared.
#' @param adj.r.squared Adjusted R-squared.
#' @param rss Residual sum of squares (scalar).
#' @param r2u Uncentered R-squared (always `1 - rss / sum(w * y^2)`).
#' @param r2c Centered R-squared (always `1 - rss / sum(w * (y - wmean)^2)`).
#' @param mss Model sum of squares (`tss - rss`).
#' @param model_f Model F-statistic (NULL until ticket F1).
#' @param model_f_p p-value for model F-test.
#' @param model_f_df1 Numerator df for model F-test.
#' @param model_f_df2 Denominator df for model F-test.
#' @param diagnostics List of diagnostic test results (NULL for OLS).
#' @param first_stage List of first-stage results (NULL for OLS).
#' @param reduced_form List of reduced-form regression results (NULL if not
#'   requested or for OLS). See [ivreg2()] `reduced_form` argument.
#' @param call The original function call.
#' @param formula The parsed Formula object.
#' @param terms List of terms objects.
#' @param nobs Number of observations (integer).
#' @param vcov_type Character: `"iid"`, `"HC0"`, `"HC1"`, or `"CL"`.
#' @param small Logical: whether small-sample corrections were applied.
#' @param cluster_var Name of cluster variable(s) (character scalar for one-way,
#'   character vector of length 2 for two-way, or NULL).
#' @param n_clusters Number of clusters: effective count used for small-sample
#'   corrections (min(M1, M2) for two-way), or NULL.
#' @param n_clusters1 Number of clusters in first dimension (two-way only), or NULL.
#' @param n_clusters2 Number of clusters in second dimension (two-way only), or NULL.
#' @param na.action Information about removed observations.
#' @param weights Weights used (or NULL).
#' @param endogenous Character vector of endogenous variable names.
#' @param instruments Character vector of excluded instrument names.
#' @param dropped_regressors Character vector of regressor names dropped due to
#'   collinearity (does not include reclassified endogenous variables).
#' @param dropped_instruments Character vector of instrument names dropped due
#'   to collinearity.
#' @param reclassified_endogenous Character vector of endogenous variable names
#'   reclassified as exogenous because they were collinear with the instruments.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @param method Character: estimation method (`"ols"`, `"2sls"`, `"liml"`,
#'   or `"kclass"`).
#' @param lambda Numeric: LIML eigenvalue (NA for OLS/2SLS/kclass).
#' @param kclass_value Numeric: the k value actually used (NA for OLS/2SLS).
#' @param fuller_parameter Numeric: Fuller modification parameter (0 if none).
#' @param coviv Logical: whether COVIV (2SLS bread) was used for VCV.
#' @param kernel Canonical kernel name (character) or NULL if no HAC/AC.
#' @param bw Numeric bandwidth or NULL if no HAC/AC.
#' @param tvar Character: name of time variable, or NULL.
#' @param kiefer Logical: whether Kiefer VCE was used.
#' @param dkraay Numeric Driscoll-Kraay bandwidth, or NULL.
#' @param ivar Character: name of panel variable, or NULL.
#' @param contrasts List of contrasts used for factor variables (or NULL).
#' @param xlevels List of factor levels (or NULL).
#' @param model Model frame (or NULL if `model = FALSE`).
#' @param x List with `X` and `Z` matrices (or NULL if `x = FALSE`).
#' @param y Response vector (or NULL if `y = FALSE`).
#' @return An object of class `"ivreg2"`.
#' @keywords internal
.new_ivreg2 <- function(coefficients, residuals, fitted.values, vcov, sigma,
                         df.residual, rank, r.squared, adj.r.squared, rss,
                         r2u, r2c, mss,
                         model_f = NULL, model_f_p = NULL,
                         model_f_df1 = NULL, model_f_df2 = NULL,
                         diagnostics = NULL, first_stage = NULL,
                         reduced_form = NULL,
                         call, formula, terms, nobs, vcov_type, small,
                         dofminus = 0L, sdofminus = 0L,
                         cluster_var = NULL, n_clusters = NULL,
                         n_clusters1 = NULL, n_clusters2 = NULL,
                         na.action = NULL, weights = NULL,
                         weight_type = "aweight", n_physical = NULL,
                         endogenous = character(0),
                         endo_colnames = character(0),
                         instruments = character(0),
                         dropped_regressors = character(0),
                         dropped_instruments = character(0),
                         reclassified_endogenous = character(0),
                         method = "ols",
                         lambda = NA_real_,
                         kclass_value = NA_real_,
                         fuller_parameter = 0,
                         coviv = FALSE,
                         kernel = NULL, bw = NULL,
                         tvar = NULL,
                         kiefer = FALSE, dkraay = NULL,
                         ivar = NULL,
                         contrasts = NULL, xlevels = NULL,
                         model = NULL, x = NULL, y = NULL) {
  structure(
    list(
      coefficients   = coefficients,
      residuals      = residuals,
      fitted.values  = fitted.values,
      vcov           = vcov,
      sigma          = sigma,
      df.residual    = df.residual,
      rank           = rank,
      r.squared      = r.squared,
      adj.r.squared  = adj.r.squared,
      rss            = rss,
      r2u            = r2u,
      r2c            = r2c,
      mss            = mss,
      model_f        = model_f,
      model_f_p      = model_f_p,
      model_f_df1    = model_f_df1,
      model_f_df2    = model_f_df2,
      diagnostics    = diagnostics,
      first_stage    = first_stage,
      reduced_form   = reduced_form,
      call           = call,
      formula        = formula,
      terms          = terms,
      nobs           = nobs,
      vcov_type      = vcov_type,
      small          = small,
      dofminus       = dofminus,
      sdofminus      = sdofminus,
      cluster_var    = cluster_var,
      n_clusters     = n_clusters,
      n_clusters1    = n_clusters1,
      n_clusters2    = n_clusters2,
      na.action      = na.action,
      weights        = weights,
      weight_type    = weight_type,
      n_physical     = n_physical,
      endogenous     = endogenous,
      endo_colnames  = endo_colnames,
      instruments    = instruments,
      dropped_regressors      = dropped_regressors,
      dropped_instruments     = dropped_instruments,
      reclassified_endogenous = reclassified_endogenous,
      method         = method,
      lambda         = lambda,
      kclass_value   = kclass_value,
      fuller_parameter = fuller_parameter,
      coviv          = coviv,
      kernel         = kernel,
      bw             = bw,
      tvar           = tvar,
      kiefer         = kiefer,
      dkraay         = dkraay,
      ivar           = ivar,
      contrasts      = contrasts,
      xlevels        = xlevels,
      model          = model,
      x              = x,
      y              = y
    ),
    class = "ivreg2"
  )
}


# --------------------------------------------------------------------------
# print.ivreg2
# --------------------------------------------------------------------------
#' Print an ivreg2 object
#'
#' @param x An object of class `"ivreg2"`.
#' @param digits Minimum number of significant digits to print.
#' @param ... Additional arguments (ignored).
#' @return `x`, invisibly.
#' @export
print.ivreg2 <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  est_type <- .estimation_label(x)
  cat("\n", est_type, "\n\n", sep = "")
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  if (length(x$reclassified_endogenous) > 0L)
    cat("Reclassified as exog:", paste(x$reclassified_endogenous, collapse = ", "), "\n")
  if (length(x$dropped_regressors) > 0L || length(x$dropped_instruments) > 0L) {
    dropped <- c(x$dropped_regressors, x$dropped_instruments)
    cat("Dropped collinear:", paste(dropped, collapse = ", "), "\n")
  }
  cat("\n")
  invisible(x)
}


# --------------------------------------------------------------------------
# coef.ivreg2
# --------------------------------------------------------------------------
#' Extract coefficients from an ivreg2 object
#'
#' @param object An object of class `"ivreg2"`.
#' @param ... Additional arguments (ignored).
#' @return Named numeric vector of coefficient estimates.
#' @export
coef.ivreg2 <- function(object, ...) {
  object$coefficients
}


# --------------------------------------------------------------------------
# vcov.ivreg2
# --------------------------------------------------------------------------
#' Extract variance-covariance matrix from an ivreg2 object
#'
#' @param object An object of class `"ivreg2"`.
#' @param ... Additional arguments (ignored).
#' @return The variance-covariance matrix of the coefficient estimates.
#' @export
vcov.ivreg2 <- function(object, ...) {
  object$vcov
}


# --------------------------------------------------------------------------
# residuals.ivreg2
# --------------------------------------------------------------------------
#' Extract residuals from an ivreg2 object
#'
#' @param object An object of class `"ivreg2"`.
#' @param ... Additional arguments (ignored).
#' @return Numeric vector of residuals.
#' @export
residuals.ivreg2 <- function(object, ...) {
  object$residuals
}


# --------------------------------------------------------------------------
# fitted.ivreg2
# --------------------------------------------------------------------------
#' Extract fitted values from an ivreg2 object
#'
#' @param object An object of class `"ivreg2"`.
#' @param ... Additional arguments (ignored).
#' @return Numeric vector of fitted values.
#' @export
fitted.ivreg2 <- function(object, ...) {
  object$fitted.values
}


# --------------------------------------------------------------------------
# nobs.ivreg2
# --------------------------------------------------------------------------
#' Extract number of observations from an ivreg2 object
#'
#' @param object An object of class `"ivreg2"`.
#' @param ... Additional arguments (ignored).
#' @return Integer: number of observations.
#' @export
nobs.ivreg2 <- function(object, ...) {
  object$nobs
}


# --------------------------------------------------------------------------
# formula.ivreg2
# --------------------------------------------------------------------------
#' Extract formula from an ivreg2 object
#'
#' @param x An object of class `"ivreg2"`.
#' @param ... Additional arguments passed to [Formula::formula.Formula()]
#'   (e.g., `rhs`, `lhs`, `collapse`).
#' @return The original model formula.
#' @export
formula.ivreg2 <- function(x, ...) {
  stats::formula(x$formula, ...)
}


# --------------------------------------------------------------------------
# confint.ivreg2
# --------------------------------------------------------------------------
#' Confidence intervals for ivreg2 coefficients
#'
#' Computes confidence intervals using the t distribution when `small = TRUE`
#' was used at estimation time, and the standard normal otherwise.
#'
#' @param object An object of class `"ivreg2"`.
#' @param parm A specification of which parameters to give intervals for,
#'   either a numeric vector of positions or a character vector of names.
#'   If missing, all parameters are included.
#' @param level The confidence level (default 0.95).
#' @param ... Additional arguments (ignored).
#' @return A matrix with columns for the lower and upper confidence limits.
#' @export
confint.ivreg2 <- function(object, parm, level = 0.95, ...) {
  cf <- coef(object)
  ses <- sqrt(diag(vcov(object)))
  a <- (1 - level) / 2
  crit <- if (object$small) {
    stats::qt(a, df = object$df.residual)
  } else {
    stats::qnorm(a)
  }
  ci <- cbind(cf + crit * ses, cf - crit * ses)
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE,
                       digits = 3), "%")
  colnames(ci) <- pct
  if (!missing(parm)) {
    if (is.character(parm)) parm <- match(parm, names(cf))
    ci <- ci[parm, , drop = FALSE]
  }
  ci
}


# --------------------------------------------------------------------------
# predict.ivreg2
# --------------------------------------------------------------------------
#' Predict from an ivreg2 model
#'
#' @param object An object of class `"ivreg2"`.
#' @param newdata An optional data frame for prediction. If omitted, fitted
#'   values from the original data are returned.
#' @param na.action Function for handling `NA`s in `newdata`.
#' @param ... Additional arguments (ignored).
#' @return Numeric vector of predicted values.
#' @export
predict.ivreg2 <- function(object, newdata, na.action = stats::na.pass, ...) {
  if (missing(newdata)) {
    return(object$fitted.values)
  }
  tt <- object$terms$regressors
  mf <- stats::model.frame(stats::delete.response(tt), newdata,
                            na.action = na.action,
                            xlev = object$xlevels)
  X <- stats::model.matrix(stats::delete.response(tt), mf,
                            contrasts.arg = object$contrasts)
  cf <- coef(object)
  # model.matrix may regenerate columns dropped for collinearity;
  # use only columns that have corresponding coefficients
  common <- intersect(colnames(X), names(cf))
  drop(X[, common, drop = FALSE] %*% cf[common])
}


# --------------------------------------------------------------------------
# summary.ivreg2
# --------------------------------------------------------------------------
#' Summary for ivreg2 objects
#'
#' Builds a coefficient table (estimates, standard errors, test statistics,
#' p-values) and collects model diagnostics for display by
#' [print.summary.ivreg2()].
#'
#' @param object An object of class `"ivreg2"`.
#' @param ... Additional arguments (ignored).
#' @return An object of class `"summary.ivreg2"`.
#' @export
summary.ivreg2 <- function(object, ...) {
  cf <- coef(object)
  se <- sqrt(diag(vcov(object)))
  stat <- cf / se
  if (object$small) {
    pval <- 2 * stats::pt(-abs(stat), df = object$df.residual)
  } else {
    pval <- 2 * stats::pnorm(-abs(stat))
  }
  coef_table <- cbind(Estimate = cf, `Std. Error` = se,
                      stat, `Pr(>|z|)` = pval)
  stat_label <- if (object$small) "t value" else "z value"
  p_label <- if (object$small) "Pr(>|t|)" else "Pr(>|z|)"
  colnames(coef_table)[3:4] <- c(stat_label, p_label)

  structure(
    c(object, list(coef_table = coef_table)),
    class = "summary.ivreg2"
  )
}


# --------------------------------------------------------------------------
# print.summary.ivreg2
# --------------------------------------------------------------------------
#' Print a summary.ivreg2 object
#'
#' Formats and displays the full estimation output, including coefficient
#' table, fit statistics, and IV diagnostic tests.
#'
#' @param x An object of class `"summary.ivreg2"`.
#' @param digits Minimum number of significant digits.
#' @param signif.stars Logical: print significance stars? Default `TRUE`.
#' @param ... Additional arguments passed to [printCoefmat()].
#' @return `x`, invisibly.
#' @export
print.summary.ivreg2 <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars", TRUE),
                                  ...) {
  is_iv <- length(x$endogenous) > 0L
  est_type <- .estimation_label(x)

  # --- Header ---
  cat("\n", est_type, "\n\n", sep = "")
  cat("Call:\n")
  print(x$call)

  # --- Meta ---
  cat("\nObservations:", format(x$nobs, big.mark = ","), "\n")
  cat("VCV type:    ", .vcov_description(x$vcov_type, x$small,
                                         x$kernel, x$bw,
                                         x$kiefer, x$dkraay), "\n")
  if (!is.null(x$n_clusters)) {
    if (!is.null(x$n_clusters1)) {
      # Two-way clustering
      cat("Clusters:     ", format(x$n_clusters1, big.mark = ","),
          " (", x$cluster_var[1L], "), ",
          format(x$n_clusters2, big.mark = ","),
          " (", x$cluster_var[2L], ")\n", sep = "")
    } else {
      cat("Clusters:    ", format(x$n_clusters, big.mark = ","),
          " (", x$cluster_var, ")\n", sep = "")
    }
  }
  if (!is.null(x$weight_type) && x$weight_type != "aweight") {
    cat("Weight type:  ", x$weight_type, "\n")
  }
  if (!is.null(x$dofminus) && x$dofminus > 0L) {
    cat("dofminus:    ", x$dofminus, "\n")
  }
  if (!is.null(x$sdofminus) && x$sdofminus > 0L) {
    cat("sdofminus:   ", x$sdofminus, "\n")
  }
  if (!is.null(x$method) && x$method %in% c("liml", "kclass")) {
    if (!is.na(x$lambda)) {
      cat("lambda:      ", formatC(x$lambda, digits = 6, format = "f"), "\n")
    }
    if (!is.na(x$kclass_value)) {
      cat("kclass:      ", formatC(x$kclass_value, digits = 6, format = "f"),
          "\n")
    }
    if (!is.null(x$fuller_parameter) && x$fuller_parameter > 0) {
      cat("fuller:      ", x$fuller_parameter, "\n")
    }
    if (isTRUE(x$coviv)) cat("coviv:        TRUE\n")
  }

  # --- Coefficient table ---
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coef_table, digits = digits,
                      signif.stars = signif.stars,
                      na.print = "NA", ...)

  # --- Fit statistics ---
  cat("---\n")
  cat("R-squared:     ", formatC(x$r.squared, digits = 4, format = "f"), "\n")
  cat("Adj. R-squared:", formatC(x$adj.r.squared, digits = 4, format = "f"), "\n")
  if (!is.null(x$model_f)) {
    if (x$small && !is.null(x$model_f_df2)) {
      cat("F(", x$model_f_df1, ", ", x$model_f_df2, "):",
          "     ", formatC(x$model_f, digits = 1, format = "f"),
          " (p ", .format_pval(x$model_f_p), ")\n", sep = "")
    } else {
      cat("Wald chi2(", x$model_f_df1, "):",
          "  ", formatC(x$model_f, digits = 1, format = "f"),
          " (p ", .format_pval(x$model_f_p), ")\n", sep = "")
    }
  }
  cat("Root MSE:      ", formatC(x$sigma, digits = 4, format = "f"), "\n")

  # --- Diagnostics (IV only) ---
  if (is_iv && !is.null(x$diagnostics)) {
    cat("\n")
    .print_iv_diagnostics(x, digits)
  }

  # --- First-stage summary (IV only) ---
  if (is_iv && !is.null(x$first_stage)) {
    .print_first_stage_table(x$first_stage, digits)
  }

  # --- Reduced-form (IV only) ---
  if (is_iv && !is.null(x$reduced_form)) {
    .print_reduced_form(x$reduced_form, digits)
  }

  # --- Footer (IV only) ---
  if (is_iv) {
    cat("\nInstrumented:         ", paste(x$endogenous, collapse = ", "), "\n")
    # Included instruments = exogenous regressors (excl intercept)
    incl <- setdiff(names(coef(x)), c(x$endo_colnames, "(Intercept)"))
    if (length(incl) > 0L) {
      cat("Included instruments: ", paste(incl, collapse = ", "), "\n")
    }
    cat("Excluded instruments: ", paste(x$instruments, collapse = ", "), "\n")
  }

  cat("\n")
  invisible(x)
}


# --------------------------------------------------------------------------
# print.summary helpers
# --------------------------------------------------------------------------

#' Generate estimation type label for display
#' @keywords internal
#' @noRd
.estimation_label <- function(x) {
  m <- x$method %||% (if (length(x$endogenous) > 0L) "2sls" else "ols")
  switch(m,
    "ols"    = "OLS Estimation",
    "2sls"   = "2SLS Estimation",
    "liml"   = if (!is.null(x$fuller_parameter) && x$fuller_parameter > 0) {
                 "Fuller LIML Estimation"
               } else {
                 "LIML Estimation"
               },
    "kclass" = "k-class Estimation",
    paste0(toupper(m), " Estimation")
  )
}

#' Format a VCV type description
#' @keywords internal
#' @noRd
.vcov_description <- function(vcov_type, small, kernel = NULL, bw = NULL,
                               kiefer = FALSE, dkraay = NULL) {
  base <- switch(vcov_type,
    "iid" = "Classical (iid)",
    "HC0" = "Robust (HC0)",
    "HC1" = "Robust (HC1)",
    "CL"  = "Cluster-robust",
    "HAC" = "HAC",
    "AC"  = "AC",
    vcov_type
  )
  if (isTRUE(kiefer)) {
    base <- "Autocorrelation-consistent (Kiefer)"
  } else if (!is.null(dkraay)) {
    base <- paste0("Driscoll-Kraay (kernel=", kernel, "; bandwidth=",
                   formatC(bw, format = "g"), ")")
  } else if (!is.null(kernel) && vcov_type == "CL") {
    # Thompson (two-way cluster + kernel)
    base <- paste0("Two-way cluster + kernel (Thompson; kernel=", kernel,
                   "; bandwidth=", formatC(bw, format = "g"), ")")
  } else if (!is.null(kernel) && vcov_type %in% c("HAC", "AC")) {
    base <- paste0(base, " (kernel=", kernel, "; bandwidth=",
                   formatC(bw, format = "g"), ")")
  }
  if (small && !vcov_type %in% c("iid", "AC") && !isTRUE(kiefer)) {
    paste0(base, ", small-sample corrected")
  } else {
    base
  }
}

#' Format a p-value for display
#' @keywords internal
#' @noRd
.format_pval <- function(p) {
  if (is.null(p) || is.na(p)) return("= NA")
  if (p < 2.2e-16) return("< 2.2e-16")
  paste0("= ", formatC(p, digits = 4, format = "f"))
}

#' Print IV diagnostic tests
#' @keywords internal
#' @noRd
.print_iv_diagnostics <- function(x, digits) {
  diag <- x$diagnostics

  # --- Weak identification ---
  if (!is.null(diag$weak_id)) {
    cat("Weak identification test:\n")
    cat("  Cragg-Donald Wald F:          ",
        formatC(diag$weak_id$stat, digits = 2, format = "f"), "\n")
  }
  if (!is.null(diag$weak_id_robust)) {
    cat("  Kleibergen-Paap rk Wald F:    ",
        formatC(diag$weak_id_robust$stat, digits = 2, format = "f"), "\n")
  }
  # Stock-Yogo critical values
  if (!is.null(diag$weak_id_sy)) {
    sy <- diag$weak_id_sy
    if (x$vcov_type != "iid") {
      cat("  (Stock-Yogo critical values are for iid errors)\n")
    }
    # IV size distortion
    iv_size <- sy[sy$type == "IV size", ]
    if (nrow(iv_size) > 0L) {
      cat("  Stock-Yogo critical values (IV size):\n")
      for (i in seq_len(nrow(iv_size))) {
        cat("    ", iv_size$threshold[i], " maximal IV size",
            "     ", formatC(iv_size$critical_value[i], digits = 2,
                             format = "f"), "\n")
      }
    }
    # IV relative bias
    iv_bias <- sy[sy$type == "IV relative bias", ]
    if (nrow(iv_bias) > 0L) {
      cat("  Stock-Yogo critical values (IV relative bias):\n")
      for (i in seq_len(nrow(iv_bias))) {
        cat("    ", iv_bias$threshold[i], " maximal IV relative bias",
            " ", formatC(iv_bias$critical_value[i], digits = 2,
                         format = "f"), "\n")
      }
    }
    # LIML size distortion
    liml_size <- sy[sy$type == "LIML size", ]
    if (nrow(liml_size) > 0L) {
      cat("  Stock-Yogo critical values (LIML size):\n")
      for (i in seq_len(nrow(liml_size))) {
        cat("    ", liml_size$threshold[i], " maximal LIML size",
            "    ", formatC(liml_size$critical_value[i], digits = 2,
                            format = "f"), "\n")
      }
    }
    # Fuller relative bias
    full_rel <- sy[sy$type == "Fuller relative bias", ]
    if (nrow(full_rel) > 0L) {
      cat("  Stock-Yogo critical values (Fuller relative bias):\n")
      for (i in seq_len(nrow(full_rel))) {
        cat("    ", full_rel$threshold[i], " maximal Fuller rel. bias",
            " ", formatC(full_rel$critical_value[i], digits = 2,
                         format = "f"), "\n")
      }
    }
    # Fuller maximum bias
    full_max <- sy[sy$type == "Fuller maximum bias", ]
    if (nrow(full_max) > 0L) {
      cat("  Stock-Yogo critical values (Fuller maximum bias):\n")
      for (i in seq_len(nrow(full_max))) {
        cat("    ", full_max$threshold[i], " Fuller maximum bias",
            "     ", formatC(full_max$critical_value[i], digits = 2,
                             format = "f"), "\n")
      }
    }
    # Fuller warning (matches Stata line 3426)
    if (nrow(full_rel) > 0L || nrow(full_max) > 0L) {
      cat("  NB: Critical values based on Fuller parameter=1\n")
    }
  } else {
    cat("  Stock-Yogo critical values:   <not available>\n")
  }

  # --- Underidentification ---
  if (!is.null(diag$underid)) {
    uid <- diag$underid
    cat("\nUnderidentification test (", uid$test_name, "):\n", sep = "")
    cat("  Chi-sq(", uid$df, ") = ",
        formatC(uid$stat, digits = 2, format = "f"),
        " (p ", .format_pval(uid$p), ")\n", sep = "")
  }

  # --- Weak-instrument-robust inference ---
  ar <- diag$anderson_rubin
  sw <- diag$stock_wright
  if (!is.null(ar) || !is.null(sw)) {
    cat("\nWeak-instrument-robust inference:\n")
    cat("  H0: B1=0 and orthogonality conditions are valid\n")
    if (!is.null(ar)) {
      cat("  Anderson-Rubin Wald F(",
          ar$f_df1, ",", ar$f_df2, ") = ",
          formatC(ar$f_stat, digits = 2, format = "f"),
          " (p ", .format_pval(ar$f_p), ")\n", sep = "")
      cat("  Anderson-Rubin Wald Chi-sq(",
          ar$chi2_df, ") = ",
          formatC(ar$chi2_stat, digits = 2, format = "f"),
          " (p ", .format_pval(ar$chi2_p), ")\n", sep = "")
    }
    if (!is.null(sw) && !is.na(sw$stat)) {
      cat("  Stock-Wright LM S Chi-sq(",
          sw$df, ") = ",
          formatC(sw$stat, digits = 2, format = "f"),
          " (p ", .format_pval(sw$p), ")\n", sep = "")
    }
  }

  # --- Overidentification ---
  if (!is.null(diag$overid)) {
    oid <- diag$overid
    cat("\nOveridentification test (", oid$test_name, "):", sep = "")
    if (oid$df == 0L) {
      cat("  excluded (exactly identified)\n")
    } else {
      cat("\n  Chi-sq(", oid$df, ") = ",
          formatC(oid$stat, digits = 2, format = "f"),
          " (p ", .format_pval(oid$p), ")\n", sep = "")
    }
  }

  # --- AR LIML overidentification (H3) ---
  if (!is.null(diag$anderson_rubin_overid)) {
    aro <- diag$anderson_rubin_overid
    cat("\nAnderson-Rubin overidentification:", sep = "")
    if (aro$df == 0L) {
      cat("  excluded (exactly identified)\n")
    } else {
      cat("\n  LR Chi-sq(", aro$df, ") = ",
          formatC(aro$lr_stat, digits = 3, format = "f"),
          " (p ", .format_pval(aro$lr_p), ")\n", sep = "")
      cat("  Linearized Chi-sq(", aro$df, ") = ",
          formatC(aro$lin_stat, digits = 3, format = "f"),
          " (p ", .format_pval(aro$lin_p), ")\n", sep = "")
    }
  }

  # --- Endogeneity ---
  if (!is.null(diag$endogeneity)) {
    endog <- diag$endogeneity
    cat("\nEndogeneity test:", sep = "")
    if (is.na(endog$stat)) {
      cat("  not computed (rank-deficient S)\n")
    } else if (endog$stat == 0 && endog$df == 0L) {
      cat("  not computed (collinearity in restricted equation)\n")
    } else {
      cat("\n  Chi-sq(", endog$df, ") = ",
          formatC(endog$stat, digits = 2, format = "f"),
          " (p ", .format_pval(endog$p), ")\n", sep = "")
      cat("  Tested: ", paste(endog$tested_vars, collapse = ", "), "\n", sep = "")
    }
  }

  # --- Orthogonality (instrument subset) ---
  if (!is.null(diag$orthog)) {
    orth <- diag$orthog
    cat("\nOrthogonality test (C-stat):", sep = "")
    if (is.na(orth$stat)) {
      cat("  not computed (rank-deficient S)\n")
    } else if (orth$stat == 0 && orth$df == 0L) {
      cat("  not computed (collinearity or underidentification)\n")
    } else {
      cat("\n  Chi-sq(", orth$df, ") = ",
          formatC(orth$stat, digits = 2, format = "f"),
          " (p ", .format_pval(orth$p), ")\n", sep = "")
      cat("  Tested: ", paste(orth$tested_vars, collapse = ", "), "\n", sep = "")
    }
  }
}

#' Print first-stage diagnostics table
#' @keywords internal
#' @noRd
.print_first_stage_table <- function(first_stage, digits) {
  cat("\nFirst-stage diagnostics:\n")

  # Build table rows
  nms <- names(first_stage)
  rows <- lapply(nms, function(nm) {
    fs <- first_stage[[nm]]
    data.frame(
      Endogenous = nm,
      F_stat = fs$f_stat,
      p_value = fs$f_p,
      Partial_R2 = fs$partial_r2,
      Shea_PR2 = fs$shea_partial_r2,
      SW_F = if (!is.null(fs$sw_f)) fs$sw_f else NA_real_,
      AP_F = if (!is.null(fs$ap_f)) fs$ap_f else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  tbl <- do.call(rbind, rows)

  # Format for display
  cat("  ", formatC("Endogenous", width = 14, flag = "-"),
      formatC("F-stat", width = 10),
      formatC("p-value", width = 10),
      formatC("Partial R2", width = 12),
      formatC("Shea PR2", width = 10),
      formatC("SW F", width = 10),
      formatC("AP F", width = 10), "\n", sep = "")
  for (i in seq_len(nrow(tbl))) {
    cat("  ", formatC(tbl$Endogenous[i], width = 14, flag = "-"),
        formatC(tbl$F_stat[i], digits = 2, format = "f", width = 10),
        formatC(tbl$p_value[i], digits = 4, format = "f", width = 10),
        formatC(tbl$Partial_R2[i], digits = 4, format = "f", width = 12),
        .fmt_or_dash(tbl$Shea_PR2[i], digits = 4, width = 10),
        .fmt_or_dash(tbl$SW_F[i], digits = 2, width = 10),
        .fmt_or_dash(tbl$AP_F[i], digits = 2, width = 10), "\n", sep = "")
  }
}

#' Format a number or display a dash for NA
#' @keywords internal
#' @noRd
.fmt_or_dash <- function(x, digits = 2, width = 10) {
  if (is.na(x)) {
    formatC("\u2014", width = width)
  } else {
    formatC(x, digits = digits, format = "f", width = width)
  }
}


#' Print reduced-form regression summary
#' @keywords internal
#' @noRd
.print_reduced_form <- function(rf, digits) {
  if (rf$mode == "rf") {
    cat("\nReduced-form regression (", rf$depvar, " ~ instruments):\n", sep = "")
    cat("  Root MSE:   ", formatC(rf$sigma, digits = 4, format = "f"), "\n")
    if (!is.na(rf$f_stat)) {
      cat("  F(", rf$f_df1, ", ", rf$f_df2, ") = ",
          formatC(rf$f_stat, digits = 2, format = "f"),
          " (p ", .format_pval(rf$f_p), ")\n", sep = "")
    }
    if (!is.na(rf$partial_r2)) {
      cat("  Partial R2: ", formatC(rf$partial_r2, digits = 4, format = "f"),
          "\n")
    }
  } else {
    # System mode
    cat("\nReduced-form system (instruments on all endogenous + depvar):\n")
    eq_names <- rf$depvar
    for (j in seq_along(eq_names)) {
      eq <- rf$equations[[j]]
      cat("  ", eq_names[j], ":",
          " RMSE=", formatC(rf$sigma[j], digits = 4, format = "f"),
          sep = "")
      if (!is.na(eq$f_stat)) {
        cat("  F(", eq$f_df1, ",", eq$f_df2, ")=",
            formatC(eq$f_stat, digits = 2, format = "f"), sep = "")
      }
      cat("\n")
    }
  }
}
