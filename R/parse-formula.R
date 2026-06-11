#' @importFrom Formula as.Formula
#' @importFrom stats model.frame model.response model.matrix model.weights
#'   terms as.formula na.omit .getXlevels
NULL

# --------------------------------------------------------------------------
# .parse_formula
# --------------------------------------------------------------------------
#' Parse a three-part IV formula and build model matrices
#'
#' @param formula A formula: `y ~ exog | endo | instruments` (3-part IV) or
#'   `y ~ exog` (1-part OLS).
#' @param data A data frame.
#' @param weights Optional weights expression (evaluated in `data`).
#' @param subset Optional subset expression (evaluated in `data`).
#' @param na.action Function for handling `NA`s (default `na.omit`).
#' @return A named list; see Details.
#' @keywords internal
.parse_formula <- function(formula, data, weights = NULL, subset = NULL,
                           na.action = na.omit) {


  # --- 1. Convert & validate ---
  formula <- Formula::as.Formula(formula)
  n_rhs <- .check_formula_parts(formula)

  # --- 2. Build model frame (match.call / eval idiom from ivreg) ---
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$formula <- formula
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # --- 3. Extract response ---
  if (length(formula)[1L] == 0L) {
    stop("Formula must have a response variable on the left-hand side.",
         call. = FALSE)
  }
  y <- model.response(mf, "numeric")
  if (is.null(y)) {
    stop("Formula must have a response variable on the left-hand side.",
         call. = FALSE)
  }
  if (!is.null(dim(y))) {
    stop("Multivariate responses (e.g. cbind()) are not supported. ",
         "The left-hand side must be a single variable.",
         call. = FALSE)
  }
  y_name <- names(mf)[1L]

  # --- 4. Extract weights / na.action from frame ---
  w <- model.weights(mf)
  na_act <- attr(mf, "na.action")

  # --- 5. Model matrices per part ---
  has_intercept <- attr(terms(formula, rhs = 1L), "intercept") == 1L

  # Part 1: exogenous regressors (includes intercept if present)
  mt_exog <- terms(formula, data = data, rhs = 1L)
  exog <- model.matrix(mt_exog, mf)
  exog_names <- .varnames_from_terms(mt_exog)

  # Track exogenous column names and assign (for orthog/partial factor support)
  exog_assign_full_raw <- attr(exog, "assign")
  icept_pos_exog <- which(colnames(exog) == "(Intercept)")
  if (length(icept_pos_exog) > 0L) {
    orig_exog_colnames <- colnames(exog)[-icept_pos_exog]
    exog_assign_raw <- exog_assign_full_raw[-icept_pos_exog]
  } else {
    orig_exog_colnames <- colnames(exog)
    exog_assign_raw <- exog_assign_full_raw
  }

  if (n_rhs == 1L) {
    # --- OLS path ---
    X <- exog
    Z <- NULL
    endo <- NULL
    excluded <- NULL
    endo_names <- character(0L)
    excluded_names <- character(0L)
    endo_colnames <- character(0L)
    excluded_colnames <- character(0L)
    endo_assign <- integer(0L)
    excluded_assign <- integer(0L)
    endo_ct <- NULL
    excluded_ct <- NULL
    mt_endo <- NULL
    mt_excl <- NULL
  } else {
    # --- IV path ---
    # Part 2: endogenous regressors (strip intercept)
    mt_endo <- terms(formula, data = data, rhs = 2L)
    endo_full <- model.matrix(mt_endo, mf)
    endo_ct <- attr(endo_full, "contrasts")
    endo_assign_full <- attr(endo_full, "assign")
    endo <- .strip_intercept(endo_full)
    # Build assign vector for non-intercept columns
    icept_pos <- which(colnames(endo_full) == "(Intercept)")
    endo_assign <- if (length(icept_pos) > 0L) {
      endo_assign_full[-icept_pos]
    } else {
      endo_assign_full
    }
    endo_names <- .varnames_from_terms(mt_endo)
    orig_endo_colnames <- colnames(endo)

    # Part 3: excluded instruments (strip intercept)
    mt_excl <- terms(formula, data = data, rhs = 3L)
    excluded_full <- model.matrix(mt_excl, mf)
    excluded_ct <- attr(excluded_full, "contrasts")
    excluded_assign_full <- attr(excluded_full, "assign")
    excluded <- .strip_intercept(excluded_full)
    icept_pos_excl <- which(colnames(excluded_full) == "(Intercept)")
    excluded_assign <- if (length(icept_pos_excl) > 0L) {
      excluded_assign_full[-icept_pos_excl]
    } else {
      excluded_assign_full
    }
    excluded_names <- .varnames_from_terms(mt_excl)
    orig_excluded_colnames <- colnames(excluded)

    # --- 6. Check duplicates (on original variable names) ---
    .check_duplicates(formula)

    # Compose X and Z
    # Endogenous before exogenous (matching Stata/ivreg convention), but
    # intercept stays first (R model.matrix convention).
    X <- cbind(endo, exog)
    if (has_intercept && "(Intercept)" %in% colnames(X)) {
      ic <- which(colnames(X) == "(Intercept)")
      if (ic != 1L) X <- X[, c(ic, setdiff(seq_len(ncol(X)), ic)), drop = FALSE]
    }
    Z <- cbind(exog, excluded)
  }

  # --- 7. Detect collinearity ---
  # 3-pass algorithm matching Stata's CheckDupsCollin (ivreg2.ado 4092-4183).
  # OLS path: single pass on X. IV path: 3-pass with reclassification.
  dropped_regressors <- character(0L)
  dropped_instruments <- character(0L)
  reclassified_endogenous <- character(0L)

  if (is.null(Z)) {
    # OLS path: single pass
    col_x <- .detect_collinearity(X, "regressor")
    X <- col_x$matrix
    dropped_regressors <- col_x$dropped
  } else {
    # IV path: 3-pass collinearity detection with reclassification.
    # Suppress intermediate warnings; emit consolidated warnings at the end.
    all_dropped_endo <- character(0L)
    all_dropped_exog <- character(0L)
    all_dropped_excluded <- character(0L)

    # Pass 1 — Intra-endogenous collinearity
    col_endo <- suppressWarnings(.detect_collinearity(endo, "endogenous"))
    if (length(col_endo$dropped) > 0L) {
      all_dropped_endo <- c(all_dropped_endo, col_endo$dropped)
      endo <- col_endo$matrix
    }

    # Pass 2 — Cross-list collinearity (endo columns last so QR drops them first)
    combined <- cbind(exog, excluded, endo)
    col_combined <- suppressWarnings(.detect_collinearity(combined, "column"))
    if (length(col_combined$dropped) > 0L) {
      exog_colnames <- colnames(exog)
      excluded_colnames <- colnames(excluded)
      endo_colnames <- colnames(endo)
      for (d in col_combined$dropped) {
        if (d %in% endo_colnames) {
          # Endogenous var collinear with instruments → reclassify as exogenous
          reclassified_endogenous <- c(reclassified_endogenous, d)
          endo <- endo[, colnames(endo) != d, drop = FALSE]
          exog <- cbind(exog, combined[, d, drop = FALSE])
        } else if (d %in% excluded_colnames) {
          all_dropped_excluded <- c(all_dropped_excluded, d)
          excluded <- excluded[, colnames(excluded) != d, drop = FALSE]
        } else if (d %in% exog_colnames) {
          all_dropped_exog <- c(all_dropped_exog, d)
          exog <- exog[, colnames(exog) != d, drop = FALSE]
        }
      }
    }

    # Pass 3 — Re-check after reclassification (only if reclassification
    # occurred). Runs even when no endogenous columns remain: a reclassified
    # variable can itself be collinear with the exogenous regressors and
    # must be dropped from the combined set (endo is a 0-column matrix then,
    # which cbind and the classification loop handle).
    if (length(reclassified_endogenous) > 0L) {
      combined2 <- cbind(exog, excluded, endo)
      col_combined2 <- suppressWarnings(.detect_collinearity(combined2, "column"))
      if (length(col_combined2$dropped) > 0L) {
        for (d in col_combined2$dropped) {
          if (d %in% colnames(endo)) {
            reclassified_endogenous <- c(reclassified_endogenous, d)
            endo <- endo[, colnames(endo) != d, drop = FALSE]
            exog <- cbind(exog, combined2[, d, drop = FALSE])
          } else if (d %in% colnames(excluded)) {
            all_dropped_excluded <- c(all_dropped_excluded, d)
            excluded <- excluded[, colnames(excluded) != d, drop = FALSE]
          } else if (d %in% colnames(exog)) {
            all_dropped_exog <- c(all_dropped_exog, d)
            exog <- exog[, colnames(exog) != d, drop = FALSE]
          }
        }
      }
    }

    # K1 = 0 cases: part 2 was empty (`y ~ exog | 0 | z`, Stata's `(=z)`), or
    # every endogenous variable was reclassified/dropped by collinearity.
    # Either way Stata classifies the model as OLS but KEEPS the surviving
    # excluded instruments as surplus moment conditions (e(model) decided
    # from the post-collinearity endogenous list, ivreg2.ado:2101-2106; the
    # overid J is not gated on the endogenous count, ivreg2.ado:1164). Only
    # when no excluded instruments survive either does the model degrade to
    # plain OLS.
    if (ncol(endo) == 0L) {
      endo <- NULL
      if (ncol(excluded) == 0L) {
        X <- exog
        Z <- NULL
        excluded <- NULL
        if (length(all_dropped_excluded) > 0L) {
          warning("All excluded instruments dropped; model is now plain OLS ",
                  "(no overidentification test).", call. = FALSE)
        }
        n_rhs <- 1L
      } else {
        X <- exog
        Z <- cbind(exog, excluded)
      }
    } else {
      X <- cbind(endo, exog)
      if (has_intercept && "(Intercept)" %in% colnames(X)) {
        ic <- which(colnames(X) == "(Intercept)")
        if (ic != 1L) X <- X[, c(ic, setdiff(seq_len(ncol(X)), ic)), drop = FALSE]
      }
      Z <- cbind(exog, excluded)
    }

    # Consolidated warnings
    if (length(reclassified_endogenous) > 0L) {
      warning("Endogenous variable(s) collinear with instruments. ",
              "Now treated as exogenous: ",
              paste(reclassified_endogenous, collapse = ", "),
              call. = FALSE)
    }
    dropped_regressors <- c(all_dropped_endo, all_dropped_exog)
    dropped_instruments <- all_dropped_excluded
    all_dropped <- c(dropped_regressors, dropped_instruments)
    n_drop_reg <- length(dropped_regressors)
    n_drop_iv <- length(dropped_instruments)
    if (n_drop_reg > 0L) {
      label_plural <- if (n_drop_reg == 1L) "regressor" else "regressors"
      warning("Dropped ", n_drop_reg, " collinear ", label_plural, ": ",
              paste(dropped_regressors, collapse = ", "),
              call. = FALSE)
    }
    if (n_drop_iv > 0L) {
      label_plural <- if (n_drop_iv == 1L) "instrument" else "instruments"
      warning("Dropped ", n_drop_iv, " collinear ", label_plural, ": ",
              paste(dropped_instruments, collapse = ", "),
              call. = FALSE)
    }
  }

  # Update name vectors after collinearity detection. The endogenous and
  # excluded updates are independent: a K1 = 0 (HOLS) model has endo == NULL
  # but keeps its excluded instruments.
  if (!is.null(endo)) {
    # Term-level names: use assign to map surviving columns → terms
    endo_colnames <- colnames(endo)
    surviving_endo_col_idx <- match(endo_colnames, orig_endo_colnames)
    surviving_endo_terms <- unique(endo_assign[surviving_endo_col_idx])
    endo_names <- attr(mt_endo, "term.labels")[surviving_endo_terms]
    # Re-index assign to map columns → positions in endo_names
    endo_assign <- match(endo_assign[surviving_endo_col_idx],
                         surviving_endo_terms)
  } else {
    # No endogenous regressors (1-part OLS, K1 = 0 form, or all
    # dropped/reclassified)
    endo_names <- character(0L)
    endo_colnames <- character(0L)
    endo_assign <- integer(0L)
  }
  if (!is.null(excluded)) {
    excluded_colnames <- colnames(excluded)
    surviving_excl_col_idx <- match(excluded_colnames, orig_excluded_colnames)
    surviving_excl_terms <- unique(excluded_assign[surviving_excl_col_idx])
    excluded_names <- attr(mt_excl, "term.labels")[surviving_excl_terms]
    excluded_assign <- match(excluded_assign[surviving_excl_col_idx],
                             surviving_excl_terms)
  } else {
    excluded_names <- character(0L)
    excluded_colnames <- character(0L)
    excluded_assign <- integer(0L)
  }
  exog_names <- setdiff(exog_names, dropped_regressors)
  # Only add reclassified vars that survived all passes (not re-dropped in pass 3)
  surviving_reclassified <- intersect(reclassified_endogenous, colnames(X))
  exog_names <- c(exog_names, surviving_reclassified)

  # Exogenous column names, assign, and term labels (original exog part only,
  # excluding reclassified endogenous).  Parallels endo_colnames/endo_assign
  # and excluded_colnames/excluded_assign.
  surviving_orig_exog <- intersect(orig_exog_colnames, colnames(X))
  if (length(surviving_orig_exog) > 0L) {
    surv_idx_exog <- match(surviving_orig_exog, orig_exog_colnames)
    exog_colnames <- surviving_orig_exog
    surv_assign_exog <- exog_assign_raw[surv_idx_exog]
    surviving_exog_terms <- unique(surv_assign_exog)
    exog_term_labels <- attr(mt_exog, "term.labels")[surviving_exog_terms]
    exog_assign <- match(surv_assign_exog, surviving_exog_terms)
  } else {
    exog_colnames <- character(0L)
    exog_assign <- integer(0L)
    exog_term_labels <- character(0L)
  }

  # Reconcile has_intercept with what survived collinearity detection.
  # Unlike Stata (which passes the constant as a separate flag to its
  # collinearity checker), R's model.matrix() puts (Intercept) into X,
  # so QR pivoting could theoretically drop it.
  has_intercept <- "(Intercept)" %in% colnames(X)

  # --- 8. Dimensions and identification checks ---
  N <- length(y)
  K <- ncol(X)                       # total regressors (after drops)
  if (!is.null(endo)) {
    endo_cols_in_X <- intersect(colnames(endo), colnames(X))
    K1 <- length(endo_cols_in_X)
    K2 <- K - K1   # exogenous regressors incl. intercept
  } else {
    K1 <- 0L
    K2 <- K
  }
  L <- if (is.null(Z)) K else ncol(Z)
  L1 <- L - K2                       # excluded instruments (after drops)

  is_iv <- !is.null(Z)

  if (is_iv) {
    if (L1 < K1) {
      stop("Model is underidentified: ", K1, " endogenous regressor(s) but only ",
           L1, " excluded instrument(s) after dropping collinear columns.",
           call. = FALSE)
    }
    if (L1 == K1 && length(dropped_instruments) > 0L) {
      warning("Dropping collinear instruments reduced the model to exact identification. ",
              "Overidentification test will not be available.",
              call. = FALSE)
    }
  }

  is_overid <- is_iv && (L1 > K1)
  overid_df <- if (is_overid) L1 - K1 else 0L

  # --- 9. Build terms objects ---
  # When reclassification or drops changed the variable lists, rebuild terms
  # from the surviving names rather than from the original formula parts.
  all_reg_labels <- c(endo_names, exog_names)
  mt_regressors <- if (length(all_reg_labels) == 0L) {
    # Intercept-only model
    mt_exog
  } else {
    reg_fml <- if (has_intercept) {
      as.formula(paste("~", paste(all_reg_labels, collapse = " + ")))
    } else {
      as.formula(paste("~ 0 +", paste(all_reg_labels, collapse = " + ")))
    }
    terms(reg_fml, data = mf)
  }

  # --- 10. Extract contrasts and xlevels (needed by predict with newdata) ---
  ct <- c(attr(exog, "contrasts"), endo_ct, excluded_ct)
  xl <- .getXlevels(mt_regressors, mf)

  # --- 11. Assemble return list ---
  structure(
    list(
      y              = y,
      X              = X,
      Z              = Z,
      y_name         = y_name,
      exog_names     = exog_names,
      endo_names     = endo_names,
      excluded_names = excluded_names,
      endo_colnames  = endo_colnames,
      excluded_colnames = excluded_colnames,
      exog_colnames  = exog_colnames,
      endo_assign    = endo_assign,
      excluded_assign = excluded_assign,
      exog_assign    = exog_assign,
      exog_term_labels = exog_term_labels,
      X_names        = colnames(X),
      Z_names        = if (!is.null(Z)) colnames(Z) else NULL,
      N              = N,
      K              = K,
      K1             = K1,
      K2             = K2,
      L              = L,
      L1             = L1,
      has_intercept  = has_intercept,
      model_frame    = mf,
      terms          = list(
        regressors  = mt_regressors,
        instruments = if (is_iv) mt_excl else NULL,
        full        = terms(formula)
      ),
      formula            = formula,
      contrasts          = ct,
      xlevels            = xl,
      weights            = w,
      na.action          = na_act,
      dropped_regressors        = dropped_regressors,
      dropped_instruments       = dropped_instruments,
      reclassified_endogenous   = reclassified_endogenous,
      is_iv              = is_iv,
      is_overid          = is_overid,
      overid_df          = overid_df
    ),
    class = "parsed_formula"
  )
}


# --------------------------------------------------------------------------
# .check_formula_parts
# --------------------------------------------------------------------------
#' Validate the number of RHS parts in a Formula object
#'
#' @param formula A `Formula` object.
#' @return Integer: 1 (OLS) or 3 (IV). Stops on invalid counts.
#' @keywords internal
.check_formula_parts <- function(formula) {
  n_rhs <- length(formula)[2L]
  if (n_rhs == 2L) {
    stop("Two-part formula detected. ivreg2() uses a three-part formula: ",
         "y ~ exog | endo | instruments. Did you mean ivreg::ivreg()?",
         call. = FALSE)
  }
  if (n_rhs >= 4L) {
    stop("Formula has ", n_rhs, " parts; ivreg2() expects 1 part (OLS) or ",
         "3 parts: y ~ exog | endo | instruments.",
         call. = FALSE)
  }
  if (n_rhs == 3L) {
    # Empty part 2 (e.g. `y ~ exog | 0 | z`) is the empty-endogenous (HOLS)
    # form, mirroring Stata's `(=z)`: no endogenous regressors, but excluded
    # instruments contributing surplus moment conditions. It requires a
    # non-empty part 3 — with both parts empty there is nothing beyond OLS.
    mt2 <- terms(formula, rhs = 2L)
    vars2 <- attr(mt2, "term.labels")
    if (length(vars2) == 0L) {
      mt3 <- terms(formula, rhs = 3L)
      vars3 <- attr(mt3, "term.labels")
      if (length(vars3) == 0L) {
        stop("Parts 2 (endogenous) and 3 (instruments) are both empty. ",
             "Use a one-part formula for OLS: y ~ x1 + x2.",
             call. = FALSE)
      }
    }
  }
  n_rhs
}


# --------------------------------------------------------------------------
# .check_duplicates
# --------------------------------------------------------------------------
#' Check for duplicate variables across formula parts
#'
#' @param formula A 3-part `Formula` object.
#' @return Invisible `NULL`. Stops with an error naming all duplicates.
#' @keywords internal
.check_duplicates <- function(formula) {
  vars1 <- all.vars(formula(formula, lhs = 0L, rhs = 1L))
  vars2 <- all.vars(formula(formula, lhs = 0L, rhs = 2L))
  vars3 <- all.vars(formula(formula, lhs = 0L, rhs = 3L))

  dup_12 <- intersect(vars1, vars2)
  dup_13 <- intersect(vars1, vars3)
  dup_23 <- intersect(vars2, vars3)
  dupes <- unique(c(dup_12, dup_13, dup_23))

  if (length(dupes) > 0L) {
    stop("Variable(s) appear in multiple formula parts: ",
         paste(dupes, collapse = ", "), ".\n",
         "Each variable must appear in exactly one part of: ",
         "y ~ exog | endo | instruments.",
         call. = FALSE)
  }
  invisible(NULL)
}


# --------------------------------------------------------------------------
# .detect_collinearity
# --------------------------------------------------------------------------
#' Detect and drop collinear columns via QR decomposition
#'
#' @param mat A numeric matrix.
#' @param label `"regressor"` or `"instrument"` (for messaging).
#' @return A list with `matrix` (cleaned), `dropped` (character names), `rank`.
#' @keywords internal
.detect_collinearity <- function(mat, label = "column") {
  if (ncol(mat) == 0L) {
    return(list(matrix = mat, dropped = character(0L), rank = 0L))
  }
  qr_obj <- qr(mat)
  r <- qr_obj$rank
  p <- ncol(mat)

  if (r < p) {
    keep <- qr_obj$pivot[seq_len(r)]
    drop_idx <- qr_obj$pivot[(r + 1L):p]
    dropped <- colnames(mat)[drop_idx]
    n_drop <- length(dropped)
    label_plural <- if (n_drop == 1L) label else paste0(label, "s")
    warning("Dropped ", n_drop, " collinear ", label_plural, ": ",
            paste(dropped, collapse = ", "),
            call. = FALSE)
    mat <- mat[, keep, drop = FALSE]
  } else {
    dropped <- character(0L)
  }
  list(matrix = mat, dropped = dropped, rank = r)
}


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

#' Strip the intercept column from a model matrix
#' @keywords internal
.strip_intercept <- function(mat) {
  icept <- which(colnames(mat) == "(Intercept)")
  if (length(icept) > 0L) mat[, -icept, drop = FALSE] else mat
}

#' Extract original variable names from a terms object
#' @keywords internal
.varnames_from_terms <- function(mt) {
  attr(mt, "term.labels")
}

#' Map term labels to column names using the assign attribute
#'
#' Expands user-supplied term labels (e.g., `"race"`) to the corresponding
#' model matrix column names (e.g., `c("race2", "race3")`) using the assign
#' attribute that maps each column to its term index.
#'
#' @param term_labels Character vector of term labels to expand.
#' @param all_term_labels Full vector of term labels (from `endo_names` or
#'   `excluded_names`).
#' @param all_colnames Full vector of column names (from `endo_colnames` or
#'   `excluded_colnames`).
#' @param assign Integer vector mapping each column to its term index in
#'   `all_term_labels`.
#' @return Character vector of column names corresponding to the given terms.
#' @keywords internal
.expand_terms_to_colnames <- function(term_labels, all_term_labels,
                                       all_colnames, assign) {
  term_idx <- match(term_labels, all_term_labels)
  col_mask <- assign %in% term_idx
  all_colnames[col_mask]
}
