# --------------------------------------------------------------------------
# ivreg2 pipeline: model-frame, weights, partialling, clusters, time index
# Called only by ivreg2(); see R/ivreg2.R for the orchestrator.
# --------------------------------------------------------------------------

#' Match a user-supplied S or W matrix to the instrument columns by name.
#'
#' Mirrors Stata's matsort (ivreg2.ado:4480): when the matrix carries
#' dimnames, its rows and columns are selected and reordered to the
#' instrument order; names that do not cover the instrument list are an
#' error (Stata r(198), "supplied matrix columns/rows do not match IV
#' list"). Extra named rows/columns are dropped by the selection, as in
#' Stata. Unnamed matrices are returned unchanged (positional use).
#'
#' @param M User-supplied square matrix.
#' @param inames Character: instrument column names (colnames(Z)).
#' @param arg Argument name for error messages.
#' @return The matrix with rows/columns in instrument order.
#' @noRd
.match_user_matrix <- function(M, inames, arg) {
  cn <- colnames(M)
  rn <- rownames(M)
  if (is.null(cn) && is.null(rn)) {
    return(M)
  }
  # Rows and columns are matched independently by their own names, as
  # Stata's matsort does (each axis selected separately) -- e.g.
  # S[rev(names), names] is valid input. A missing dimension's names
  # default to the other dimension's.
  match_axis <- function(nm, axis) {
    if (anyDuplicated(nm)) {
      stop("`", arg, "` has duplicated ", axis, " names.", call. = FALSE)
    }
    perm <- match(inames, nm)
    if (anyNA(perm)) {
      stop("`", arg, "` ", axis, " names do not match the instrument ",
           "list; missing: ",
           paste0("'", inames[is.na(perm)], "'", collapse = ", "),
           ". Supply the matrix with instrument names (or unnamed, in ",
           "instrument column order).", call. = FALSE)
    }
    perm
  }
  rperm <- match_axis(if (!is.null(rn)) rn else cn, "row")
  cperm <- match_axis(if (!is.null(cn)) cn else rn, "column")
  M[rperm, cperm, drop = FALSE]
}

#' Stop with a too-few-observations diagnosis.
#'
#' Leads with the observation count relative to the model size (the usual
#' cause). Each of `dofminus`/`sdofminus` is named only when it is actually
#' positive, so a plain one-row fit is not misdiagnosed as a
#' degrees-of-freedom argument being too large, and an unset argument is
#' never presented as if the user supplied it. `sdofminus` may be positive
#' without the user setting it (partialling raises it by the number of
#' partialled-out regressors); the caller passes `context` to say so.
#'
#' @param N Observation count: rows remaining after listwise deletion, or the
#'   weighted count when `weighted = TRUE`.
#' @param count Model dimension being checked (`K` parameters or `L`
#'   instruments).
#' @param dofminus,sdofminus Degrees-of-freedom adjustments.
#' @param count_label `"parameter"` or `"instrument"`.
#' @param weighted Logical: is `N` a sum of weights rather than a row count?
#' @param context Optional sentence appended after the lead diagnosis (for
#'   example, explaining that partialling raised `sdofminus`).
#' @noRd
.stop_too_few_obs <- function(N, count, dofminus, sdofminus, count_label,
                              weighted = FALSE, context = NULL) {
  obs_phrase <- if (weighted) {
    paste0("the weighted observation count (the sum of the weights) is ",
           format(N))
  } else if (N == 1L) {
    "1 observation remains after listwise deletion"
  } else {
    paste0(N, " observations remain after listwise deletion")
  }
  count_phrase <- paste0(count, " ", count_label, if (count == 1L) "" else "s")
  msg <- paste0("Too few observations: ", obs_phrase,
                " but the model has ", count_phrase, ".")
  if (!is.null(context)) {
    msg <- paste0(msg, " ", context)
  }
  if (dofminus > 0L || sdofminus > 0L) {
    named <- c(
      if (dofminus > 0L) paste0("`dofminus` = ", dofminus),
      if (sdofminus > 0L) paste0("`sdofminus` = ", sdofminus)
    )
    letter <- if (identical(count_label, "instrument")) "L" else "K"
    remaining <- N - count - dofminus - sdofminus
    msg <- paste0(msg, " With ", paste(named, collapse = " and "),
                  ", N - ", letter, " - dofminus - sdofminus = ",
                  format(remaining), " (must be > 0).")
  }
  stop(msg, call. = FALSE)
}

#' Stop because a cluster variable contains NA values.
#'
#' One wording shared by the one-way and two-way cluster NA checks.
#' @param var_name Name of the offending cluster variable.
#' @noRd
.stop_cluster_na <- function(var_name) {
  stop("Cluster variable '", var_name, "' contains NA values. ",
       "Remove or handle rows with missing cluster values explicitly ",
       "(e.g. by filtering them out) before calling `ivreg2()`.",
       call. = FALSE)
}

#' Prepare model matrices, weights, clusters, and time-index for estimation.
#'
#' Post-parse validation, weight normalization, FWL partialling, b0 validation,
#' wmatrix/smatrix name matching and dimension checks, cluster parsing,
#' DK auto-clustering, time-index construction, and HAC/AC sorting.
#'
#' @return Named list with parsed (modified), sdofminus, b0, cluster_vec,
#'   cluster_var_name, M/M1/M2, time_index, unsort_order, partial_ct,
#'   partial_names, partialcons, n_physical, w_raw, bw.
#' @noRd
.prepare_model <- function(parsed, opts, data, clusters) {
  method <- opts$method
  vcov <- opts$vcov
  kernel <- opts$kernel
  bw <- opts$bw
  kiefer <- opts$kiefer
  dkraay <- opts$dkraay
  tvar <- opts$tvar
  ivar <- opts$ivar
  dofminus <- opts$dofminus
  sdofminus <- opts$sdofminus
  weight_type <- opts$weight_type
  partial <- opts$partial
  nopartialsmall <- opts$nopartialsmall
  fuller <- opts$fuller
  b0 <- opts$b0
  wmatrix <- opts$wmatrix
  smatrix <- opts$smatrix

  # --- Validate dofminus/sdofminus against model dimensions ---
  if (dofminus >= parsed$N) {
    stop("`dofminus` (", dofminus, ") must be less than N (", parsed$N, ").",
         call. = FALSE)
  }
  if (parsed$N - parsed$K - dofminus - sdofminus <= 0L) {
    .stop_too_few_obs(parsed$N, parsed$K, dofminus, sdofminus, "parameter")
  }
  if (parsed$is_iv && parsed$N - parsed$L - dofminus - sdofminus <= 0L) {
    .stop_too_few_obs(parsed$N, parsed$L, dofminus, sdofminus, "instrument")
  }

  # --- Validate method against parsed model ---
  if (method %in% c("liml", "kclass", "gmm2s", "cue") && !parsed$is_iv) {
    stop('`method = "', method, '"` requires an IV model (3-part formula).',
         call. = FALSE)
  }
  # Empty-endogenous (K1 = 0) models: every method is valid. gmm2s/cue with
  # a non-iid VCE give Cragg's (1983) HOLS (iid weighting collapses to OLS);
  # liml/kclass are well-defined too — with X contained
  # in Z, X'M_Z = 0, so the k-class normal equations reduce to X'X for ANY
  # k and the coefficients equal OLS exactly, while the LIML lambda is the
  # 1x1 eigenvalue RSS(y~X)/RSS(y~Z) feeding the Anderson-Rubin LR overid
  # statistic e(arubin) (Stata runs this; verified vs e(lambda) on mroz).
  if (fuller > 0 && parsed$is_iv && fuller >= (parsed$N - parsed$L)) {
    stop("`fuller` (", fuller, ") must be less than N - L (",
         parsed$N - parsed$L, ").", call. = FALSE)
  }

  # --- Validate wmatrix/smatrix (IV-only; symmetry checked after the
  #     name matching below, mirroring Stata: matsort first, issymmetric
  #     on the sorted matrix, ivreg2.ado:4281-4292) ---
  if (!is.null(wmatrix) && !parsed$is_iv) {
    stop("`wmatrix` requires an IV model (endogenous variables).",
         call. = FALSE)
  }
  if (!is.null(smatrix) && !parsed$is_iv) {
    stop("`smatrix` requires an IV model (endogenous variables).",
         call. = FALSE)
  }

  # --- Validate and normalize weights ---
  w_raw <- parsed$weights
  if (!is.null(w_raw) && any(!is.finite(w_raw)))
    stop("Weights must be finite and non-missing.", call. = FALSE)
  if (!is.null(w_raw) && any(w_raw <= 0))
    stop("All weights must be strictly positive.", call. = FALSE)
  if (weight_type == "fweight" && !is.null(w_raw)) {
    if (any(abs(w_raw - round(w_raw)) > sqrt(.Machine$double.eps)))
      stop('Frequency weights (`weight_type = "fweight"`) must be integers.',
           call. = FALSE)
  }

  # Weight normalization dispatch.
  # Raw (user-supplied) weights are stored in the return object.
  n_physical <- parsed$N
  if (!is.null(w_raw)) {
    if (weight_type == "fweight") {
      # fweight: no normalization; N = sum(w) as integer
      parsed$weights <- w_raw
      N_eff <- sum(w_raw)
      if (N_eff > .Machine$integer.max)
        stop("Sum of frequency weights exceeds integer limit.", call. = FALSE)
      parsed$N <- as.integer(round(N_eff))
    } else if (weight_type == "iweight") {
      # iweight: no normalization; N = sum(w) as FLOAT (Stata v3.0.00+)
      # All intermediate calculations (sigma, R², F) use float N.
      # floor(N) applied only to posted nobs and df_r.
      parsed$weights <- w_raw
      parsed$N <- sum(w_raw)
    } else {
      # aweight/pweight: normalize to sum = N (Stata convention)
      parsed$weights <- w_raw * (parsed$N / sum(w_raw))
    }
  }

  # Re-validate dofminus and fuller against (possibly updated) N for
  # fweight/iweight, which redefine parsed$N after the initial validation.
  if (weight_type %in% c("fweight", "iweight") && !is.null(w_raw)) {
    N_check <- if (weight_type == "iweight") floor(parsed$N) else parsed$N
    if (dofminus >= N_check) {
      stop("`dofminus` (", dofminus, ") must be less than N (", N_check, ").",
           call. = FALSE)
    }
    if (N_check - parsed$K - dofminus - sdofminus <= 0L) {
      .stop_too_few_obs(N_check, parsed$K, dofminus, sdofminus, "parameter",
                        weighted = TRUE)
    }
    if (parsed$is_iv && N_check - parsed$L - dofminus - sdofminus <= 0L) {
      .stop_too_few_obs(N_check, parsed$L, dofminus, sdofminus, "instrument",
                        weighted = TRUE)
    }
    # Re-validate fuller against weighted N (Stata validates inside Mata
    # using the already-weighted N; our pre-weight check at line ~447 uses
    # the physical row count).
    if (fuller > 0 && parsed$is_iv && fuller >= (parsed$N - parsed$L)) {
      stop("`fuller` (", fuller, ") must be less than N - L (",
           parsed$N - parsed$L, ").", call. = FALSE)
    }
  }

  # ===== Design-matrix work (weights, partialling, b0) =====

  # --- Partial out exogenous regressors (FWL projection) ---
  partial_ct <- 0L
  partial_names <- character(0L)
  partialcons <- FALSE

  if (!is.null(partial)) {
    # Time/panel variables cannot appear in the model when partialling
    # (Stata ivreg2.ado:562-577: "cannot use time variable ... in
    # combination with -partial- option").
    model_vars <- all.vars(parsed$formula)
    if (!is.null(tvar) && tvar %in% model_vars) {
      stop("Cannot use time variable '", tvar, "' as dependent variable, ",
           "regressor or instrument in combination with `partial`.",
           call. = FALSE)
    }
    if (!is.null(ivar) && ivar %in% model_vars) {
      stop("Cannot use panel variable '", ivar, "' as dependent variable, ",
           "regressor or instrument in combination with `partial`.",
           call. = FALSE)
    }

    # Resolve special strings
    if (length(partial) == 1L && partial == "_all") {
      partial <- parsed$exog_term_labels
    }

    # Handle _cons / (Intercept)
    cons_strings <- c("_cons", "(Intercept)")
    has_cons_request <- any(partial %in% cons_strings)
    if (has_cons_request) {
      if (!parsed$has_intercept) {
        stop('Cannot partial `"_cons"` from a model without an intercept.',
             call. = FALSE)
      }
      partialcons <- TRUE
      partial <- setdiff(partial, cons_strings)
    }

    # Validate remaining names against exogenous term labels
    if (length(partial) > 0L) {
      bad <- setdiff(partial, parsed$exog_names)
      if (length(bad) > 0L) {
        stop("`partial` contains variables not in the exogenous regressor list: ",
             paste0("'", bad, "'", collapse = ", "), ".", call. = FALSE)
      }
      partial_names <- partial

      # When partialling variables, constant is always included (Stata behavior)
      if (parsed$has_intercept) {
        partialcons <- TRUE
      }
    }

    # Nothing to do if no vars and no cons
    if (length(partial_names) == 0L && !partialcons) {
      partial <- NULL
    }
  }

  if (!is.null(partial) || partialcons) {
    exog_mt <- stats::terms(parsed$formula, data = data, rhs = 1L)
    exog_term_labels <- attr(exog_mt, "term.labels")
    exog_colnames <- setdiff(colnames(parsed$X),
                              c("(Intercept)", parsed$endo_colnames))
    x_all_names <- colnames(parsed$X)
    exog_col_mask <- x_all_names %in% exog_colnames
    exog_full_mm <- stats::model.matrix(exog_mt, parsed$model_frame)
    exog_full_assign <- attr(exog_full_mm, "assign")
    icept_pos <- which(colnames(exog_full_mm) == "(Intercept)")
    if (length(icept_pos) > 0L) {
      exog_full_colnames <- colnames(exog_full_mm)[-icept_pos]
      exog_full_assign <- exog_full_assign[-icept_pos]
    } else {
      exog_full_colnames <- colnames(exog_full_mm)
    }
    surv_idx <- match(exog_colnames, exog_full_colnames)
    surv_idx <- surv_idx[!is.na(surv_idx)]
    exog_assign <- exog_full_assign[surv_idx]

    if (length(partial_names) > 0L) {
      partial_colnames <- .expand_terms_to_colnames(
        partial_names, exog_term_labels, exog_full_colnames, exog_full_assign
      )
      partial_colnames <- intersect(partial_colnames, exog_colnames)
    } else {
      partial_colnames <- character(0L)
    }

    if (partialcons && "(Intercept)" %in% colnames(parsed$X)) {
      partial_colnames <- c("(Intercept)", partial_colnames)
    }

    partial_ct <- length(partial_colnames)
    parsed <- .partial_out(parsed, partial_colnames, partialcons)

    # Stata recounts regressors AFTER partialling and errors if none remain
    # (CheckMisc rhs1_ct check, ivreg2.ado:633/641/4246-4249, exit 102).
    # Reachable when every regressor is exogenous (1-part OLS or
    # empty-endogenous models) and partial = "_all".
    if (parsed$K == 0L) {
      stop("No regressors remain after partialling: all regressors were ",
           "partialled out via `partial`.", call. = FALSE)
    }

    if (!nopartialsmall) {
      sdofminus <- sdofminus + as.integer(partial_ct)
    }

    # Re-validate dimensions after partialling. Partialling raised sdofminus
    # by partial_ct (unless nopartialsmall), so tell the user that is where
    # the degrees of freedom went.
    partial_context <- if (!nopartialsmall && partial_ct > 0L) {
      paste0("The ", partial_ct, " regressor",
             if (partial_ct == 1L) "" else "s",
             " partialled out via `partial` also count against the degrees ",
             "of freedom (via `sdofminus`).")
    }
    if (parsed$N - parsed$K - dofminus - sdofminus <= 0L) {
      .stop_too_few_obs(parsed$N, parsed$K, dofminus, sdofminus, "parameter",
                        context = partial_context)
    }
    if (parsed$is_iv && parsed$N - parsed$L - dofminus - sdofminus <= 0L) {
      .stop_too_few_obs(parsed$N, parsed$L, dofminus, sdofminus, "instrument",
                        context = partial_context)
    }
  }

  # --- Validate b0 dimensions and reorder (after partialling) ---
  if (!is.null(b0)) {
    if (length(b0) != parsed$K) {
      stop("`b0` length (", length(b0), ") must equal the number of ",
           "regressors K (", parsed$K, ").", call. = FALSE)
    }
    if (!is.null(names(b0))) {
      xnames <- colnames(parsed$X)
      bad <- setdiff(names(b0), xnames)
      if (length(bad) > 0L) {
        stop("`b0` has names not matching model columns: ",
             paste0("'", bad, "'", collapse = ", "), ".", call. = FALSE)
      }
      missing_names <- setdiff(xnames, names(b0))
      if (length(missing_names) > 0L) {
        stop("`b0` is missing names for model columns: ",
             paste0("'", missing_names, "'", collapse = ", "), ".",
             call. = FALSE)
      }
      b0 <- b0[xnames]
    } else {
      names(b0) <- colnames(parsed$X)
    }
  }

  # --- Normalize and validate wmatrix/smatrix (after partialling) ---
  # Named matrices are matched to the instrument columns BY NAME, mirroring
  # Stata's matsort (ivreg2.ado:4281, 4307: rows/columns are selected and
  # reordered against cnZ1; missing names error with r(198)). Unnamed
  # matrices are used positionally (an R extension -- Stata matrices always
  # carry names), with the same convention as named/unnamed `b0` above.
  if (!is.null(wmatrix)) {
    wmatrix <- .match_user_matrix(wmatrix, colnames(parsed$Z), "wmatrix")
    if (nrow(wmatrix) != parsed$L || ncol(wmatrix) != parsed$L)
      stop("`wmatrix` dimensions (", nrow(wmatrix), "x", ncol(wmatrix),
           ") do not match the number of instruments (", parsed$L, ").",
           call. = FALSE)
    if (!isSymmetric(unname(wmatrix), tol = sqrt(.Machine$double.eps)))
      stop("`wmatrix` is not symmetric.", call. = FALSE)
  }
  if (!is.null(smatrix)) {
    smatrix <- .match_user_matrix(smatrix, colnames(parsed$Z), "smatrix")
    if (nrow(smatrix) != parsed$L || ncol(smatrix) != parsed$L)
      stop("`smatrix` dimensions (", nrow(smatrix), "x", ncol(smatrix),
           ") do not match the number of instruments (", parsed$L, ").",
           call. = FALSE)
    if (!isSymmetric(unname(smatrix), tol = sqrt(.Machine$double.eps)))
      stop("`smatrix` is not symmetric.", call. = FALSE)
  }

  # ===== Dependence-structure work (clusters, time-index, sorting) =====

  # --- Parse clusters ---
  cluster_vec <- NULL
  cluster_var_name <- NULL
  M <- NULL
  M1 <- NULL
  M2 <- NULL
  if (!is.null(clusters)) {
    cl_terms <- attr(stats::terms(clusters), "term.labels")
    has_interaction <- any(grepl(":", cl_terms, fixed = TRUE))
    if (has_interaction)
      stop("`clusters` must use `+` for two-way clustering (e.g. ~a + b), ",
           "not `:` or `*`. For one-way clustering on an interaction, ",
           "create the variable first (e.g. interaction(a, b)).", call. = FALSE)
    cluster_var_names <- cl_terms
    n_clvars <- length(cluster_var_names)
    if (n_clvars < 1L || n_clvars > 2L)
      stop("`clusters` must reference one or two variables.", call. = FALSE)

    mf_rows <- match(rownames(parsed$model_frame), rownames(data))

    if (n_clvars == 1L) {
      cluster_var_name <- cluster_var_names
      if (!cluster_var_name %in% names(data))
        stop("Cluster variable '", cluster_var_name, "' not found in data.",
             call. = FALSE)
      cluster_vec <- data[[cluster_var_name]][mf_rows]
      if (anyNA(cluster_vec))
        .stop_cluster_na(cluster_var_name)
      M <- length(unique(cluster_vec))
      if (M < 2L)
        stop("At least 2 clusters required; found ", M, ".", call. = FALSE)
    } else {
      cluster_var_name <- cluster_var_names
      for (cvn in cluster_var_name) {
        if (!cvn %in% names(data))
          stop("Cluster variable '", cvn, "' not found in data.",
               call. = FALSE)
      }
      cv1 <- data[[cluster_var_name[1L]]][mf_rows]
      cv2 <- data[[cluster_var_name[2L]]][mf_rows]
      if (anyNA(cv1))
        .stop_cluster_na(cluster_var_name[1L])
      if (anyNA(cv2))
        .stop_cluster_na(cluster_var_name[2L])
      M1 <- length(unique(cv1))
      M2 <- length(unique(cv2))
      if (M1 < 2L)
        stop("At least 2 clusters required in '", cluster_var_name[1L],
             "'; found ", M1, ".", call. = FALSE)
      if (M2 < 2L)
        stop("At least 2 clusters required in '", cluster_var_name[2L],
             "'; found ", M2, ".", call. = FALSE)
      M <- min(M1, M2)
      cluster_vec <- list(cv1, cv2)
    }
  }

  # --- DK auto-clustering on tvar ---
  if (!is.null(dkraay) && is.null(cluster_vec)) {
    mf_rows_dk <- match(rownames(parsed$model_frame), rownames(data))
    cluster_var_name <- tvar
    cluster_vec <- data[[tvar]][mf_rows_dk]
    M <- length(unique(cluster_vec))
    if (M < 2L) {
      stop("At least 2 time periods required for dkraay; found ", M, ".",
           call. = FALSE)
    }
  }

  # --- Validate cluster+kernel combinations ---
  if (!is.null(kernel) && !is.null(cluster_vec)) {
    if (is.list(cluster_vec)) {
      if (is.null(tvar) || is.null(ivar)) {
        stop("cluster+kernel requires both `tvar` and `ivar`.", call. = FALSE)
      }
      cv_names <- cluster_var_name
      cv_set <- sort(cv_names)
      iv_set <- sort(c(ivar, tvar))
      if (!identical(cv_set, iv_set)) {
        stop("cluster+kernel requires cluster variables to match ivar and tvar.",
             call. = FALSE)
      }
      if (cv_names[1L] != ivar) {
        cluster_vec <- list(cluster_vec[[2L]], cluster_vec[[1L]])
        cluster_var_name <- c(ivar, tvar)
        M1_old <- M1
        M2_old <- M2
        M1 <- M2_old
        M2 <- M1_old
      }
    } else {
      if (is.null(tvar)) {
        stop("cluster+kernel requires `tvar`.", call. = FALSE)
      }
      if (!identical(cluster_var_name, tvar)) {
        stop("cluster+kernel requires clustering on the time variable.",
             call. = FALSE)
      }
    }
  }

  time_index <- NULL
  unsort_order <- NULL
  if (!is.null(kernel)) {
    mf_rows_ti <- match(rownames(parsed$model_frame), rownames(data))
    if (!tvar %in% names(data)) {
      stop("Time variable '", tvar, "' not found in data.", call. = FALSE)
    }
    tvar_vec <- data[[tvar]][mf_rows_ti]
    if (anyNA(tvar_vec)) {
      stop("Time variable '", tvar, "' contains NA values.", call. = FALSE)
    }
    if (!is.numeric(tvar_vec)) {
      stop("Time variable '", tvar, "' must be numeric.", call. = FALSE)
    }
    if (any(!is.finite(tvar_vec))) {
      stop("Time variable '", tvar, "' must contain only finite values.",
           call. = FALSE)
    }

    ivar_vec <- NULL
    if (!is.null(ivar)) {
      if (!ivar %in% names(data)) {
        stop("Panel variable '", ivar, "' not found in data.", call. = FALSE)
      }
      ivar_vec <- data[[ivar]][mf_rows_ti]
      if (anyNA(ivar_vec)) {
        stop("Panel variable '", ivar, "' contains NA values.", call. = FALSE)
      }
    }

    time_index <- .build_time_index(tvar_vec, ivar_vec)

    if (isTRUE(kiefer)) {
      bw <- time_index$T_span
    }

    if (time_index$n_gaps > 0L) {
      warning("Time variable '", tvar, "' has ", time_index$n_gaps,
              " gap(s) in relevant range.", call. = FALSE)
    }

    if (is.numeric(bw) && !isTRUE(kiefer)) {
      max_bw <- (time_index$T_span - 1) / time_index$tdelta
      if (bw > max_bw) {
        stop("The bandwidth `bw` (", bw, ") cannot exceed the time span of ",
             "the data (", max_bw, ").", call. = FALSE)
      }
    }

    # Sort all matrices by time-index sort order
    so <- time_index$sort_order
    unsort_order <- order(so)
    parsed$X <- parsed$X[so, , drop = FALSE]
    parsed$Z <- parsed$Z[so, , drop = FALSE]
    parsed$y <- parsed$y[so]
    if (!is.null(parsed$weights)) {
      parsed$weights <- parsed$weights[so]
    }
    if (!is.null(cluster_vec)) {
      if (is.list(cluster_vec)) {
        cluster_vec[[1L]] <- cluster_vec[[1L]][so]
        cluster_vec[[2L]] <- cluster_vec[[2L]][so]
      } else {
        cluster_vec <- cluster_vec[so]
      }
    }
  }

  # --- Gap warning for ts operators without a kernel ---
  # Stata warns about gaps whenever ts operators OR a kernel are present
  # (tsreport if touse, ivreg2.ado:410-418). The kernel block above already
  # warns via time_index; cover the operators-without-kernel case here,
  # on the estimation sample (Stata's touse — post lag-induced shrinkage).
  if (isTRUE(parsed$has_ts_ops) && is.null(kernel)) {
    mf_rows_gap <- match(rownames(parsed$model_frame), rownames(data))
    ivar_vec_gap <- if (!is.null(ivar)) data[[ivar]][mf_rows_gap] else NULL
    ti_gap <- .build_time_index(data[[tvar]][mf_rows_gap], ivar_vec_gap)
    if (ti_gap$n_gaps > 0L) {
      warning("Time variable '", tvar, "' has ", ti_gap$n_gaps,
              " gap(s) in relevant range.", call. = FALSE)
    }
  }

  # --- Extract ivar_vec for SW (when kernel is not used, ivar_vec not yet built) ---
  sw_ivar_vec <- NULL
  if (isTRUE(opts$sw) && !is.null(ivar)) {
    mf_rows_sw <- match(rownames(parsed$model_frame), rownames(data))
    if (!ivar %in% names(data)) {
      stop("Panel variable '", ivar, "' not found in data.", call. = FALSE)
    }
    sw_ivar_vec <- data[[ivar]][mf_rows_sw]
    if (anyNA(sw_ivar_vec)) {
      stop("Panel variable '", ivar, "' contains NA values.", call. = FALSE)
    }
  }

  list(
    parsed = parsed, sdofminus = sdofminus, b0 = b0,
    wmatrix = wmatrix, smatrix = smatrix,
    cluster_vec = cluster_vec, cluster_var_name = cluster_var_name,
    M = M, M1 = M1, M2 = M2,
    time_index = time_index, unsort_order = unsort_order,
    partial_ct = partial_ct, partial_names = partial_names,
    partialcons = partialcons, n_physical = n_physical, w_raw = w_raw,
    bw = bw, ivar_vec = sw_ivar_vec
  )
}
