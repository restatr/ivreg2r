# --------------------------------------------------------------------------
# ivreg2 pipeline: argument validation and normalization
# Called only by ivreg2(); see R/ivreg2.R for the orchestrator.
# --------------------------------------------------------------------------

#' Validate and normalize all user arguments to ivreg2().
#'
#' Performs all type/value checks, option routing (kiefer→kernel+vcov,
#' dkraay→bw+kernel, pweight→robust/HAC, kernel→HAC/AC), and
#' method/kclass/fuller mutual-exclusion checks.
#'
#' @return Named list of normalized options.
#' @noRd
.validate_and_normalize_args <- function(vcov, clusters, endog, orthog,
                                          redundant, method, kclass, fuller,
                                          coviv, small, dofminus, sdofminus,
                                          weight_type, kernel, bw, tvar, ivar,
                                          kiefer, dkraay, wmatrix, smatrix,
                                          b0, partial, nopartialsmall,
                                          center, psd, reduced_form,
                                          first_stage, noid, sw,
                                          model_flag, x_flag, y_flag) {

  # --- Validate vcov ---
  if (!is.character(vcov) || length(vcov) != 1L || is.na(vcov)) {
    stop("`vcov` must be a single character string.", call. = FALSE)
  }
  # Redirects for known-wrong values, matched case-insensitively so that
  # "hc1" and "Cluster" get the same guidance as "HC1" and "cluster".
  if (toupper(vcov) %in% c("HC0", "HC1")) {
    stop('vcov = "', vcov, '" is not supported. ',
         'Use vcov = "robust" instead. ',
         'The `small` argument controls the finite-sample correction: ',
         'vcov = "robust" matches Stata `, robust`; ',
         'vcov = "robust", small = TRUE matches Stata `, robust small`.',
         call. = FALSE)
  }
  if (tolower(vcov) == "cluster") {
    stop('vcov = "', vcov, '" is not a VCE type: request clustering with the ',
         '`clusters` argument (e.g. clusters = ~ firmid), which combines with ',
         'any `vcov`. With `clusters` supplied, the default vcov = "iid" ',
         'already gives cluster-robust standard errors.', call. = FALSE)
  }
  # Match case-insensitively (so "hac", "Hac", "IID", "Robust" are all
  # accepted), then normalize to the canonical internal spelling used
  # throughout the rest of the pipeline and in printed output. An
  # unrecognized value is left untouched so the error below echoes exactly
  # what the user typed.
  valid_vcov <- c("iid", "robust", "HAC", "AC")
  vcov_idx <- match(tolower(vcov), tolower(valid_vcov))
  if (!is.na(vcov_idx)) {
    vcov <- valid_vcov[vcov_idx]
  }
  if (!vcov %in% valid_vcov) {
    stop('vcov = "', vcov, '" is not supported. ',
         'Supported values: ', paste0('"', valid_vcov, '"', collapse = ", "),
         '.', call. = FALSE)
  }
  if (!is.null(clusters) && !inherits(clusters, "formula")) {
    stop("`clusters` must be a one-sided formula (e.g. ~firmid).", call. = FALSE)
  }
  if (!is.logical(small) || length(small) != 1L || is.na(small)) {
    stop("`small` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(endog) && !is.character(endog)) {
    stop("`endog` must be a character vector or NULL.", call. = FALSE)
  }
  if (!is.null(endog) && anyDuplicated(endog)) {
    endog <- unique(endog)
    warning("`endog` contained duplicate entries; duplicates removed.",
            call. = FALSE)
  }
  if (!is.null(orthog) && !is.character(orthog)) {
    stop("`orthog` must be a character vector or NULL.", call. = FALSE)
  }
  if (!is.null(orthog) && anyDuplicated(orthog)) {
    orthog <- unique(orthog)
    warning("`orthog` contained duplicate entries; duplicates removed.",
            call. = FALSE)
  }
  if (!is.null(redundant) && !is.character(redundant)) {
    stop("`redundant` must be a character vector or NULL.", call. = FALSE)
  }
  if (!is.null(redundant) && anyDuplicated(redundant)) {
    redundant <- unique(redundant)
    warning("`redundant` contained duplicate entries; duplicates removed.",
            call. = FALSE)
  }
  if (!is.numeric(dofminus) || length(dofminus) != 1L || is.na(dofminus) ||
      !is.finite(dofminus) || dofminus < 0 || dofminus != trunc(dofminus) ||
      dofminus > .Machine$integer.max) {
    stop("`dofminus` must be a non-negative integer.", call. = FALSE)
  }
  if (!is.numeric(sdofminus) || length(sdofminus) != 1L || is.na(sdofminus) ||
      !is.finite(sdofminus) || sdofminus < 0 || sdofminus != trunc(sdofminus) ||
      sdofminus > .Machine$integer.max) {
    stop("`sdofminus` must be a non-negative integer.", call. = FALSE)
  }
  dofminus <- as.integer(dofminus)
  sdofminus <- as.integer(sdofminus)

  valid_wt <- c("aweight", "fweight", "pweight", "iweight")
  if (!is.character(weight_type) || length(weight_type) != 1L ||
      !weight_type %in% valid_wt) {
    stop('`weight_type` must be one of "aweight", "fweight", "pweight", or "iweight".',
         call. = FALSE)
  }

  if (!is.logical(coviv) || length(coviv) != 1L || is.na(coviv)) {
    stop("`coviv` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("`center` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(noid) || length(noid) != 1L || is.na(noid)) {
    stop("`noid` must be TRUE or FALSE.", call. = FALSE)
  }

  # --- Validate psd ---
  if (!is.null(psd)) {
    psd <- match.arg(psd, c("psd0", "psda"))
  }

  # --- Validate partial / nopartialsmall (type checks only) ---
  if (!is.null(partial) && !is.character(partial)) {
    stop("`partial` must be a character vector or NULL.", call. = FALSE)
  }
  if (!is.logical(nopartialsmall) || length(nopartialsmall) != 1L ||
      is.na(nopartialsmall)) {
    stop("`nopartialsmall` must be TRUE or FALSE.", call. = FALSE)
  }

  valid_rf <- c("none", "rf", "system")
  if (!is.character(reduced_form) || length(reduced_form) != 1L ||
      !reduced_form %in% valid_rf) {
    stop('`reduced_form` must be one of "none", "rf", or "system".',
         call. = FALSE)
  }

  if (!is.logical(first_stage) || length(first_stage) != 1L ||
      is.na(first_stage)) {
    stop("`first_stage` must be TRUE or FALSE.", call. = FALSE)
  }

  # --- Normalize kernel name early (before kiefer/dkraay checks) ---
  if (!is.null(kernel)) {
    kernel <- .validate_kernel(kernel)
  }

  # --- Validate kiefer ---
  if (!is.logical(kiefer) || length(kiefer) != 1L || is.na(kiefer)) {
    stop("`kiefer` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(kiefer)) {
    if (!is.null(dkraay)) {
      stop("incompatible options: `kiefer` and `dkraay`.", call. = FALSE)
    }
    if (is.null(tvar) || is.null(ivar)) {
      stop("kiefer requires panel data (both `tvar` and `ivar`).",
           call. = FALSE)
    }
    if (vcov == "robust" || !is.null(clusters)) {
      stop("kiefer is incompatible with robust VCE or clustering.",
           call. = FALSE)
    }
    if (!is.null(kernel) && kernel != "Truncated") {
      stop("kiefer is incompatible with kernel='", kernel, "'.",
           call. = FALSE)
    }
    if (!is.null(bw)) {
      stop("kiefer is incompatible with explicit `bw`.", call. = FALSE)
    }
    # Set kernel and vcov; bw deferred until after time_index is built
    kernel <- "Truncated"
    vcov <- "AC"
  }

  # --- Validate dkraay ---
  if (!is.null(dkraay)) {
    if (!is.numeric(dkraay) || length(dkraay) != 1L || !is.finite(dkraay) ||
        dkraay <= 0) {
      stop("`dkraay` must be a positive numeric scalar.", call. = FALSE)
    }
    if (is.null(tvar) || is.null(ivar)) {
      stop("dkraay requires panel data (both `tvar` and `ivar`).",
           call. = FALSE)
    }
    if (!is.null(bw)) {
      stop("cannot specify both `dkraay` and `bw`.", call. = FALSE)
    }
    if (!is.null(clusters)) {
      stop("cannot specify both `dkraay` and explicit `clusters`.",
           call. = FALSE)
    }
    bw <- dkraay
    if (is.null(kernel)) kernel <- "Bartlett"
    # DK clustering on tvar is handled after cluster construction
  }

  # pweight forces robust VCE (Stata: [pw=weight] → robust, unconditionally;
  # ivreg2.ado lines 353-357 set the robust flag before any VCE routing, so
  # [pw=] + bw() yields HAC and [pw=] + kiefer yields a heteroskedasticity-
  # robust Kiefer VCE). Must run BEFORE kernel routing so kernel + pweight +
  # iid → HAC (not AC), and must also promote an explicit (or kiefer-implied)
  # "AC" to "HAC".
  if (weight_type == "pweight" && is.null(clusters)) {
    if (vcov == "iid") {
      .warn_vce_promotion("pweight", "iid", "robust")
      vcov <- "robust"
    } else if (vcov == "AC") {
      .warn_vce_promotion("pweight", "AC", "HAC")
      vcov <- "HAC"
    }
  }

  # --- Validate bw / tvar / ivar ---
  if (!is.null(bw) && is.null(kernel)) {
    # bw specified without kernel: default to Bartlett
    kernel <- "Bartlett"
  }
  if (!is.null(kernel)) {
    if (is.null(bw) && !isTRUE(kiefer)) {
      stop("bandwidth `bw` is required when `kernel` is specified.",
           call. = FALSE)
    }
    if (!is.null(bw)) .validate_bandwidth(bw, kernel)

    # VCE inference: kernel + iid → AC; kernel + robust → HAC
    if (vcov == "iid") {
      vcov <- "AC"
    } else if (vcov == "robust") {
      vcov <- "HAC"
    }
  }
  # HAC/AC without kernel is an error
  if (vcov %in% c("HAC", "AC") && is.null(kernel)) {
    stop('vcov = "', vcov, '" requires `kernel` to be specified.',
         call. = FALSE)
  }
  # tvar required for HAC/AC
  if (vcov %in% c("HAC", "AC") && is.null(tvar)) {
    stop('Time variable `tvar` is required for vcov = "', vcov, '".',
         call. = FALSE)
  }
  # tvar/ivar must be character
  if (!is.null(tvar) && (!is.character(tvar) || length(tvar) != 1L)) {
    stop("`tvar` must be a single character string.", call. = FALSE)
  }
  if (!is.null(ivar) && (!is.character(ivar) || length(ivar) != 1L)) {
    stop("`ivar` must be a single character string.", call. = FALSE)
  }

  # --- Validate wmatrix / smatrix (type checks) ---
  # Must come BEFORE method/kclass/fuller promotion so parameter-specific
  # error messages fire before method has been changed.
  if (!is.null(wmatrix)) {
    if (!is.matrix(wmatrix) || !is.numeric(wmatrix))
      stop("`wmatrix` must be a numeric matrix.", call. = FALSE)
    if (any(!is.finite(wmatrix)))
      stop("`wmatrix` must not contain NA, NaN, or infinite values.",
           call. = FALSE)
  }
  if (!is.null(smatrix)) {
    if (!is.matrix(smatrix) || !is.numeric(smatrix))
      stop("`smatrix` must be a numeric matrix.", call. = FALSE)
    if (any(!is.finite(smatrix)))
      stop("`smatrix` must not contain NA, NaN, or infinite values.",
           call. = FALSE)
  }
  # --- Normalize method early (needed for b0 checks below) ---
  if (!is.character(method) || length(method) != 1L) {
    stop("`method` must be a single character string.", call. = FALSE)
  }
  method <- tolower(method)

  # iweight restricts to IID VCE only (Stata: ivreg2.ado lines 343-347)
  # Placed after tolower(method) and bw/kernel routing so all promotions
  # have already happened.
  if (weight_type == "iweight") {
    if (!is.null(clusters))
      stop("iweights not allowed with cluster VCE.", call. = FALSE)
    if (!is.null(kernel) || !is.null(bw))
      stop("iweights not allowed with kernel-based VCE.", call. = FALSE)
    if (method %in% c("gmm2s", "cue"))
      stop("iweights not allowed with GMM estimation.", call. = FALSE)
    if (!is.null(wmatrix))
      stop("iweights not allowed with GMM estimation (`wmatrix`).", call. = FALSE)
    if (!is.null(smatrix))
      stop("iweights not allowed with GMM estimation (`smatrix`).", call. = FALSE)
    if (vcov != "iid")
      stop("iweights not allowed with robust or HAC VCE.", call. = FALSE)
  }

  # --- Validate sw ---
  if (!is.logical(sw) || length(sw) != 1L || is.na(sw)) {
    stop("`sw` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(sw)) {
    if (is.null(ivar)) {
      stop("Stock-Watson VCE requires panel data (`ivar`).", call. = FALSE)
    }
    if (!is.null(clusters)) {
      stop("Stock-Watson VCE not supported with clustering.", call. = FALSE)
    }
    if (!is.null(kernel) || !is.null(bw)) {
      stop("Stock-Watson VCE not supported with kernel-based VCE.", call. = FALSE)
    }
    if (isTRUE(kiefer)) {
      stop("Stock-Watson VCE not supported with Kiefer VCE.", call. = FALSE)
    }
    if (!is.null(dkraay)) {
      stop("Stock-Watson VCE not supported with Driscoll-Kraay VCE.", call. = FALSE)
    }
    if (weight_type == "fweight") {
      stop("fweights not supported with Stock-Watson VCE.", call. = FALSE)
    }
    if (weight_type == "iweight") {
      stop("iweights not supported with Stock-Watson VCE.", call. = FALSE)
    }
    # SW forces robust VCE (Stata ivreg2.ado line 3758)
    if (vcov == "iid") {
      .warn_vce_promotion("sw", "iid", "robust")
      vcov <- "robust"
    }
  }

  # --- Validate b0 (type checks) ---
  if (!is.null(b0)) {
    if (!is.numeric(b0) || !is.null(dim(b0)))
      stop("`b0` must be a numeric vector.", call. = FALSE)
    if (any(!is.finite(b0)))
      stop("`b0` must contain only finite values.", call. = FALSE)
    # b0 implies method = "cue" if method is default
    if (method == "2sls") method <- "cue"
    # b0 is incompatible with gmm2s, liml, kclass, wmatrix, fuller
    if (method == "gmm2s")
      stop("Cannot specify `b0` with `method = \"gmm2s\"`.", call. = FALSE)
    if (method == "liml")
      stop("Cannot specify `b0` with `method = \"liml\"`.", call. = FALSE)
    if (method %in% c("kclass"))
      stop("Cannot specify `b0` with `method = \"kclass\"`.", call. = FALSE)
    if (!is.null(kclass))
      stop("Cannot specify `b0` with `kclass`.", call. = FALSE)
    if (fuller > 0)
      stop("Cannot specify `b0` with `fuller`.", call. = FALSE)
    if (!is.null(wmatrix))
      stop("Cannot specify `b0` with `wmatrix`.", call. = FALSE)
  }

  # Mutual exclusion with kclass/fuller (check raw params, not promoted method)
  if (!is.null(wmatrix) && !is.null(kclass))
    stop("Cannot specify `wmatrix` with `kclass`.", call. = FALSE)
  if (!is.null(smatrix) && !is.null(kclass))
    stop("Cannot specify `smatrix` with `kclass`.", call. = FALSE)
  if (!is.null(wmatrix) && fuller > 0)
    stop("Cannot specify `wmatrix` with `fuller`.", call. = FALSE)
  if (!is.null(smatrix) && fuller > 0)
    stop("Cannot specify `smatrix` with `fuller`.", call. = FALSE)

  # --- Validate method / kclass / fuller ---
  valid_methods <- c("2sls", "liml", "kclass", "gmm2s", "cue")
  if (!method %in% valid_methods) {
    stop('`method` must be one of "2sls", "liml", "kclass", "gmm2s", or ',
         '"cue"; you supplied "', method, '".', call. = FALSE)
  }
  if (!is.numeric(fuller) || length(fuller) != 1L || !is.finite(fuller)) {
    stop("`fuller` must be a single finite numeric value.", call. = FALSE)
  }
  if (fuller < 0) {
    stop("`fuller` must be non-negative.", call. = FALSE)
  }
  if (!is.null(kclass)) {
    if (!is.numeric(kclass) || length(kclass) != 1L || !is.finite(kclass)) {
      stop("`kclass` must be a single finite numeric value.", call. = FALSE)
    }
    if (kclass < 0) {
      stop("`kclass` must be non-negative.", call. = FALSE)
    }
  }
  # Mutual exclusion: fuller and kclass cannot both be specified
  if (fuller > 0 && !is.null(kclass)) {
    stop("Cannot specify both `fuller` and `kclass`.", call. = FALSE)
  }
  # GMM2S is incompatible with fuller and kclass
  if (method == "gmm2s" && fuller > 0) {
    stop("Cannot specify `fuller` with `method = \"gmm2s\"`.", call. = FALSE)
  }
  if (method == "gmm2s" && !is.null(kclass)) {
    stop("Cannot specify `kclass` with `method = \"gmm2s\"`.", call. = FALSE)
  }
  # CUE is incompatible with fuller, kclass, wmatrix, smatrix
  if (method == "cue" && fuller > 0) {
    stop("Cannot specify `fuller` with `method = \"cue\"`.", call. = FALSE)
  }
  if (method == "cue" && !is.null(kclass)) {
    stop("Cannot specify `kclass` with `method = \"cue\"`.", call. = FALSE)
  }
  if (method == "cue" && !is.null(wmatrix)) {
    stop("Cannot specify `wmatrix` with `method = \"cue\"`.", call. = FALSE)
  }
  if (method == "cue" && !is.null(smatrix)) {
    stop("Cannot specify `smatrix` with `method = \"cue\"`.", call. = FALSE)
  }
  if (method == "cue" && !is.null(partial) && is.null(b0)) {
    warning("FWL invariance does not hold for CUE: results with `partial` may ",
            "differ from the equivalent unpartialled specification ",
            "(Baum, Schaffer & Stillman, 2007, p. 484).", call. = FALSE)
  }
  # wmatrix/smatrix mutual exclusion with liml (check after method validation)
  if (!is.null(wmatrix) && method %in% c("liml", "kclass"))
    stop("Cannot specify `wmatrix` with method = \"", method, "\".",
         call. = FALSE)
  if (!is.null(smatrix) && method %in% c("liml", "kclass"))
    stop("Cannot specify `smatrix` with method = \"", method, "\".",
         call. = FALSE)
  # fuller implies liml (but not for gmm2s — already guarded above)
  if (fuller > 0 && method != "liml" && method != "gmm2s") {
    method <- "liml"
  }
  # kclass supplied implies method = "kclass" (but not for gmm2s — already guarded)
  if (!is.null(kclass) && method != "kclass" && method != "gmm2s") {
    if (method == "liml") {
      stop('Cannot specify `kclass` with `method = "liml"`.', call. = FALSE)
    }
    method <- "kclass"
  }
  # method = "kclass" requires an explicit kclass value
  if (method == "kclass" && is.null(kclass)) {
    stop('`method = "kclass"` requires a numeric `kclass` value.',
         call. = FALSE)
  }
  # coviv is only meaningful for LIML/kclass — warn and reset otherwise
  # (this also covers gmm2s and cue: neither is in the liml/kclass set, so
  # a separate gmm2s-specific check is unreachable and has been removed)
  if (coviv && !method %in% c("liml", "kclass")) {
    warning('`coviv` is ignored: it applies only to LIML/k-class estimation ',
            '(method = "', method, '").', call. = FALSE)
    coviv <- FALSE
  }

  # fweight + kernel incompatible (Stata ivreg2.ado:335)
  if (weight_type == "fweight" && !is.null(kernel)) {
    stop("fweights not allowed with kernel-based VCE.", call. = FALSE)
  }

  # Canonicalize ts-operator entries in option varlists so they match the
  # canonical term labels produced by the formula pass — "l(unem,1)" matches
  # "l(unem, 1)", and ranges like "l(unem, 1:2)" expand to component terms
  # (Stata analogue: fvexpand normalization before option-varlist matching,
  # ivreg2.ado:511-530). No-ops for entries without ts operators.
  endog <- unique(.canonicalize_ts_labels(endog))
  orthog <- unique(.canonicalize_ts_labels(orthog))
  redundant <- unique(.canonicalize_ts_labels(redundant))
  partial <- unique(.canonicalize_ts_labels(partial))

  list(
    method = method, vcov = vcov, kernel = kernel, bw = bw,
    kiefer = kiefer, dkraay = dkraay, tvar = tvar, ivar = ivar,
    coviv = coviv, small = small, dofminus = dofminus, sdofminus = sdofminus,
    endog = endog, orthog = orthog, redundant = redundant,
    weight_type = weight_type, center = center, psd = psd,
    partial = partial, nopartialsmall = nopartialsmall,
    reduced_form = reduced_form,
    first_stage_flag = first_stage,
    kclass = kclass, fuller = fuller, b0 = b0,
    wmatrix = wmatrix, smatrix = smatrix,
    noid = noid, sw = sw,
    model_flag = model_flag, x_flag = x_flag, y_flag = y_flag
  )
}

#' Warn that an option forces a more general VCE, overriding the requested
#' `vcov`.
#'
#' Shared by the pweight -> robust/HAC and sw -> robust promotions: both
#' silently substitute a more general `vcov` than the user asked for, and
#' planning/31 R4 rules that this substitution gets a fit-time `warning()`
#' (not `message()`), so it is independently testable with
#' `expect_warning()`.
#' @param option Character: the option name forcing the promotion (e.g.
#'   `"pweight"`, `"sw"`).
#' @param from Character: the `vcov` value being overridden.
#' @param to Character: the `vcov` value substituted for it.
#' @noRd
.warn_vce_promotion <- function(option, from, to) {
  warning(option, ' implies robust VCE; overriding vcov = "', from,
          '" to vcov = "', to, '".', call. = FALSE)
}
