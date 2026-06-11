# --------------------------------------------------------------------------
# ivreg2 — pipeline helpers defined first, then the main function
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Pipeline helpers (internal, called only by ivreg2)
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
  if (!is.character(vcov) || length(vcov) != 1L) {
    stop("`vcov` must be a single character string.", call. = FALSE)
  }
  if (vcov %in% c("HC0", "HC1")) {
    stop('vcov = "', vcov, '" is no longer supported. ',
         'Use vcov = "robust" instead. ',
         'The `small` argument controls the finite-sample correction: ',
         'vcov = "robust" matches Stata `, robust`; ',
         'vcov = "robust", small = TRUE matches Stata `, robust small`.',
         call. = FALSE)
  }
  valid_vcov <- c("iid", "robust", "HAC", "AC")
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
      message('pweight implies robust VCE; overriding vcov = "iid" to vcov = "robust".')
      vcov <- "robust"
    } else if (vcov == "AC") {
      message('pweight implies robust VCE; overriding vcov = "AC" to vcov = "HAC".')
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
  }
  if (!is.null(smatrix)) {
    if (!is.matrix(smatrix) || !is.numeric(smatrix))
      stop("`smatrix` must be a numeric matrix.", call. = FALSE)
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
    stop('`method` must be one of "2sls", "liml", "kclass", "gmm2s", or "cue".',
         call. = FALSE)
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
  # coviv is only meaningful for LIML/kclass — silently ignored otherwise
  if (coviv && !method %in% c("liml", "kclass")) {
    coviv <- FALSE
  }
  # GMM2S is incompatible with coviv
  if (coviv && method == "gmm2s") {
    coviv <- FALSE
  }

  # fweight + kernel incompatible (Stata ivreg2.ado:335)
  if (weight_type == "fweight" && !is.null(kernel)) {
    stop("fweights not allowed with kernel-based VCE.", call. = FALSE)
  }

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
    stop("`dofminus` + `sdofminus` too large: N - K - dofminus - sdofminus = ",
         parsed$N - parsed$K - dofminus - sdofminus,
         " (must be > 0).", call. = FALSE)
  }
  if (parsed$is_iv && parsed$N - parsed$L - dofminus - sdofminus <= 0L) {
    stop("`dofminus` + `sdofminus` too large: N - L - dofminus - sdofminus = ",
         parsed$N - parsed$L - dofminus - sdofminus,
         " (must be > 0).", call. = FALSE)
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
      stop("`dofminus` + `sdofminus` too large: N - K - dofminus - sdofminus = ",
           N_check - parsed$K - dofminus - sdofminus,
           " (must be > 0).", call. = FALSE)
    }
    if (parsed$is_iv && N_check - parsed$L - dofminus - sdofminus <= 0L) {
      stop("`dofminus` + `sdofminus` too large: N - L - dofminus - sdofminus = ",
           N_check - parsed$L - dofminus - sdofminus,
           " (must be > 0).", call. = FALSE)
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

    # Re-validate dimensions after partialling
    if (parsed$N - parsed$K - dofminus - sdofminus <= 0L) {
      stop("After partialling: N - K - dofminus - sdofminus = ",
           parsed$N - parsed$K - dofminus - sdofminus,
           " (must be > 0).", call. = FALSE)
    }
    if (parsed$is_iv && parsed$N - parsed$L - dofminus - sdofminus <= 0L) {
      stop("After partialling: N - L - dofminus - sdofminus = ",
           parsed$N - parsed$L - dofminus - sdofminus,
           " (must be > 0).", call. = FALSE)
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
        stop("Cluster variable '", cluster_var_name, "' contains NA values.",
             call. = FALSE)
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
        stop("Cluster variable '", cluster_var_name[1L],
             "' contains NA values.", call. = FALSE)
      if (anyNA(cv2))
        stop("Cluster variable '", cluster_var_name[2L],
             "' contains NA values.", call. = FALSE)
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
        stop("invalid bandwidth in option bw() - cannot exceed timespan of data",
             call. = FALSE)
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


#' Determine estimation path, build omega_fn, dispatch to fitter.
#'
#' Handles wmatrix/smatrix routing, auto-bandwidth (before estimation for GMM,
#' after for non-GMM), omega_fn closure construction, and estimation dispatch
#' (7 branches: OLS, 2SLS, LIML/kclass, GMM2S, CUE, gmmw, smatrix-gmm2s).
#'
#' @return Named list with fit, method (possibly promoted), bw (resolved),
#'   wmatrix (possibly nulled), omega_fn.
#' @noRd
.resolve_and_fit <- function(parsed, prep, opts) {
  method <- opts$method
  vcov <- opts$vcov
  kernel <- opts$kernel
  small <- opts$small
  dofminus <- opts$dofminus
  fuller <- opts$fuller
  kclass <- opts$kclass
  weight_type <- opts$weight_type
  center <- opts$center
  psd <- opts$psd
  ivar <- opts$ivar
  # name-normalized in .prepare_model (.match_user_matrix)
  wmatrix <- prep$wmatrix
  smatrix <- prep$smatrix

  sdofminus   <- prep$sdofminus
  b0          <- prep$b0
  cluster_vec <- prep$cluster_vec
  time_index  <- prep$time_index
  bw          <- prep$bw

  # --- Determine effective estimation path for wmatrix/smatrix ---
  use_wmatrix <- !is.null(wmatrix)
  use_smatrix <- !is.null(smatrix)

  # wmatrix without robust/cluster/kernel/gmm2s → warn & ignore (Stata behavior)
  if (use_wmatrix && is.null(cluster_vec) && vcov == "iid" &&
      is.null(kernel) && method != "gmm2s") {
    warning("`wmatrix` is ignored without robust VCE, clustering, or HAC; ",
            "using standard IV estimation.", call. = FALSE)
    use_wmatrix <- FALSE
    wmatrix <- NULL
  }

  # smatrix alone implies efficient GMM if method is default.
  if (use_smatrix && !use_wmatrix && method == "2sls") {
    method <- "gmm2s"
  }

  # --- Build omega_fn closure ---
  omega_fn <- NULL
  needs_omega_fn <- method %in% c("gmm2s", "cue") || (use_wmatrix && !use_smatrix)

  if (needs_omega_fn) {
    # Resolve auto-bw using step-1 residuals (Stata ivreg2.ado:970-980)
    if (is.character(bw) && tolower(bw) == "auto") {
      if (!is.null(ivar)) {
        stop("Automatic bandwidth selection not available for panel data.",
             call. = FALSE)
      }
      if (use_wmatrix) {
        abw_resid <- .wmatrix_first_step_resid(parsed, wmatrix)
      } else {
        fit_step1 <- .fit_2sls(parsed, small = FALSE, dofminus = dofminus,
                               sdofminus = sdofminus)
        abw_resid <- fit_step1$residuals
      }
      bw <- .auto_bandwidth(
        resid = abw_resid,
        Z = parsed$Z,
        time_index = time_index,
        kernel = kernel,
        has_intercept = parsed$has_intercept,
        N = parsed$N
      )
      max_bw <- (time_index$T_span - 1) / time_index$tdelta
      if (bw > max_bw) {
        warning("Automatic bandwidth (", formatC(bw, format = "g"),
                ") exceeds time-span limit (", formatC(max_bw, format = "g"),
                "); capping at ", formatC(max_bw, format = "g"), ".",
                call. = FALSE)
        bw <- max_bw
      }
    }

    # Capture VCE parameters in closure
    gmm_Z <- parsed$Z
    gmm_w <- parsed$weights
    gmm_N <- parsed$N
    gmm_cv <- cluster_vec
    gmm_dm <- dofminus
    gmm_wt <- weight_type
    gmm_k  <- kernel
    gmm_bw <- bw
    gmm_ti <- time_index
    gmm_is_iid <- is.null(cluster_vec) && vcov == "iid" && is.null(kernel)
    gmm_center <- center
    gmm_psd <- psd
    gmm_vcov_type <- vcov
    gmm_sw <- isTRUE(opts$sw)
    gmm_ivar_vec <- prep$ivar_vec
    if (vcov == "AC") {
      if (!is.null(parsed$weights)) {
        gmm_ZwZ <- crossprod(parsed$Z, parsed$weights * parsed$Z)
      } else {
        gmm_ZwZ <- crossprod(parsed$Z)
      }
    } else {
      gmm_ZwZ <- NULL
    }
    omega_fn <- function(resid) {
      .compute_moment_cov(gmm_Z, resid, gmm_w, gmm_cv, gmm_N,
                          iid = gmm_is_iid,
                          dofminus = gmm_dm, weight_type = gmm_wt,
                          kernel = gmm_k, bw = gmm_bw, time_index = gmm_ti,
                          center = gmm_center, psd = gmm_psd,
                          vcov_type = gmm_vcov_type, ZwZ = gmm_ZwZ,
                          sw = gmm_sw, ivar_vec = gmm_ivar_vec)
    }
  }

  # --- Estimation dispatch ---
  if (use_smatrix && !use_wmatrix) {
    omega_fn <- function(resid) smatrix
    fit <- .fit_gmm2s(parsed, small = small, dofminus = dofminus,
                      sdofminus = sdofminus, omega_fn = omega_fn)

  } else if (use_wmatrix && method == "gmm2s") {
    wstep_resid <- .wmatrix_first_step_resid(parsed, wmatrix)
    if (use_smatrix) {
      omega_fn <- function(resid) smatrix
    } else {
      omega_fn_base <- omega_fn
      omega_fn <- function(resid) omega_fn_base(wstep_resid)
    }
    fit <- .fit_gmm2s(parsed, small = small, dofminus = dofminus,
                      sdofminus = sdofminus, omega_fn = omega_fn)

  } else if (use_wmatrix) {
    if (use_smatrix) {
      omega_fn <- function(resid) smatrix
    }
    fit <- .fit_gmm_wmatrix(parsed, small = small, dofminus = dofminus,
                             sdofminus = sdofminus,
                             W = wmatrix, omega_fn = omega_fn)
    method <- "gmmw"

  } else if (method == "cue") {
    fit <- .fit_cue(parsed, small = small, dofminus = dofminus,
                    sdofminus = sdofminus, omega_fn = omega_fn, b0 = b0,
                    iid = gmm_is_iid)

  } else if (method == "gmm2s") {
    fit <- .fit_gmm2s(parsed, small = small, dofminus = dofminus,
                      sdofminus = sdofminus, omega_fn = omega_fn)

  } else if (method %in% c("liml", "kclass")) {
    fit <- .fit_kclass(parsed, method = method, kclass = kclass, fuller = fuller,
                small = small, dofminus = dofminus, sdofminus = sdofminus)
  } else if (parsed$is_iv && parsed$K1 > 0L) {
    fit <- .fit_2sls(parsed, small = small, dofminus = dofminus,
              sdofminus = sdofminus)
  } else {
    # 1-part OLS, or empty-endogenous (K1 = 0): estimation is exact OLS —
    # the excluded instruments enter only the diagnostics (Stata posts
    # e(model) = "ols" for `(=z)` models, ivreg2.ado:2101-2106).
    fit <- .fit_ols(parsed, small = small, dofminus = dofminus,
             sdofminus = sdofminus)
  }

  # --- Resolve bw = "auto" (Newey-West 1994) ---
  # For non-GMM: resolve AFTER estimation (need residuals).
  # For GMM/CUE: already resolved above.
  if (!method %in% c("gmm2s", "gmmw", "cue") && is.character(bw) && tolower(bw) == "auto") {
    if (!is.null(ivar)) {
      stop("Automatic bandwidth selection not available for panel data.",
           call. = FALSE)
    }
    abw_Z <- if (is.null(parsed$Z)) parsed$X else parsed$Z
    bw <- .auto_bandwidth(
      resid = fit$residuals,
      Z = abw_Z,
      time_index = time_index,
      kernel = kernel,
      has_intercept = parsed$has_intercept,
      N = parsed$N
    )
    max_bw <- (time_index$T_span - 1) / time_index$tdelta
    if (bw > max_bw) {
      warning("Automatic bandwidth (", formatC(bw, format = "g"),
              ") exceeds time-span limit (", formatC(max_bw, format = "g"),
              "); capping at ", formatC(max_bw, format = "g"), ".",
              call. = FALSE)
      bw <- max_bw
    }
  }

  list(fit = fit, method = method, bw = bw, wmatrix = wmatrix,
       omega_fn = omega_fn)
}


#' Invert a symmetric positive-definite matrix with graceful failure.
#'
#' Tries Cholesky (chol2inv), falls back to qr.solve, and returns NULL with
#' a warning if the matrix is computationally singular. Used for the stored
#' weighting matrix W, where Stata's hard error (exit 506) is relaxed to the
#' package's warn-and-skip convention.
#'
#' @return Inverse matrix, or NULL if singular.
#' @noRd
.safe_inverse <- function(A, what = "matrix") {
  inv <- tryCatch(chol2inv(chol(A)), error = function(e) NULL)
  if (is.null(inv)) {
    inv <- tryCatch(qr.solve(A), error = function(e) NULL)
  }
  if (is.null(inv)) {
    warning("Could not invert ", what, " (computationally singular); ",
            "storing NULL.", call. = FALSE)
  }
  inv
}


#' Compute the L x L moment covariance S from the final residuals.
#'
#' Thin wrapper that assembles the [.compute_moment_cov()] arguments from the
#' fit/parsed/prep/opts objects (including the AC path's Z'WZ precompute).
#' Shared by the stored `fit$S` and the psd-corrected VCV assembly in
#' [.compute_vcov()], so the two are the same matrix by construction.
#' psd correction happens inside [.compute_moment_cov()] when `opts$psd`
#' is set.
#'
#' Must be called BEFORE the unsort step: HAC/AC paths index residuals by the
#' sorted `prep$time_index`.
#'
#' @return L x L symmetric matrix S (no dimnames).
#' @noRd
.fresh_moment_cov <- function(fit, parsed, prep, opts, bw, center) {
  Zmat <- if (parsed$is_iv) parsed$Z else parsed$X
  is_iid <- is.null(prep$cluster_vec) && opts$vcov == "iid" &&
    is.null(opts$kernel)
  ZwZ <- NULL
  if (opts$vcov == "AC") {
    ZwZ <- if (is.null(parsed$weights)) {
      crossprod(Zmat)
    } else {
      crossprod(Zmat, parsed$weights * Zmat)
    }
  }
  .compute_moment_cov(
    Z = Zmat, residuals = fit$residuals, weights = parsed$weights,
    cluster_vec = prep$cluster_vec, N = parsed$N, iid = is_iid,
    dofminus = opts$dofminus, weight_type = opts$weight_type,
    kernel = opts$kernel, bw = bw, time_index = prep$time_index,
    center = center, psd = opts$psd, vcov_type = opts$vcov,
    ZwZ = ZwZ, sw = isTRUE(opts$sw), ivar_vec = prep$ivar_vec
  )
}


#' Compute the stored moment-covariance S and weighting matrix W.
#'
#' Builds the matrices posted on the fitted object as `fit$S` and `fit$W`,
#' matching Stata's `e(S)` and `e(W)` (ivreg2.ado:1905-1918):
#'
#' - S: user `smatrix` is echoed (never recomputed); GMM methods reuse the
#'   defining Omega already computed by the fitter (`fit$omega`); all other
#'   paths (OLS, 2SLS, LIML/k-class, any VCE) compute S fresh from the final
#'   residuals via [.fresh_moment_cov()].
#' - W: NULL for LIML/k-class ("No weighting matrix defined", ivreg2.ado:1912);
#'   user `wmatrix` is echoed (even under gmm2s, where it is the first-step
#'   W); gmm2s/CUE/b0/smatrix-implied-GMM use `solve(S)`; 1-step OLS/2SLS use
#'   `(1/sigma^2) (Z'WZ/N)^{-1}` from the final residuals with the
#'   large-sample `sigma^2 = RSS/(N - dofminus)` (s_gmm1s, ivreg2.ado:5287-5371;
#'   independent of `small`).
#'
#' Must be called BEFORE the unsort step: HAC/AC paths index residuals by the
#' sorted `prep$time_index`.
#'
#' @return list(S = matrix, W = matrix or NULL), dimnames = instrument names.
#' @noRd
.compute_stored_sw <- function(fit, parsed, prep, opts, est, bw, center) {
  method <- est$method
  Zmat <- if (parsed$is_iv) parsed$Z else parsed$X
  inames <- colnames(Zmat)

  # --- S ---
  if (!is.null(prep$smatrix)) {
    S <- prep$smatrix
  } else if (method %in% c("gmm2s", "gmmw", "cue")) {
    S <- fit$omega
  } else if (!is.null(fit$psd_moment_cov)) {
    # Reuse the corrected Omega the psd VCV assembly already computed
    # (.compute_vcov): same residuals, same correction — recomputing would
    # only duplicate the psd warning.
    S <- fit$psd_moment_cov
  } else {
    S <- .fresh_moment_cov(fit, parsed, prep, opts, bw, center)
  }
  dimnames(S) <- list(inames, inames)

  # --- W ---
  if (method %in% c("liml", "kclass")) {
    W <- NULL
  } else if (!is.null(est$wmatrix)) {
    W <- est$wmatrix
  } else if (method %in% c("gmm2s", "gmmw", "cue")) {
    W <- .safe_inverse(S, what = "moment covariance S for fit$W")
  } else {
    ZWZ <- if (is.null(parsed$weights)) {
      crossprod(Zmat)
    } else {
      crossprod(Zmat, parsed$weights * Zmat)
    }
    sigma2 <- fit$rss / (parsed$N - opts$dofminus)
    QZZinv <- .safe_inverse(ZWZ / parsed$N, what = "Z'Z/N for fit$W")
    W <- if (is.null(QZZinv)) NULL else QZZinv / sigma2
  }
  if (!is.null(W)) dimnames(W) <- list(inames, inames)

  list(S = S, W = W)
}


#' Apply VCV corrections and compute final variance-covariance matrix.
#'
#' GMM path: apply small-sample corrections to fit$vcov, recompute sigma.
#' Non-GMM path: select bread (k-class vs 2SLS for COVIV), dispatch to VCV
#' helpers. When `psd` is set, plain robust-family paths assemble the VCV
#' from the psd-corrected L x L moment covariance via [.vcov_from_omega()]
#' (Stata m_omega parity); the iid VCV is never psd-corrected.
#'
#' @return Updated fit with final vcov, sigma, df.residual.
#' @noRd
.compute_vcov <- function(fit, parsed, prep, opts, method, bw) {
  vcov <- opts$vcov
  kernel <- opts$kernel
  small <- opts$small
  dofminus <- opts$dofminus
  coviv <- opts$coviv
  weight_type <- opts$weight_type
  center <- opts$center
  psd <- opts$psd

  sdofminus   <- prep$sdofminus
  cluster_vec <- prep$cluster_vec
  time_index  <- prep$time_index
  M           <- prep$M

  if (method %in% c("gmm2s", "gmmw", "cue")) {
    needs_vcov_correction <- small
    if (needs_vcov_correction) {
      if (!is.null(cluster_vec)) {
        fit$vcov <- fit$vcov * ((parsed$N - 1) / (parsed$N - parsed$K - sdofminus)) *
          (M / (M - 1))
        fit$df.residual <- as.integer(M - 1L)
      } else {
        fit$vcov <- fit$vcov * ((parsed$N - dofminus) /
                                  (parsed$N - parsed$K - dofminus - sdofminus))
      }
    } else if (!is.null(cluster_vec)) {
      fit$df.residual <- as.integer(M - 1L)
    }
    if (small) {
      fit$sigma <- sqrt(fit$rss / (parsed$N - parsed$K - dofminus - sdofminus))
    }
    colnames(fit$vcov) <- rownames(fit$vcov) <- names(fit$coefficients)
  } else {

  # For LIML/kclass, select bread: k-class bread by default, 2SLS bread if coviv
  bread_vcov <- if (method %in% c("liml", "kclass") && !coviv) {
    fit$bread_kclass
  } else {
    fit$bread
  }

  # COVIV + IID: override the k-class IID VCV with 2SLS-bread IID VCV
  if (method %in% c("liml", "kclass") && coviv &&
      vcov == "iid" && is.null(cluster_vec)) {
    fit$vcov <- fit$sigma^2 * fit$bread
    colnames(fit$vcov) <- rownames(fit$vcov) <- names(fit$coefficients)
  }

  if (isTRUE(opts$sw) || !is.null(cluster_vec) ||
      vcov %in% c("robust", "HAC", "AC")) {
    # K1 = 0 (empty-endogenous) models are estimated by OLS: the meat basis
    # is X itself, like the 1-part OLS path (Stata's robust VCV for `(=z)`
    # models never involves Z).
    X_hat_vcov <- if (parsed$is_iv && parsed$K1 > 0L) fit$X_hat else parsed$X
    resid_vcov <- fit$residuals
  }

  if (!is.null(psd) &&
      (isTRUE(opts$sw) || !is.null(cluster_vec) ||
       vcov %in% c("robust", "HAC", "AC"))) {
    # F2 (D9b): Stata applies the psd0/psda correction to the L x L moment
    # covariance S inside m_omega (livreg2.do:607-617) and assembles the
    # plain non-GMM VCV from the corrected S by congruence with the
    # first-stage map A (s_iegmm, ivreg2.ado:5556). Correcting the K x K
    # meat or the final VCV diverges whenever the correction binds, so with
    # psd set we route every plain robust-family path (HC, cluster, two-way,
    # HAC, AC, DK, cluster+kernel, SW) through the corrected S. This S is
    # the same computation as the stored fit$S.
    omega_vcov <- .fresh_moment_cov(fit, parsed, prep, opts, bw, center)
    is_cluster_family <- !is.null(cluster_vec)
    # First-stage map A with X_hat = Z A. For K1 = 0 (OLS estimation with
    # surplus instruments), X is a subset of Z's columns, so A is the exact
    # 0/1 selection matrix — no solve needed.
    A_vcov <- if (!parsed$is_iv) {
      NULL
    } else if (parsed$K1 > 0L) {
      fit$proj_coef
    } else {
      A_sel <- matrix(0, nrow = parsed$L, ncol = parsed$K,
                      dimnames = list(colnames(parsed$Z), colnames(parsed$X)))
      A_sel[cbind(match(colnames(parsed$X), colnames(parsed$Z)),
                  seq_len(parsed$K))] <- 1
      A_sel
    }
    fit$vcov <- .vcov_from_omega(
      bread = bread_vcov,
      A = A_vcov,
      omega = omega_vcov,
      N = parsed$N, K = parsed$K, M = M, small = small,
      dofminus = dofminus, sdofminus = sdofminus,
      cluster = is_cluster_family
    )
    if (is_cluster_family) {
      fit$df.residual <- as.integer(M - 1L)
    }
    # Stash for .compute_stored_sw(): fit$S must be this exact matrix, and
    # recomputing it would emit a duplicate psd-correction warning.
    fit$psd_moment_cov <- omega_vcov
  } else if (isTRUE(opts$sw)) {
    fit$vcov <- .compute_sw_vcov(
      bread = bread_vcov, X_hat = X_hat_vcov, resid = resid_vcov,
      ivar_vec = prep$ivar_vec, N = parsed$N, K = parsed$K,
      small = small, dofminus = dofminus, sdofminus = sdofminus,
      weights = parsed$weights, weight_type = weight_type,
      center = center
    )
  } else if (!is.null(cluster_vec) && !is.null(kernel)) {
    is_twoway <- is.list(cluster_vec)
    fit$vcov <- .compute_cluster_kernel_vcov(
      bread = bread_vcov, X_hat = X_hat_vcov, resid = resid_vcov,
      cluster_vec = cluster_vec, time_index = time_index,
      kernel = kernel, bw = bw,
      N = parsed$N, K = parsed$K, M = M, small = small,
      dofminus = dofminus, sdofminus = sdofminus,
      weights = parsed$weights, weight_type = weight_type,
      is_twoway = is_twoway,
      center = center
    )
    fit$df.residual <- as.integer(M - 1L)
  } else if (vcov == "HAC") {
    fit$vcov <- .compute_hac_vcov(bread_vcov, X_hat_vcov, resid_vcov,
                                   time_index, kernel, bw,
                                   parsed$N, parsed$K,
                                   dofminus = dofminus, sdofminus = sdofminus,
                                   small = small,
                                   weights = parsed$weights,
                                   weight_type = weight_type,
                                   center = center)
  } else if (vcov == "AC") {
    fit$vcov <- .compute_ac_vcov(bread_vcov, X_hat_vcov, resid_vcov,
                                  time_index, kernel, bw,
                                  parsed$N, parsed$K,
                                  dofminus = dofminus, sdofminus = sdofminus,
                                  small = small,
                                  weights = parsed$weights,
                                  weight_type = weight_type)
  } else if (!is.null(cluster_vec)) {
    fit$vcov <- .compute_cl_vcov(bread_vcov, X_hat_vcov, resid_vcov,
                                  cluster_vec, parsed$N, parsed$K, M, small,
                                  dofminus = dofminus, sdofminus = sdofminus,
                                  weights = parsed$weights,
                                  weight_type = weight_type,
                                  center = center)
    fit$df.residual <- as.integer(M - 1L)
  } else if (vcov == "robust") {
    fit$vcov <- .compute_hc_vcov(bread_vcov, X_hat_vcov, resid_vcov,
                                  parsed$N, parsed$K,
                                  small = small, dofminus = dofminus,
                                  sdofminus = sdofminus,
                                  weights = parsed$weights,
                                  weight_type = weight_type,
                                  center = center)
  }

  }  # end of non-GMM VCV block

  # No final-VCV psd pass: the correction lives at the S level (the
  # psd branch above for plain paths; omega_fn for GMM paths), and the
  # iid VCV is never corrected, matching Stata (m_omega is not involved
  # in the iid VCV; psd is silently inert there).

  fit
}


#' Compute all diagnostic tests, reduced-form regression, and model F.
#'
#' Derives effective_vcov_type, resets center if no effect, runs 11 diagnostic
#' tests (gated by b0), reduced-form regression, and model F-test.
#'
#' @return Named list with diagnostics, first_stage, reduced_form_result,
#'   model_f_result, effective_vcov_type, center.
#' @noRd
.compute_diagnostics <- function(fit, parsed, prep, opts, method, bw) {
  vcov <- opts$vcov
  kernel <- opts$kernel
  small <- opts$small
  dofminus <- opts$dofminus
  fuller <- opts$fuller
  kiefer <- opts$kiefer
  endog <- opts$endog
  orthog <- opts$orthog
  redundant <- opts$redundant
  weight_type <- opts$weight_type
  center <- opts$center
  psd <- opts$psd
  reduced_form <- opts$reduced_form
  ivar <- opts$ivar
  b0 <- opts$b0
  noid <- opts$noid
  first_stage_flag <- opts$first_stage_flag

  sdofminus   <- prep$sdofminus
  cluster_vec <- prep$cluster_vec
  time_index  <- prep$time_index
  M           <- prep$M
  sw_flag     <- isTRUE(opts$sw)
  ivar_vec    <- prep$ivar_vec

  # Derive effective VCE type
  effective_vcov_type <- if (!is.null(cluster_vec)) {
    "CL"
  } else if (vcov == "HAC") {
    "HAC"
  } else if (vcov == "AC") {
    "AC"
  } else {
    vcov
  }

  # Warn and reset center if it has no effect
  if (center && effective_vcov_type %in% c("iid", "AC")) {
    warning("`center = TRUE` has no effect with ", effective_vcov_type,
            " VCE (centering only applies to robust/cluster/HAC).",
            call. = FALSE)
    center <- FALSE
  }

  diagnostics <- list()
  first_stage <- NULL
  if (parsed$is_iv) {

    # Overidentification test (D1)
    if (method %in% c("gmm2s", "gmmw", "cue")) {
      overid_test_name <- if (method %in% c("gmmw", "cue")) {
        "Hansen J"
      } else if (effective_vcov_type %in% c("iid", "AC")) {
        "Sargan"
      } else {
        "Hansen J"
      }
      diagnostics$overid <- list(
        stat = fit$j_stat, p = fit$j_p,
        df = fit$j_df, test_name = overid_test_name
      )
    } else {
    diagnostics$overid <- .compute_overid_test(
      Z = parsed$Z, X = parsed$X, y = parsed$y,
      residuals = fit$residuals, rss = fit$rss,
      weights = parsed$weights, cluster_vec = cluster_vec,
      vcov_type = effective_vcov_type, is_iv = parsed$is_iv,
      N = parsed$N, K = parsed$K, L = parsed$L,
      overid_df = parsed$overid_df, dofminus = dofminus,
      weight_type = weight_type,
      kernel = kernel, bw = bw, time_index = time_index,
      center = center, psd = psd,
      sw = sw_flag, ivar_vec = ivar_vec
    )
    }

    # --- Validate endog/orthog/redundant names unconditionally ---
    # (Stata validates at lines 480-531, before any b0/noid gating)
    endog_cols <- NULL
    if (!is.null(endog)) {
      bad <- setdiff(endog, parsed$endo_names)
      if (length(bad) > 0L) {
        stop("`endog` contains variables not in the endogenous list: ",
             paste0("'", bad, "'", collapse = ", "), ".", call. = FALSE)
      }
      endog_cols <- .expand_terms_to_colnames(
        endog, parsed$endo_names, parsed$endo_colnames, parsed$endo_assign
      )
    }

    orthog_cols <- NULL
    if (!is.null(orthog) && length(orthog) > 0L) {
      valid_terms <- c(parsed$excluded_names, parsed$exog_term_labels)
      bad <- setdiff(orthog, valid_terms)
      if (length(bad) > 0L) {
        stop("`orthog` contains variables not in the instrument list: ",
             paste0("'", bad, "'", collapse = ", "),
             ". Must be excluded or exogenous instruments (not endogenous ",
             "regressors or the intercept).", call. = FALSE)
      }
      orthog_in_excl <- intersect(orthog, parsed$excluded_names)
      orthog_in_exog <- intersect(orthog, parsed$exog_term_labels)
      orthog_cols <- c(
        .expand_terms_to_colnames(orthog_in_excl, parsed$excluded_names,
                                   parsed$excluded_colnames,
                                   parsed$excluded_assign),
        .expand_terms_to_colnames(orthog_in_exog, parsed$exog_term_labels,
                                   parsed$exog_colnames, parsed$exog_assign)
      )
    }

    if (!is.null(redundant) && length(redundant) > 0L) {
      bad <- setdiff(redundant, parsed$excluded_names)
      if (length(bad) > 0L) {
        stop("`redundant` contains variables not in the excluded instrument list: ",
             paste0("'", bad, "'", collapse = ", "), ".", call. = FALSE)
      }
    }

    # b0 suppresses all identification diagnostics (Stata line 3819)
    if (is.null(b0)) {

    # AR LIML overidentification (H3)
    if (method == "liml" && effective_vcov_type == "iid") {
      diagnostics$anderson_rubin_overid <- .compute_ar_liml_overid(
        lambda = fit$lambda, N = parsed$N,
        overid_df = parsed$overid_df, dofminus = dofminus
      )
    }

    # Identification tests (D2) + Stock-Yogo critical values (D3).
    # K1 = 0 (empty-endogenous): skipped entirely, matching Stata's master
    # gate `if endo1_ct > 0 & noid==""` (ivreg2.ado:1625) — idstat/widstat
    # are not posted for `(=z)` models.
    if (!noid && parsed$K1 > 0L) {
      # Kiefer fits: Stata's ranktest never receives the truncated kernel
      # (kiefer assigns bw after the parser builds bwopt/kernopt, so the
      # kernel options are empty), but it DOES receive the robust flag when
      # one is set -- which happens when pweights force robust. So plain
      # kiefer id tests are numerically iid, while pweight+kiefer id tests
      # are heteroskedasticity-robust (KP rk) with no kernel. Verified
      # against Stata 2026-06-10 (plain: Anderson/CD match rk output to
      # ~1e-11; pw+kiefer: KP rk LM/Wald differ from iid and match Stata).
      id_vcov_type <- if (isTRUE(kiefer)) {
        if (effective_vcov_type == "HAC") "robust" else "iid"
      } else {
        effective_vcov_type
      }
      id_kernel <- if (isTRUE(kiefer)) NULL else kernel
      id_tests <- .compute_id_tests(
        X = parsed$X, Z = parsed$Z, y = parsed$y,
        residuals = fit$residuals, weights = parsed$weights,
        cluster_vec = cluster_vec, vcov_type = id_vcov_type,
        N = parsed$N, K = parsed$K, L = parsed$L,
        K1 = parsed$K1, L1 = parsed$L1, M = M,
        endo_names = parsed$endo_colnames,
        excluded_names = parsed$excluded_colnames,
        has_intercept = parsed$has_intercept,
        dofminus = dofminus, sdofminus = sdofminus,
        weight_type = weight_type,
        kernel = id_kernel, bw = bw, time_index = time_index
      )
      diagnostics$underid        <- id_tests$underid
      diagnostics$weak_id        <- id_tests$weak_id
      diagnostics$weak_id_robust <- id_tests$weak_id_robust
      diagnostics$ccev           <- id_tests$ccev
      diagnostics$cdev           <- id_tests$cdev

      sy_method <- if (method %in% c("gmm2s", "gmmw")) {
        "2sls"
      } else if (method == "cue") {
        "liml"
      } else {
        method
      }
      diagnostics$weak_id_sy <- .stock_yogo_lookup(
        parsed$K1, parsed$L1, method = sy_method, fuller = fuller
      )
    }

    # First-stage diagnostics (E1), Anderson-Rubin (E3), Stock-Wright (J2),
    # and the default endogeneity test (E4) all concern the endogenous
    # regressors — skipped for K1 = 0, matching Stata (first-stage gated at
    # ivreg2.ado:1256/2071; AR and Stock-Wright flow through ranktest, which
    # is never called when endo1_ct == 0).
    if (parsed$K1 > 0L) {

    # First-stage diagnostics (E1)
    first_stage <- .compute_first_stage(
      X = parsed$X, Z = parsed$Z,
      weights = parsed$weights, cluster_vec = cluster_vec,
      vcov_type = effective_vcov_type,
      endo_names = parsed$endo_colnames,
      excluded_names = parsed$excluded_colnames,
      N = parsed$N, K = parsed$K, L = parsed$L,
      K1 = parsed$K1, L1 = parsed$L1, M = M,
      bread_2sls = fit$bread,
      dofminus = dofminus, sdofminus = sdofminus,
      weight_type = weight_type,
      kernel = kernel, bw = bw, time_index = time_index,
      center = center,
      sw = sw_flag, ivar_vec = ivar_vec,
      psd = psd
    )

    # Anderson-Rubin test (E3)
    diagnostics$anderson_rubin <- .compute_anderson_rubin(
      Z = parsed$Z, X = parsed$X, y = parsed$y,
      weights = parsed$weights, cluster_vec = cluster_vec,
      vcov_type = effective_vcov_type,
      N = parsed$N, K = parsed$K, L = parsed$L,
      K1 = parsed$K1, L1 = parsed$L1, M = M,
      endo_names = parsed$endo_colnames,
      excluded_names = parsed$excluded_colnames,
      dofminus = dofminus, sdofminus = sdofminus,
      weight_type = weight_type,
      kernel = kernel, bw = bw, time_index = time_index,
      center = center,
      sw = sw_flag, ivar_vec = ivar_vec,
      psd = psd
    )

    # Stock-Wright S statistic (J2)
    diagnostics$stock_wright <- .compute_stock_wright(
      Z = parsed$Z, X = parsed$X, y = parsed$y,
      weights = parsed$weights, cluster_vec = cluster_vec,
      vcov_type = effective_vcov_type,
      N = parsed$N, K1 = parsed$K1, L1 = parsed$L1,
      endo_names = parsed$endo_colnames, dofminus = dofminus,
      weight_type = weight_type,
      kernel = kernel, bw = bw, time_index = time_index,
      center = center, psd = psd,
      sw = sw_flag, ivar_vec = ivar_vec
    )

    # Endogeneity test / C-statistic (E4)
    diagnostics$endogeneity <- .compute_endogeneity_test(
      Z = parsed$Z, X = parsed$X, y = parsed$y,
      residuals = fit$residuals, rss = fit$rss,
      weights = parsed$weights, cluster_vec = cluster_vec,
      vcov_type = effective_vcov_type,
      N = parsed$N, K = parsed$K, L = parsed$L,
      K1 = parsed$K1, endo_names = parsed$endo_colnames,
      endog_vars = endog_cols, dofminus = dofminus,
      weight_type = weight_type,
      kernel = kernel, bw = bw, time_index = time_index,
      psd = psd,
      sw = sw_flag, ivar_vec = ivar_vec
    )

    }  # end of if (parsed$K1 > 0L) — endogenous-regressor diagnostics

    # Orthogonality test (J1)
    if (!is.null(orthog_cols) && length(orthog_cols) > 0L) {
      partialled_out <- setdiff(orthog_cols, colnames(parsed$Z))
      if (length(partialled_out) > 0L) {
        stop("Cannot test orthogonality of variables that were partialled out: ",
             paste0("'", partialled_out, "'", collapse = ", "),
             ". Remove these from `orthog` or `partial`.", call. = FALSE)
      }
      diagnostics$orthog <- .compute_orthog_test(
        Z = parsed$Z, X = parsed$X, y = parsed$y,
        residuals = fit$residuals, rss = fit$rss,
        weights = parsed$weights, cluster_vec = cluster_vec,
        vcov_type = effective_vcov_type,
        N = parsed$N, K = parsed$K, L = parsed$L,
        orthog_vars = orthog_cols, dofminus = dofminus,
        weight_type = weight_type,
        kernel = kernel, bw = bw, time_index = time_index,
        center = center, psd = psd,
        omega = fit$omega,
        sw = sw_flag, ivar_vec = ivar_vec
      )
    }

    # Redundancy test (P1) — computation only if !noid (validation already
    # done above). K1 = 0: nothing to test redundancy against — Stata
    # silently ignores redundant() here (gated at ivreg2.ado:1741); we warn
    # instead, following the package's explicit-request-ignored convention.
    if (!is.null(redundant) && length(redundant) > 0L && parsed$K1 == 0L) {
      warning("`redundant` ignored: the model has no endogenous regressors, ",
              "so there is nothing to test instrument redundancy against.",
              call. = FALSE)
    }
    if (!noid && !is.null(redundant) && length(redundant) > 0L &&
        parsed$K1 > 0L) {
      redundant_cols <- .expand_terms_to_colnames(
        redundant, parsed$excluded_names, parsed$excluded_colnames,
        parsed$excluded_assign
      )
      diagnostics$redundancy <- .compute_redundancy_test(
        X = parsed$X, Z = parsed$Z,
        weights = parsed$weights, cluster_vec = cluster_vec,
        vcov_type = effective_vcov_type,
        N = parsed$N, K1 = parsed$K1,
        endo_colnames = parsed$endo_colnames,
        excluded_colnames = parsed$excluded_colnames,
        redundant_vars = redundant_cols, dofminus = dofminus,
        weight_type = weight_type,
        kernel = kernel, bw = bw, time_index = time_index,
        center = center
      )
    }

    }  # end of if (is.null(b0)) — identification diagnostics block
  }
  if (length(diagnostics) == 0L) diagnostics <- NULL

  # Reduced-form regression. K1 = 0: skipped, matching Stata — the entire
  # RF/first-stage estimation block is gated on `endo1_ct > 0`
  # (ivreg2.ado:1256), and saverf on a `(=z)` model stores nothing
  # (e(rfeq) empty; verified by probe 2026-06-11). Note the y-on-Z reduced
  # form WOULD differ from the main y-on-X model (Z carries the surplus
  # instruments), but Stata never computes it, so there is no ground truth
  # to validate an extension against.
  reduced_form_result <- NULL
  if (parsed$is_iv && parsed$K1 > 0L && reduced_form != "none") {
    rf_depvar <- parsed$y_name
    reduced_form_result <- .compute_reduced_form(
      mode           = reduced_form,
      Z              = parsed$Z,
      X              = parsed$X,
      y              = parsed$y,
      weights        = parsed$weights,
      cluster_vec    = cluster_vec,
      vcov_type      = effective_vcov_type,
      N              = parsed$N,
      K              = parsed$K,
      L              = parsed$L,
      K1             = parsed$K1,
      L1             = parsed$L1,
      M              = M,
      endo_names     = parsed$endo_colnames,
      excluded_names = parsed$excluded_colnames,
      depvar_name    = rf_depvar,
      dofminus       = dofminus,
      sdofminus      = sdofminus,
      weight_type    = weight_type,
      kernel         = kernel,
      bw             = bw,
      time_index     = time_index,
      center         = center,
      sw             = sw_flag,
      ivar_vec       = ivar_vec,
      psd            = psd
    )
  }

  # Build extractable first-stage model objects (Ticket S1)
  first_stage_models <- NULL
  if (first_stage_flag && parsed$is_iv && !is.null(first_stage)) {
    first_stage_models <- .build_first_stage_models(
      fs_results   = first_stage,
      Z            = parsed$Z,
      weights      = parsed$weights,
      cluster_vec  = cluster_vec,
      vcov_type    = effective_vcov_type,
      N            = parsed$N,
      L            = parsed$L,
      M            = M,
      dofminus     = dofminus,
      sdofminus    = sdofminus,
      weight_type  = weight_type,
      kernel       = kernel,
      bw           = bw,
      time_index   = time_index,
      center       = center,
      cluster_var  = prep$cluster_var_name
    )
  }

  # Model F-test
  model_f_result <- .compute_model_f(
    coefficients  = fit$coefficients,
    vcov          = fit$vcov,
    N             = parsed$N,
    K             = parsed$K,
    has_intercept = parsed$has_intercept,
    vcov_type     = effective_vcov_type,
    small         = small,
    M             = M,
    dofminus      = dofminus,
    sdofminus     = sdofminus
  )

  list(
    diagnostics = diagnostics, first_stage = first_stage,
    first_stage_models = first_stage_models,
    reduced_form_result = reduced_form_result,
    model_f_result = model_f_result,
    effective_vcov_type = effective_vcov_type, center = center
  )
}


#' Extended Instrumental Variables Estimation
#'
#' Estimate models by OLS, two-stage least squares (2SLS), LIML, Fuller, or
#' k-class with automatic diagnostic tests. Uses a three-part formula for IV:
#' `y ~ exog | endo | instruments`.
#'
#' @param formula A formula: `y ~ exog` (OLS),
#'   `y ~ exog | endo | instruments` (IV), or
#'   `y ~ exog | 0 | instruments` (no endogenous regressors, with the
#'   excluded instruments contributing surplus moment conditions —
#'   Stata's `(=z1 z2)` form). The latter is estimated by OLS; the surplus
#'   conditions drive the Sargan/Hansen J overidentification test and
#'   `orthog` C-tests, and `method = "gmm2s"` with `vcov = "robust"` gives
#'   Cragg's (1983) heteroskedastic OLS (HOLS) estimator (with the default
#'   iid VCE the two-step weighting matrix is proportional to
#'   \eqn{(Z'Z)^{-1}} and the estimates equal OLS exactly).
#' @param data A data frame containing the variables in the formula.
#' @param weights Optional analytic weights expression (evaluated in `data`),
#'   equivalent to Stata's `[aw=varname]`. Must be strictly positive.
#'   For aweights and pweights, the weights are normalized internally to sum
#'   to N, following Stata's convention. This makes sigma (RMSE)
#'   scale-invariant for those types: multiplying all weights by a constant
#'   changes no coefficients, standard errors, or test statistics. Frequency
#'   and importance weights are *not* normalized --- they redefine N as
#'   `sum(weights)` (see `weight_type`), so their scale matters by
#'   construction.
#'
#'   **Note:** for aweights and pweights, sigma will differ from
#'   [lm()]`(..., weights = w)` by a factor
#'   of `sqrt(N / sum(w))` because `lm()` uses raw (unnormalized) weights.
#'   Coefficients are identical to `lm()` for OLS; SEs and the VCV matrix
#'   additionally match `lm()` when `small = TRUE` (the default
#'   `small = FALSE` reports Stata-convention SEs without the N − K
#'   finite-sample correction).
#' @param subset Optional subset expression (evaluated in `data`).
#' @param na.action Function for handling `NA`s (default [na.omit]).
#'   Use [na.exclude] to drop incomplete cases during estimation but
#'   reinsert `NA`s at the omitted positions in [fitted()], [residuals()],
#'   and [predict()] output, so results align with the original data frame.
#'   Stata drops missing observations silently (equivalent to [na.omit])
#'   and has no analogue of [na.exclude].
#' @param vcov Character: covariance type. One of `"iid"` (classical),
#'   `"robust"` (White heteroskedasticity-robust), `"HAC"`
#'   (heteroskedasticity and autocorrelation consistent), or `"AC"`
#'   (autocorrelation consistent). Use `small = TRUE` to apply
#'   finite-sample corrections. To match Stata's `ivreg2, robust`:
#'   use `vcov = "robust"`. To match Stata's `ivreg2, robust small`:
#'   use `vcov = "robust", small = TRUE`.
#' @param clusters One-sided formula specifying one or two cluster variables
#'   (e.g. `~ firmid` for one-way, `~ firmid + year` for two-way).
#'   Two-way clustering uses the Cameron-Gelbach-Miller (2006) formula.
#'   The effective cluster count is `min(M1, M2)` per Stata convention.
#'   The `small` argument controls whether the finite-sample correction
#'   `(N-1)/(N-K) * M/(M-1)` is applied (matching Stata's
#'   `cluster() small` combination).
#' @param endog Character vector of endogenous regressor names to test for
#'   exogeneity (endogeneity test / C-statistic). If `NULL` (default), tests
#'   all endogenous regressors. Names must match variables in the endogenous
#'   part of the formula. Ignored for OLS models.
#'
#'   **Note:** Unlike Stata's `ivreg2`, which only computes the endogeneity
#'   test when the `endog()` option is explicitly specified, `ivreg2r` computes
#'   it automatically for all IV models.
#' @param orthog Character vector of instrument names to test for
#'   orthogonality (instrument-subset C-statistic). Names must be included
#'   or excluded instruments (not endogenous regressors or the intercept).
#'   If `NULL` (default), no orthogonality test is computed. Ignored for
#'   OLS models. Equivalent to Stata's `orthog()` option.
#' @param redundant Character vector of excluded instrument names to test
#'   for redundancy (zero first-stage explanatory power). The test is a
#'   KP rk LM test of H0: rank=0 on the first-stage coefficient matrix
#'   for the tested instruments, conditional on maintained instruments.
#'   If `NULL` (default), no redundancy test is computed. Ignored for OLS
#'   models. Equivalent to Stata's `redundant()` option.
#' @param small Logical: if `TRUE`, use small-sample corrections
#'   (t/F instead of z/chi-squared, `N-K` denominator for sigma).
#' @param weight_type Character: type of weights. One of `"aweight"`
#'   (analytic weights, default), `"fweight"` (frequency weights),
#'   `"pweight"` (probability/sampling weights), or `"iweight"`
#'   (importance weights).
#'
#'   **aweight**: Normalized to sum to N. Standard for WLS.
#'
#'   **fweight**: Integer-valued; N is redefined as `sum(weights)`.
#'   HC meat uses `w * e^2` (linear, not quadratic).
#'
#'   **pweight**: Normalized to sum to N. Forces robust VCE
#'   (overrides `vcov = "iid"` to `"robust"` and `vcov = "AC"` ---
#'   explicit or kiefer-implied --- to `"HAC"`).
#'
#'   **iweight**: Not normalized; N is redefined as `sum(weights)` (float).
#'   All intermediate calculations use float N; posted `nobs` and `df.residual`
#'   use `floor(N)`. Restricted to IID VCE only — incompatible with robust,
#'   cluster, HAC, and GMM estimation.
#' @param dofminus Non-negative integer: large-sample degrees-of-freedom
#'   adjustment. Subtracted from N in large-sample variance formulas
#'   (e.g., sigma = rss/(N-dofminus)). Useful when fixed effects have been
#'   partialled out. Equivalent to Stata's `dofminus()` option.
#' @param sdofminus Non-negative integer: small-sample degrees-of-freedom
#'   adjustment. Subtracted from the residual degrees of freedom alongside K
#'   (e.g., df.residual = N - K - dofminus - sdofminus). Useful when
#'   partialling out regressors. Equivalent to Stata's `sdofminus()` option.
#' @param method Character: estimation method. One of `"2sls"` (default),
#'   `"liml"`, `"kclass"`, `"gmm2s"` (two-step efficient GMM), or `"cue"`
#'   (continuously updated GMM estimator).
#'   For OLS models (1-part formula), this is ignored. When `fuller > 0`
#'   is specified, method is automatically promoted to `"liml"`.
#'   `"gmm2s"` uses the inverse of the moment covariance matrix as the
#'   optimal weighting matrix, yielding efficient estimates under the
#'   specified error structure. Incompatible with `fuller` and `kclass`.
#'   `"cue"` continuously updates both moment conditions and moment
#'   covariance at each candidate beta during optimization. Under iid errors
#'   (no robust/cluster/kernel VCE), CUE with no other options is equivalent
#'   to LIML with `coviv = TRUE`. Under non-iid errors, CUE generalizes
#'   LIML to produce efficient estimates.
#'   Incompatible with `fuller`, `kclass`, `wmatrix`, and `smatrix`.
#'
#'   **CUE optimizer note:** This package uses R's `optim()` (BFGS +
#'   Nelder-Mead) while Stata uses Mata's `optimize()` (modified
#'   Newton-Raphson). The CUE objective is non-convex and can have
#'   multiple local minima, so the two optimizers may converge to
#'   different solutions. In rare cases, Stata's Newton-Raphson can
#'   traverse indefinite-Hessian regions and land in pathological basins
#'   (e.g., negative R-squared). When results differ, compare the
#'   J-statistic and R-squared to assess which solution is more sensible.
#' @param kclass Numeric scalar: user-supplied k value for k-class
#'   estimation. When supplied, `method` is automatically set to `"kclass"`.
#'   Must be non-negative. Cannot be combined with `method = "liml"` or
#'   `fuller`.
#' @param fuller Numeric scalar: Fuller (1977) modification parameter.
#'   Must be positive. When supplied, `method` is automatically set to
#'   `"liml"` and `k = lambda - fuller / (N - L)`. `fuller = 1` gives the
#'   bias-corrected LIML estimator; `fuller = 4` targets MSE. Cannot be
#'   combined with `kclass`.
#' @param coviv Logical: if `TRUE`, use the 2SLS bread `(X_hat'X_hat)^{-1}`
#'   instead of the k-class bread for VCV computation in LIML/k-class
#'   estimation. This gives the "COVIV" (covariance at the IV estimates)
#'   VCV that is robust to misspecification of the LIML model. Silently
#'   ignored for OLS and 2SLS. Default `FALSE`.
#' @param kernel Character: kernel function for HAC/AC standard errors.
#'   One of `"bartlett"`, `"parzen"`, `"truncated"`, `"tukey-hanning"`,
#'   `"tukey-hamming"`, `"qs"` (quadratic spectral), `"daniell"`, or
#'   `"tent"`. When specified without an explicit `vcov` change, `vcov`
#'   is automatically set: `"iid"` becomes `"AC"`, `"robust"` becomes
#'   `"HAC"`. Requires `bw` and `tvar`.
#'
#'   **Note:** Kleibergen-Paap identification tests (underidentification and
#'   weak identification) always use the Bartlett kernel internally, regardless
#'   of the user-specified kernel. This matches a known behavior in Stata's
#'   `ivreg2`, where `ranktest` hard-codes the Bartlett kernel.
#' @param bw Numeric or `"auto"`: bandwidth for kernel estimation. Must be
#'   positive numeric, or `"auto"` for automatic selection via Newey-West
#'   (1994). Auto selection is available for Bartlett, Parzen, and Quadratic
#'   Spectral kernels only, and is not supported for panel data (`ivar`).
#'   Required when `kernel` is specified.
#' @param tvar Character: name of the time variable in `data`. Required
#'   for HAC/AC estimation.
#' @param ivar Character: name of the panel identifier variable in `data`.
#'   Optional; only needed for panel data with HAC/AC estimation.
#' @param kiefer Logical: if `TRUE`, use the Kiefer (1980) VCE —
#'   autocorrelation-consistent with kernel = Truncated and bandwidth = T
#'   (the full time span). Requires panel data (`tvar` + `ivar`).
#'   Incompatible with robust VCE, clustering, explicit kernel or bandwidth.
#'   Equivalent to Stata's `ivreg2 ..., kiefer`.
#' @param dkraay Positive numeric scalar: bandwidth for Driscoll-Kraay (1998)
#'   VCE. When specified, clusters on the time variable and applies kernel
#'   smoothing across time lags, producing standard errors robust to
#'   cross-sectional dependence. Requires panel data (`tvar` + `ivar`).
#'   Incompatible with explicit `bw`. If `kernel` is not specified, defaults
#'   to Bartlett. Equivalent to Stata's `ivreg2 ..., dkraay(3)`.
#' @param wmatrix Numeric matrix: user-supplied L x L weighting matrix for
#'   GMM estimation, where L is the number of instruments (including exogenous
#'   regressors). When supplied without `method = "gmm2s"`, produces
#'   inefficient GMM with full sandwich VCV (`method` becomes `"gmmw"`).
#'   When combined with `method = "gmm2s"`, used as the first-step weighting
#'   matrix (second step uses the optimal weighting matrix). Must be symmetric.
#'   Incompatible with `method = "liml"`, `kclass`, and `fuller`.
#'   When `method` is not `"gmm2s"`, ignored without robust VCE, clustering,
#'   or HAC (with a warning).
#'   Equivalent to Stata's `wmatrix()` option. When used (not ignored), it is
#'   echoed back as `fit$W`, matching Stata's `e(W)`. A matrix with dimnames
#'   is matched to the instrument columns by name (reordered as needed,
#'   mirroring Stata's matsort; names that do not cover the instrument list
#'   are an error); an unnamed matrix is used positionally, in instrument
#'   column order.
#' @param smatrix Numeric matrix: user-supplied L x L moment covariance matrix
#'   for GMM estimation. When supplied, the GMM estimation uses this matrix
#'   instead of computing Omega from residuals. Implies efficient GMM
#'   (`method` is promoted to `"gmm2s"` if `"2sls"`). Must be symmetric.
#'   Incompatible with `method = "liml"`, `kclass`, and `fuller`.
#'   Equivalent to Stata's `smatrix()` option. Echoed back as `fit$S`
#'   (never recomputed), matching Stata's `e(S)`. A matrix with dimnames is
#'   matched to the instrument columns by name (reordered as needed,
#'   mirroring Stata's matsort); an unnamed matrix is used positionally, in
#'   instrument column order.
#' @param b0 Numeric vector: evaluate the CUE objective at this fixed
#'   parameter vector without optimization. When supplied, `method` is
#'   promoted to `"cue"` and the J(b0) statistic is computed and stored.
#'   Identification diagnostics (underid, weak-id, first-stage, AR, SW,
#'   endogeneity, orthog) are suppressed (matching Stata's `b0()` option).
#'   Length must equal the number of regressors K. Named vectors are
#'   reordered to match model matrix columns; unnamed vectors are used
#'   as-is in model matrix column order.
#'   Incompatible with `method = "gmm2s"`, `"liml"`, `kclass`, `fuller`,
#'   and `wmatrix`.
#' @param noid Logical: if `TRUE`, suppress computation of underidentification,
#'   weak identification, and redundancy test statistics. Overidentification,
#'   Anderson-Rubin, Stock-Wright, endogeneity, and first-stage diagnostics
#'   are still computed. Default `FALSE`. Useful for speeding up simulation
#'   loops where rank-based identification tests are expensive.
#'   Equivalent to Stata's `noid` option.
#' @param partial Character vector: exogenous regressors to partial out
#'   via Frisch-Waugh-Lovell projection before estimation. Coefficients on
#'   partialled variables are not recoverable and are not reported.
#'   Special values:
#'   - `"_cons"` or `"(Intercept)"`: partial only the constant (demean).
#'   - `"_all"`: partial all included exogenous regressors.
#'   - `NULL` (default): no partialling.
#'
#'   By default, `sdofminus` is automatically incremented by the number of
#'   partialled variables (including the constant). Use `nopartialsmall = TRUE`
#'   to suppress this adjustment.
#'
#'   **Note:** FWL invariance holds for OLS, 2SLS, LIML, and two-step GMM, but
#'   **not** for CUE. The CUE objective recomputes the moment covariance at each
#'   iteration, so partialled residuals produce a different optimization surface.
#'   A warning is issued when `partial` is used with `method = "cue"`.
#'
#'   After partialling, `predict()` is restricted to residuals only (no
#'   `newdata`). Summary output notes that total SS, model F, and R-squared
#'   are partial-model values.
#'
#'   Equivalent to Stata's `partial()` option.
#' @param nopartialsmall Logical: if `TRUE`, suppress the automatic
#'   `sdofminus` adjustment from partialling. Default `FALSE`.
#'   Equivalent to Stata's `nopartialsmall` option.
#' @param center Logical: if `TRUE`, subtract the mean of moment conditions
#'   (scores) before computing the S matrix (meat of the sandwich VCE).
#'   Default `FALSE`. Centering only affects non-homoskedastic VCE types
#'   (HC, cluster, HAC); a warning is issued if used with IID or AC VCE.
#'
#'   **Note:** Centering is applied to the main model's VCE and to
#'   diagnostic tests that use the main model's S matrix (overidentification,
#'   orthogonality). However, the endogeneity test (`endog`) computes its own
#'   S matrix from the restricted model *without* centering, even when
#'   `center = TRUE`. This matches Stata's `ivreg2`, where `center` is not
#'   forwarded to the recursive call for the endogeneity test.
#' @param psd Character or NULL: PSD correction for the moment covariance
#'   matrix S. `NULL` (default) applies no correction. `"psd0"` zeroes
#'   negative eigenvalues (Politis 2007). `"psda"` replaces negative
#'   eigenvalues with their absolute values (Stock & Watson 2008). The
#'   correction is applied to S *before* the variance-covariance matrix is
#'   assembled from it, matching Stata's `m_omega`; the corrected S is
#'   stored as `$S`. The conventional (`vcov = "iid"`) VCV and the
#'   Kleibergen-Paap identification statistics are never corrected, also
#'   matching Stata (its `ranktest` does not receive the psd option). A
#'   warning is emitted whenever negative eigenvalues are detected and
#'   corrected (Stata corrects silently). Ignored when a user-supplied
#'   `smatrix` is given: the user matrix is used as-is for estimation, the
#'   VCV, and `$S`, matching Stata (whose `m_omega` is never invoked for a
#'   supplied S). Equivalent to Stata's `psd0` and `psda` options.
#' @param sw Logical: if `TRUE`, compute the Stock-Watson (2008,
#'   Econometrica) panel-robust VCE. Requires panel data (`ivar`).
#'   Incompatible with clustering, HAC kernels, kiefer, dkraay,
#'   fweights, and iweights. Forces `vcov = "robust"` automatically.
#'   Equivalent to Stata's `sw` option (labeled "BETA VERSION" in Stata).
#' @param reduced_form Character: what reduced-form output to store.
#'   `"none"` (default) stores nothing. `"rf"` stores the y ~ Z regression
#'   (equivalent to Stata's `saverf`; Stata's `rf` option *displays* the
#'   reduced form without storing it). `"system"` stores the full system of
#'   y + all endogenous variables regressed on Z, with cross-equation VCV
#'   (equivalent to Stata's `savesfirst`, an option present in `ivreg2.ado`
#'   but absent from its help file; Stata's `sfirst` displays the system
#'   without storing it). Silently ignored for OLS models and for
#'   empty-endogenous (`y ~ exog | 0 | instruments`) models, matching Stata,
#'   which skips reduced-form estimation whenever the endogenous list is
#'   empty.
#' @param first_stage Logical: if `TRUE`, store extractable first-stage
#'   regression objects on the fitted model. Access them via
#'   [first_stage()]. Each object supports [coef()], [vcov()],
#'   [summary()], [tidy()], and [glance()]. Equivalent to Stata's
#'   `savefirst` option. Default `FALSE`. Silently ignored for OLS models.
#' @param model Logical: if `TRUE` (default), store the model frame in the
#'   return object.
#' @param x Logical: if `TRUE`, store model matrices (`X`, `Z`) in the
#'   return object.
#' @param y Logical: if `TRUE` (default), store the response vector in the
#'   return object.
#'
#' @return An object of class `"ivreg2"`. Beyond the usual components
#'   (`coefficients`, `residuals`, `vcov`, `diagnostics`, ...), the object
#'   stores the estimated moment-condition covariance matrix as `$S` and the
#'   GMM weighting matrix as `$W`, corresponding to Stata's `e(S)` and
#'   `e(W)`:
#'
#'   - `$S` is present for every fit. Normalization matches Stata's
#'     `m_omega`: iid `sigma^2 (Z'Z)/N` with `sigma^2 = RSS/(N - dofminus)`;
#'     robust and HAC, meat divided by `N - dofminus`; cluster, meat divided
#'     by `N`; psd-corrected when `psd` is set; built from centered moments
#'     when `center = TRUE`. A user-supplied `smatrix` is echoed, never
#'     recomputed and never psd-corrected.
#'   - `$W` is `NULL` for `method = "liml"` and `"kclass"` (as in Stata, no
#'     weighting matrix is defined for k-class estimators). For two-step
#'     GMM, CUE, and `b0` fits it is the inverse of the defining `$S`; for
#'     one-step OLS/2SLS it is `(1/sigma^2) (Z'Z/N)^{-1}`; a user-supplied
#'     `wmatrix` is echoed (it is the first-step weighting matrix under
#'     `method = "gmm2s"`).
#'
#'   Rows and columns are named by the R instrument columns — `(Intercept)`
#'   rather than Stata's `_cons`, and R column order (exogenous regressors
#'   before excluded instruments) rather than Stata's. Passing `fit$S` back
#'   via `smatrix` reproduces the corresponding efficient-GMM fit, mirroring
#'   the `e(S)` reuse workflow in Stata's `ivreg2` help file.
#'
#' @examples
#' data(mroz)
#' mroz_work <- subset(mroz, inlf == 1)
#'
#' # --- Just-identified IV: return to schooling ---
#' # Wooldridge (2020), Example 15.1: father's education instruments
#' # education in a simple wage equation (replicates the published
#' # estimates: educ = 0.059, SE = 0.035). With one instrument for one
#' # endogenous variable, no overid test is possible; instrument validity
#' # must be defended on substantive grounds.
#' fit <- ivreg2(lwage ~ 1 | educ | fatheduc,
#'               data = mroz_work)
#' summary(fit)
#'
#' # --- Overidentified IV: testing instrument validity ---
#' # The Stata ivreg2 help file's illustrative baseline (line 1274):
#' # age, kidslt6, and kidsge6 give two overidentifying restrictions,
#' # so the Sargan test can detect misspecification. Note its weak
#' # first stage (Cragg-Donald F ~ 4.3) -- see the LIML example below.
#' fit_overid <- ivreg2(lwage ~ exper + expersq | educ |
#'                        age + kidslt6 + kidsge6, data = mroz_work)
#' summary(fit_overid)
#'
#' # --- Robust standard errors ---
#' # With heteroskedasticity, Kleibergen-Paap diagnostics replace
#' # Anderson/Cragg-Donald, and Hansen J replaces Sargan.
#' fit_robust <- ivreg2(lwage ~ exper + expersq | educ |
#'                        age + kidslt6 + kidsge6, data = mroz_work,
#'                        vcov = "robust")
#' summary(fit_robust)
#'
#' \donttest{
#' # --- LIML ---
#' # Same baseline as above (a documented variation on the help-file
#' # spec). Its weak first stage (Cragg-Donald F ~ 4.3) is the setting
#' # where 2SLS bias toward OLS grows with the overidentification degree
#' # relative to instrument strength, and LIML's approximate
#' # median-unbiasedness -- established within the iid Stock-Yogo
#' # framework -- is attractive. LIML equals 2SLS when just-identified.
#' fit_liml <- ivreg2(lwage ~ exper + expersq | educ |
#'                      age + kidslt6 + kidsge6, data = mroz_work,
#'                      method = "liml")
#' summary(fit_liml)
#'
#' # --- Endogeneity test ---
#' # Is education actually endogenous? The C-statistic tests
#' # H0: OLS is consistent (educ is exogenous).
#' fit_endog <- ivreg2(lwage ~ exper + expersq | educ |
#'                       age + kidslt6 + kidsge6, data = mroz_work,
#'                       endog = "educ")
#' fit_endog$diagnostics$endogeneity
#'
#' # --- Clustering ---
#' # Griliches (1976) wage equation, cluster on year.
#' # Adapted from Stata `ivreg2` help-file example (line 1250): the help
#' # command also includes year dummies (xi i.year), omitted here, and no
#' # small option. Either way the fit warns -- 7 year clusters < the
#' # instrument count makes the moment covariance singular, which is the
#' # help file's own lead-in to partialling; see
#' # vignette("time-series-gmm").
#' data(griliches)
#' fit_cl <- ivreg2(lw ~ s + expr + tenure + rns + smsa |
#'                    iq | med + kww + age,
#'                  data = griliches, clusters = ~year, small = TRUE)
#' summary(fit_cl)
#' }
#'
#' @export
ivreg2 <- function(formula, data, weights, subset, na.action = stats::na.omit,
                   vcov = "iid", clusters = NULL, endog = NULL,
                   orthog = NULL,
                   redundant = NULL,
                   method = "2sls", kclass = NULL, fuller = 0,
                   coviv = FALSE,
                   small = FALSE,
                   dofminus = 0L, sdofminus = 0L,
                   weight_type = "aweight",
                   kernel = NULL, bw = NULL, tvar = NULL, ivar = NULL,
                   kiefer = FALSE, dkraay = NULL,
                   wmatrix = NULL, smatrix = NULL,
                   b0 = NULL,
                   noid = FALSE,
                   partial = NULL,
                   nopartialsmall = FALSE,
                   center = FALSE,
                   psd = NULL,
                   sw = FALSE,
                   reduced_form = "none",
                   first_stage = FALSE,
                   model = TRUE, x = FALSE, y = TRUE) {

  # --- 1. Capture call ---
  cl <- match.call()

  # --- 2. Validate and normalize arguments ---
  opts <- .validate_and_normalize_args(
    vcov = vcov, clusters = clusters, endog = endog, orthog = orthog,
    redundant = redundant, method = method, kclass = kclass, fuller = fuller,
    coviv = coviv, small = small, dofminus = dofminus, sdofminus = sdofminus,
    weight_type = weight_type, kernel = kernel, bw = bw, tvar = tvar,
    ivar = ivar, kiefer = kiefer, dkraay = dkraay,
    wmatrix = wmatrix, smatrix = smatrix, b0 = b0,
    partial = partial, nopartialsmall = nopartialsmall,
    center = center, psd = psd, reduced_form = reduced_form,
    first_stage = first_stage, noid = noid, sw = sw,
    model_flag = model, x_flag = x, y_flag = y
  )
  # Unpack options needed in ivreg2() scope (for .new_ivreg2() assembly).
  # Helpers read from `opts` directly.
  small       <- opts$small
  dofminus    <- opts$dofminus
  coviv       <- opts$coviv
  weight_type <- opts$weight_type
  kernel      <- opts$kernel
  kiefer      <- opts$kiefer
  dkraay      <- opts$dkraay
  tvar        <- opts$tvar
  ivar        <- opts$ivar
  psd         <- opts$psd
  smatrix     <- opts$smatrix

  # --- 3. Forward to parser ---
  # Build a call to .parse_formula() using the NSE arguments from ivreg2().
  # Evaluate in an environment where .parse_formula is visible (it's unexported)
  # with the user's frame as parent so data/subset/weights expressions resolve.
  pf_call <- cl[c(1L, match(c("formula", "data", "weights", "subset",
                               "na.action"), names(cl), 0L))]
  pf_call[[1L]] <- quote(.parse_formula)
  pf_env <- list2env(list(.parse_formula = .parse_formula),
                     parent = parent.frame())
  parsed <- eval(pf_call, pf_env)

  # --- 3a. Prepare model ---
  prep <- .prepare_model(parsed, opts, data, clusters)
  parsed       <- prep$parsed
  sdofminus    <- prep$sdofminus
  b0           <- prep$b0
  # NOTE: `smatrix` stays the RAW user matrix for the backward-compat
  # fit$smatrix echo field; the name-normalized copy (prep$smatrix) is what
  # estimation and fit$S consume.
  cluster_vec  <- prep$cluster_vec
  cluster_var_name <- prep$cluster_var_name
  M            <- prep$M
  M1           <- prep$M1
  M2           <- prep$M2
  time_index   <- prep$time_index
  unsort_order <- prep$unsort_order
  partial_ct   <- prep$partial_ct
  partial_names <- prep$partial_names
  partialcons  <- prep$partialcons
  n_physical   <- prep$n_physical
  w_raw        <- prep$w_raw
  bw           <- prep$bw

  # --- 4. Resolve estimation path and fit ---
  est <- .resolve_and_fit(parsed, prep, opts)
  fit      <- est$fit
  method   <- est$method
  bw       <- est$bw
  # Backward-compat fit$wmatrix echo keeps the RAW user matrix (NULL when
  # the wmatrix was ignored-with-warning); the name-normalized copy
  # (est$wmatrix) feeds estimation and fit$W.
  wmatrix  <- if (is.null(est$wmatrix)) NULL else opts$wmatrix
  omega_fn <- est$omega_fn

  # --- 5. VCV ---
  fit <- .compute_vcov(fit, parsed, prep, opts, method, bw)

  # --- 5b. Diagnostics, reduced-form, model F ---
  diag_result <- .compute_diagnostics(fit, parsed, prep, opts, method, bw)
  diagnostics         <- diag_result$diagnostics
  first_stage         <- diag_result$first_stage
  first_stage_models  <- diag_result$first_stage_models
  reduced_form_result <- diag_result$reduced_form_result
  model_f_result      <- diag_result$model_f_result
  effective_vcov_type <- diag_result$effective_vcov_type
  center              <- diag_result$center

  # --- 5c. Stored moment-covariance S and weighting matrix W ---
  # Must precede the unsort: HAC/AC moment covariances index the residuals
  # via the sorted prep$time_index. Uses the effective center from the
  # diagnostics step (reset for iid/AC, matching m_omega).
  stored_sw <- .compute_stored_sw(fit, parsed, prep, opts, est,
                                  bw = bw, center = center)

  # --- 5d. Unsort for user-facing output ---
  # If data was sorted for HAC/AC, restore original row order for user-facing
  # vectors and matrices (residuals, fitted values, y, X, Z).
  if (!is.null(unsort_order)) {
    fit$residuals <- fit$residuals[unsort_order]
    fit$fitted.values <- fit$fitted.values[unsort_order]
    parsed$y <- parsed$y[unsort_order]
    parsed$X <- parsed$X[unsort_order, , drop = FALSE]
    parsed$Z <- parsed$Z[unsort_order, , drop = FALSE]
    # Also unsort first-stage model residuals/fitted values
    if (!is.null(first_stage_models)) {
      for (nm in names(first_stage_models)) {
        first_stage_models[[nm]]$residuals <-
          first_stage_models[[nm]]$residuals[unsort_order]
        first_stage_models[[nm]]$fitted.values <-
          first_stage_models[[nm]]$fitted.values[unsort_order]
      }
    }
  }

  # --- 5e. Compute additional stored results (R2) ---
  # Condition numbers: Stata's cond(XX) = max(eigenvalue)/min(eigenvalue)
  # which equals kappa(XX, exact = TRUE) in R.
  if (is.null(parsed$weights)) {
    XX <- crossprod(parsed$X)
  } else {
    XX <- crossprod(sqrt(parsed$weights) * parsed$X)
  }
  condxx <- kappa(XX, exact = TRUE)

  if (parsed$is_iv) {
    if (is.null(parsed$weights)) {
      ZZ <- crossprod(parsed$Z)
    } else {
      ZZ <- crossprod(sqrt(parsed$weights) * parsed$Z)
    }
    condzz <- kappa(ZZ, exact = TRUE)
  } else {
    # OLS: Z = X, so condzz = condxx (matches Stata behavior)
    condzz <- condxx
  }

  # Gaussian log-likelihood (Stata line 2009)
  ll <- -0.5 * (parsed$N * log(2 * pi) + parsed$N * log(fit$rss / parsed$N) +
                 parsed$N)

  # --- 6. Assemble return object ---
  # Determine effective method for the return object
  est_method <- if (method == "gmmw") {
    "gmmw"
  } else if (method %in% c("liml", "kclass", "gmm2s", "cue")) {
    method
  } else if (parsed$is_iv && parsed$K1 > 0L) {
    "2sls"
  } else {
    # Includes the empty-endogenous (K1 = 0) form: Stata posts
    # e(model) = "ols" (ivreg2.ado:2101-2106)
    "ols"
  }

  .new_ivreg2(
    coefficients  = fit$coefficients,
    residuals     = fit$residuals,
    fitted.values = fit$fitted.values,
    vcov          = fit$vcov,
    sigma         = fit$sigma,
    df.residual   = fit$df.residual,
    rank          = fit$rank,
    r.squared     = fit$r.squared,
    adj.r.squared = fit$adj.r.squared,
    rss           = fit$rss,
    r2u           = fit$r2u,
    r2c           = fit$r2c,
    mss           = fit$mss,
    model_f       = model_f_result$model_f,
    model_f_p     = model_f_result$model_f_p,
    model_f_df1   = model_f_result$model_f_df1,
    model_f_df2   = model_f_result$model_f_df2,
    diagnostics   = diagnostics,
    first_stage   = first_stage,
    first_stage_models = first_stage_models,
    reduced_form  = reduced_form_result,
    call          = cl,
    formula       = parsed$formula,
    terms         = parsed$terms,
    nobs          = if (weight_type == "iweight") as.integer(floor(parsed$N)) else parsed$N,
    vcov_type     = effective_vcov_type,
    small         = small,
    dofminus      = dofminus,
    sdofminus     = sdofminus,
    cluster_var   = cluster_var_name,
    n_clusters    = M,
    n_clusters1   = M1,
    n_clusters2   = M2,
    na.action     = parsed$na.action,
    weights       = w_raw,
    weight_type   = weight_type,
    n_physical    = if (weight_type %in% c("fweight", "iweight")) n_physical else NULL,
    endogenous    = parsed$endo_names,
    endo_colnames = parsed$endo_colnames,
    instruments   = parsed$excluded_names,
    dropped_regressors      = parsed$dropped_regressors,
    dropped_instruments     = parsed$dropped_instruments,
    reclassified_endogenous = parsed$reclassified_endogenous,
    contrasts     = parsed$contrasts,
    xlevels       = parsed$xlevels,
    method            = est_method,
    lambda            = fit$lambda %||% NA_real_,
    kclass_value      = fit$kclass_value %||% NA_real_,
    fuller_parameter  = fit$fuller_param %||% 0,
    coviv             = coviv,
    S                 = stored_sw$S,
    W                 = stored_sw$W,
    wmatrix           = wmatrix,
    smatrix           = smatrix,
    b0                = b0,
    noid              = noid,
    cue_convergence   = fit$convergence,   # NULL for non-CUE
    cue_message       = fit$cue_message,   # NULL for non-CUE
    kernel            = kernel,
    bw                = bw,
    tvar              = tvar,
    kiefer            = kiefer,
    dkraay            = dkraay,
    ivar              = ivar,
    center            = center,
    psd               = psd,
    sw                = isTRUE(opts$sw),
    partial_ct        = partial_ct,
    partial_names     = partial_names,
    partialcons       = partialcons,
    yy                = fit$tss_u,
    yyc               = fit$tss_c,
    rankzz            = if (parsed$is_iv) parsed$L else fit$rank,
    condxx            = condxx,
    condzz            = condzz,
    ll                = ll,
    model         = if (model) parsed$model_frame else NULL,
    x             = if (x) list(X = parsed$X, Z = parsed$Z) else NULL,
    y             = if (y) parsed$y else NULL
  )
}
