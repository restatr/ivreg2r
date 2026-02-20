# --------------------------------------------------------------------------
# ivreg2
# --------------------------------------------------------------------------
#' Extended Instrumental Variables Estimation
#'
#' Estimate models by OLS, two-stage least squares (2SLS), LIML, Fuller, or
#' k-class with automatic diagnostic tests. Uses a three-part formula for IV:
#' `y ~ exog | endo | instruments`.
#'
#' @param formula A formula: `y ~ exog` (OLS) or
#'   `y ~ exog | endo | instruments` (IV).
#' @param data A data frame containing the variables in the formula.
#' @param weights Optional analytic weights expression (evaluated in `data`),
#'   equivalent to Stata's `[aw=varname]`. Must be strictly positive.
#'   Weights are normalized internally to sum to N, following Stata's
#'   convention. This makes sigma (RMSE) scale-invariant: multiplying all
#'   weights by a constant does not change sigma. Coefficients, standard
#'   errors, and all test statistics are unaffected by scale.
#'
#'   **Note:** sigma will differ from [lm()]`(..., weights = w)` by a factor
#'   of `sqrt(N / sum(w))` because `lm()` uses raw (unnormalized) weights.
#'   Coefficients, SEs, and the VCV matrix are identical to `lm()` for OLS.
#' @param subset Optional subset expression (evaluated in `data`).
#' @param na.action Function for handling `NA`s (default [na.omit]).
#' @param vcov Character: covariance type. One of `"iid"` (classical),
#'   `"HC0"` (White robust, no finite-sample correction), or `"HC1"`
#'   (White robust with N/(N-K) finite-sample correction). To match
#'   Stata's `ivreg2, robust`: use `"HC0"`. To match
#'   Stata's `ivreg2, robust small`: use `"HC1"` with `small = TRUE`.
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
#' @param orthog Character vector of instrument names to test for
#'   orthogonality (instrument-subset C-statistic). Names must be included
#'   or excluded instruments (not endogenous regressors or the intercept).
#'   If `NULL` (default), no orthogonality test is computed. Ignored for
#'   OLS models. Equivalent to Stata's `orthog()` option.
#' @param small Logical: if `TRUE`, use small-sample corrections
#'   (t/F instead of z/chi-squared, `N-K` denominator for sigma).
#' @param dofminus Non-negative integer: large-sample degrees-of-freedom
#'   adjustment. Subtracted from N in large-sample variance formulas
#'   (e.g., sigma = rss/(N-dofminus)). Useful when fixed effects have been
#'   partialled out. Equivalent to Stata's `dofminus()` option.
#' @param sdofminus Non-negative integer: small-sample degrees-of-freedom
#'   adjustment. Subtracted from the residual degrees of freedom alongside K
#'   (e.g., df.residual = N - K - dofminus - sdofminus). Useful when
#'   partialling out regressors. Equivalent to Stata's `sdofminus()` option.
#' @param method Character: estimation method. One of `"2sls"` (default),
#'   `"liml"`, or `"kclass"`. For OLS models (1-part formula), this is
#'   ignored. When `fuller > 0` is specified, method is automatically
#'   promoted to `"liml"`.
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
#' @param reduced_form Character: what reduced-form output to store.
#'   `"none"` (default) stores nothing. `"rf"` stores the y ~ Z regression
#'   (equivalent to Stata's `saverf`). `"system"` stores the full system of
#'   y + all endogenous variables regressed on Z, with cross-equation VCV
#'   (equivalent to Stata's `savesfirst`). Silently ignored for OLS models.
#' @param model Logical: if `TRUE` (default), store the model frame in the
#'   return object.
#' @param x Logical: if `TRUE`, store model matrices (`X`, `Z`) in the
#'   return object.
#' @param y Logical: if `TRUE` (default), store the response vector in the
#'   return object.
#'
#' @return An object of class `"ivreg2"`.
#'
#' @examples
#' fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
#' print(fit)
#' coef(fit)
#'
#' @export
ivreg2 <- function(formula, data, weights, subset, na.action = stats::na.omit,
                   vcov = "iid", clusters = NULL, endog = NULL,
                   orthog = NULL,
                   method = "2sls", kclass = NULL, fuller = 0,
                   coviv = FALSE,
                   small = FALSE,
                   dofminus = 0L, sdofminus = 0L,
                   reduced_form = "none",
                   model = TRUE, x = FALSE, y = TRUE) {

  # --- 1. Capture call ---
  cl <- match.call()

  # --- 2. Validate arguments ---
  if (!is.character(vcov) || length(vcov) != 1L) {
    stop("`vcov` must be a single character string.", call. = FALSE)
  }
  valid_vcov <- c("iid", "HC0", "HC1")
  if (!vcov %in% valid_vcov) {
    stop('vcov = "', vcov, '" is not yet implemented. ',
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

  if (!is.logical(coviv) || length(coviv) != 1L || is.na(coviv)) {
    stop("`coviv` must be TRUE or FALSE.", call. = FALSE)
  }

  valid_rf <- c("none", "rf", "system")
  if (!is.character(reduced_form) || length(reduced_form) != 1L ||
      !reduced_form %in% valid_rf) {
    stop('`reduced_form` must be one of "none", "rf", or "system".',
         call. = FALSE)
  }

  # --- 2b. Validate method / kclass / fuller ---
  if (!is.character(method) || length(method) != 1L) {
    stop("`method` must be a single character string.", call. = FALSE)
  }
  method <- tolower(method)
  valid_methods <- c("2sls", "liml", "kclass")
  if (!method %in% valid_methods) {
    stop('`method` must be one of "2sls", "liml", or "kclass".',
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
  # fuller implies liml
  if (fuller > 0 && method != "liml") {
    method <- "liml"
  }
  # kclass supplied implies method = "kclass"
  if (!is.null(kclass) && method != "kclass") {
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

  # --- 3a. Validate dofminus/sdofminus against model dimensions ---
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

  # --- 3b. Validate method against parsed model ---
  if (method %in% c("liml", "kclass") && !parsed$is_iv) {
    stop('`method = "', method, '"` requires an IV model (3-part formula).',
         call. = FALSE)
  }
  if (fuller > 0 && parsed$is_iv && fuller >= (parsed$N - parsed$L)) {
    stop("`fuller` (", fuller, ") must be less than N - L (",
         parsed$N - parsed$L, ").", call. = FALSE)
  }
  # --- 3c. Validate and normalize weights ---
  w_raw <- parsed$weights
  if (!is.null(w_raw) && any(!is.finite(w_raw)))
    stop("Weights must be finite and non-missing.", call. = FALSE)
  if (!is.null(w_raw) && any(w_raw <= 0))
    stop("All weights must be strictly positive.", call. = FALSE)
  # Normalize weights to sum to N (Stata aweight convention).
  # Makes sigma scale-invariant; the factor cancels in VCV so SEs are unchanged.
  # Raw (user-supplied) weights are stored in the return object.
  if (!is.null(w_raw)) {
    parsed$weights <- w_raw * (parsed$N / sum(w_raw))
  }

  # --- 3b. Parse clusters ---
  cluster_vec <- NULL
  cluster_var_name <- NULL
  M <- NULL
  M1 <- NULL
  M2 <- NULL
  if (!is.null(clusters)) {
    # Use terms() to distinguish ~a+b (two additive terms) from ~a:b (one
    # interaction term).  all.vars() strips operators, making them ambiguous.
    cl_terms <- attr(stats::terms(clusters), "term.labels")
    # Reject interaction terms (`:` or `*`) — user should pre-compute the
    # interaction variable for one-way clustering, or use `+` for two-way.
    has_interaction <- any(grepl(":", cl_terms, fixed = TRUE))
    if (has_interaction)
      stop("`clusters` must use `+` for two-way clustering (e.g. ~a + b), ",
           "not `:` or `*`. For one-way clustering on an interaction, ",
           "create the variable first (e.g. interaction(a, b)).", call. = FALSE)
    cluster_var_names <- cl_terms
    n_clvars <- length(cluster_var_names)
    if (n_clvars < 1L || n_clvars > 2L)
      stop("`clusters` must reference one or two variables.", call. = FALSE)

    # Align with model frame (respects subset + na.action).
    mf_rows <- match(rownames(parsed$model_frame), rownames(data))

    if (n_clvars == 1L) {
      # --- One-way clustering ---
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
      # --- Two-way clustering (Cameron-Gelbach-Miller) ---
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
      M <- min(M1, M2)  # Stata convention: effective M = min(M1, M2)
      cluster_vec <- list(cv1, cv2)
    }
  }

  # --- 4. Dispatch ---
  fit <- if (method %in% c("liml", "kclass")) {
    .fit_kclass(parsed, method = method, kclass = kclass, fuller = fuller,
                small = small, dofminus = dofminus, sdofminus = sdofminus)
  } else if (parsed$is_iv) {
    .fit_2sls(parsed, small = small, dofminus = dofminus,
              sdofminus = sdofminus)
  } else {
    .fit_ols(parsed, small = small, dofminus = dofminus,
             sdofminus = sdofminus)
  }

  # --- 5. VCV ---
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

  # For HC/CL VCV with weights: transform X_hat and residuals by sqrt(w)
  # at the call site so VCV functions remain weight-agnostic.
  # scores = sqrt(w)*X_hat * sqrt(w)*e = w*X_hat*e (correct weighted sandwich)
  if (!is.null(cluster_vec) || vcov %in% c("HC0", "HC1")) {
    X_hat_vcov <- if (parsed$is_iv) fit$X_hat else parsed$X
    resid_vcov <- fit$residuals
    if (!is.null(w_raw)) {
      sqrt_w <- sqrt(parsed$weights)  # normalized weights
      X_hat_vcov <- sqrt_w * X_hat_vcov
      resid_vcov <- sqrt_w * resid_vcov
    }
  }

  if (!is.null(cluster_vec)) {
    fit$vcov <- .compute_cl_vcov(bread_vcov, X_hat_vcov, resid_vcov,
                                  cluster_vec, parsed$N, parsed$K, M, small,
                                  dofminus = dofminus, sdofminus = sdofminus)
    fit$df.residual <- as.integer(M - 1L)
  } else if (vcov %in% c("HC0", "HC1")) {
    fit$vcov <- .compute_hc_vcov(bread_vcov, X_hat_vcov, resid_vcov,
                                  parsed$N, parsed$K, vcov,
                                  small = small, dofminus = dofminus,
                                  sdofminus = sdofminus)
  }

  # --- 5b. Diagnostics ---
  effective_vcov_type <- if (!is.null(cluster_vec)) "CL" else vcov

  diagnostics <- list()
  first_stage <- NULL
  if (parsed$is_iv) {

    # Overidentification test (D1)
    diagnostics$overid <- .compute_overid_test(
      Z = parsed$Z, X = parsed$X, y = parsed$y,
      residuals = fit$residuals, rss = fit$rss,
      weights = parsed$weights, cluster_vec = cluster_vec,
      vcov_type = effective_vcov_type, is_iv = parsed$is_iv,
      N = parsed$N, K = parsed$K, L = parsed$L,
      overid_df = parsed$overid_df, dofminus = dofminus
    )

    # Identification tests (D2)
    id_tests <- .compute_id_tests(
      X = parsed$X, Z = parsed$Z, y = parsed$y,
      residuals = fit$residuals, weights = parsed$weights,
      cluster_vec = cluster_vec, vcov_type = effective_vcov_type,
      N = parsed$N, K = parsed$K, L = parsed$L,
      K1 = parsed$K1, L1 = parsed$L1, M = M,
      endo_names = parsed$endo_names,
      excluded_names = parsed$excluded_names,
      has_intercept = parsed$has_intercept,
      dofminus = dofminus, sdofminus = sdofminus
    )
    diagnostics$underid        <- id_tests$underid
    diagnostics$weak_id        <- id_tests$weak_id
    diagnostics$weak_id_robust <- id_tests$weak_id_robust

    # Stock-Yogo critical values (D3)
    diagnostics$weak_id_sy <- .stock_yogo_lookup(parsed$K1, parsed$L1)

    # First-stage diagnostics (E1)
    first_stage <- .compute_first_stage(
      X = parsed$X, Z = parsed$Z,
      weights = parsed$weights, cluster_vec = cluster_vec,
      vcov_type = effective_vcov_type,
      endo_names = parsed$endo_names,
      excluded_names = parsed$excluded_names,
      N = parsed$N, K = parsed$K, L = parsed$L,
      K1 = parsed$K1, L1 = parsed$L1, M = M,
      bread_2sls = fit$bread,
      dofminus = dofminus, sdofminus = sdofminus
    )

    # Anderson-Rubin test (E3)
    diagnostics$anderson_rubin <- .compute_anderson_rubin(
      Z = parsed$Z, X = parsed$X, y = parsed$y,
      weights = parsed$weights, cluster_vec = cluster_vec,
      vcov_type = effective_vcov_type,
      N = parsed$N, K = parsed$K, L = parsed$L,
      K1 = parsed$K1, L1 = parsed$L1, M = M,
      endo_names = parsed$endo_names,
      excluded_names = parsed$excluded_names,
      dofminus = dofminus, sdofminus = sdofminus
    )

    # Stock-Wright S statistic (J2)
    diagnostics$stock_wright <- .compute_stock_wright(
      Z = parsed$Z, X = parsed$X, y = parsed$y,
      weights = parsed$weights, cluster_vec = cluster_vec,
      vcov_type = effective_vcov_type,
      N = parsed$N, K1 = parsed$K1, L1 = parsed$L1,
      endo_names = parsed$endo_names, dofminus = dofminus
    )

    # Endogeneity test / C-statistic (E4)
    if (!is.null(endog)) {
      bad <- setdiff(endog, parsed$endo_names)
      if (length(bad) > 0L) {
        stop("`endog` contains variables not in the endogenous list: ",
             paste0("'", bad, "'", collapse = ", "), ".", call. = FALSE)
      }
    }
    diagnostics$endogeneity <- .compute_endogeneity_test(
      Z = parsed$Z, X = parsed$X, y = parsed$y,
      residuals = fit$residuals, rss = fit$rss,
      weights = parsed$weights, cluster_vec = cluster_vec,
      vcov_type = effective_vcov_type,
      N = parsed$N, K = parsed$K, L = parsed$L,
      K1 = parsed$K1, endo_names = parsed$endo_names,
      endog_vars = endog, dofminus = dofminus
    )

    # Orthogonality test / instrument-subset C-statistic (J1)
    if (!is.null(orthog) && length(orthog) > 0L) {
      # Validate against actual Z column names (not term labels, which can
      # diverge for factor variables). Exclude intercept and endogenous
      # regressor columns — only instrument columns are testable.
      endo_cols <- if (parsed$K1 > 0L) parsed$endo_names else character(0L)
      valid_orthog <- setdiff(colnames(parsed$Z),
                              c("(Intercept)", endo_cols))
      bad <- setdiff(orthog, valid_orthog)
      if (length(bad) > 0L) {
        stop("`orthog` contains variables not in the instrument list: ",
             paste0("'", bad, "'", collapse = ", "),
             ". Must be excluded or exogenous instruments (not endogenous ",
             "regressors or the intercept).", call. = FALSE)
      }
      diagnostics$orthog <- .compute_orthog_test(
        Z = parsed$Z, X = parsed$X, y = parsed$y,
        residuals = fit$residuals, rss = fit$rss,
        weights = parsed$weights, cluster_vec = cluster_vec,
        vcov_type = effective_vcov_type,
        N = parsed$N, K = parsed$K, L = parsed$L,
        orthog_vars = orthog, dofminus = dofminus
      )
    }
  }
  if (length(diagnostics) == 0L) diagnostics <- NULL

  # --- 5b2. Reduced-form regression ---
  reduced_form_result <- NULL
  if (parsed$is_iv && reduced_form != "none") {
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
      endo_names     = parsed$endo_names,
      excluded_names = parsed$excluded_names,
      depvar_name    = rf_depvar,
      dofminus       = dofminus,
      sdofminus      = sdofminus
    )
  }

  # --- 5c. Model F-test ---
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

  # --- 6. Assemble return object ---
  # Determine effective method for the return object
  est_method <- if (method %in% c("liml", "kclass")) {
    method
  } else if (parsed$is_iv) {
    "2sls"
  } else {
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
    reduced_form  = reduced_form_result,
    call          = cl,
    formula       = parsed$formula,
    terms         = parsed$terms,
    nobs          = parsed$N,
    vcov_type     = if (!is.null(cluster_vec)) "CL" else vcov,
    small         = small,
    dofminus      = dofminus,
    sdofminus     = sdofminus,
    cluster_var   = cluster_var_name,
    n_clusters    = M,
    n_clusters1   = M1,
    n_clusters2   = M2,
    na.action     = parsed$na.action,
    weights       = w_raw,
    endogenous    = parsed$endo_names,
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
    model         = if (model) parsed$model_frame else NULL,
    x             = if (x) list(X = parsed$X, Z = parsed$Z) else NULL,
    y             = if (y) parsed$y else NULL
  )
}
