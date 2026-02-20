# --------------------------------------------------------------------------
# ivreg2
# --------------------------------------------------------------------------
#' Extended Instrumental Variables Estimation
#'
#' Estimate models by OLS or two-stage least squares (2SLS) with automatic
#' diagnostic tests. Uses a three-part formula for IV:
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
                   small = FALSE,
                   dofminus = 0L, sdofminus = 0L,
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

  # --- 3b. Validate and normalize weights ---
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
    cluster_var_names <- all.vars(clusters)
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
  fit <- if (parsed$is_iv) {
    .fit_2sls(parsed, small = small, dofminus = dofminus,
              sdofminus = sdofminus)
  } else {
    .fit_ols(parsed, small = small, dofminus = dofminus,
             sdofminus = sdofminus)
  }

  # --- 5. VCV ---
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
    fit$vcov <- .compute_cl_vcov(fit$bread, X_hat_vcov, resid_vcov,
                                  cluster_vec, parsed$N, parsed$K, M, small,
                                  dofminus = dofminus, sdofminus = sdofminus)
    fit$df.residual <- as.integer(M - 1L)
  } else if (vcov %in% c("HC0", "HC1")) {
    fit$vcov <- .compute_hc_vcov(fit$bread, X_hat_vcov, resid_vcov,
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
  }
  if (length(diagnostics) == 0L) diagnostics <- NULL

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
    model         = if (model) parsed$model_frame else NULL,
    x             = if (x) list(X = parsed$X, Z = parsed$Z) else NULL,
    y             = if (y) parsed$y else NULL
  )
}
