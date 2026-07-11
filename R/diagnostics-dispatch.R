# --------------------------------------------------------------------------
# ivreg2 pipeline: diagnostic tests, reduced form, first stage, model F
# Called only by ivreg2(); see R/ivreg2.R for the orchestrator.
# --------------------------------------------------------------------------

#' Warn that an IV-only request is ignored on a model with no endogenous
#' regressors.
#'
#' Stata silently drops these requests; this package's convention is that an
#' explicitly requested computation never disappears silently. One wording
#' per model shape: `model_desc` distinguishes the plain OLS form from the
#' IV form with zero endogenous regressors, and `tail` completes the
#' sentence for the specific request.
#' @noRd
.warn_iv_request_ignored <- function(arg, model_desc, tail) {
  warning("`", arg, "` ignored: ", model_desc, ", so ", tail, call. = FALSE)
}

.OLS_MODEL_DESC <- paste0("this is an OLS model (no endogenous regressors ",
                          "or excluded instruments)")
.K10_MODEL_DESC <- "the model has no endogenous regressors"
.B0_MODEL_DESC <- "b0 evaluation suppresses identification diagnostics"
.NOID_MODEL_DESC <- "identification tests are suppressed (`noid = TRUE`)"

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

    # Kiefer/kernel precedence (Stata's ranktest never receives the kernel
    # when kiefer is set — see the comment inside the id-test block below).
    # Hoisted outside the `!noid` gate immediately below, even though only
    # that gate's id tests consume it, because the disclosure-note trigger
    # further down also reads id_kernel and runs unconditionally
    # whenever K1 > 0 -- if noid = TRUE skipped the block that used to define
    # id_kernel, that trigger would hit an undefined-variable error.
    id_kernel <- if (isTRUE(kiefer)) NULL else kernel

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
      # against Stata ivreg2 4.1.12 (plain: Anderson/CD match rk output to
      # ~1e-11; pw+kiefer: KP rk LM/Wald differ from iid and match Stata).
      id_vcov_type <- if (isTRUE(kiefer)) {
        if (effective_vcov_type == "HAC") "robust" else "iid"
      } else {
        effective_vcov_type
      }
      # id_kernel is now hoisted above the `!noid` gate (see comment there).
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
        j_full = if (!is.null(diagnostics$overid)) diagnostics$overid$stat,
        sw = sw_flag, ivar_vec = ivar_vec
      )
    }

    # Redundancy test (P1) — computation only if !noid (validation already
    # done above). K1 = 0: nothing to test redundancy against — Stata
    # silently ignores redundant() here (gated at ivreg2.ado:1741); we warn
    # instead, following the package's explicit-request-ignored convention.
    if (!is.null(redundant) && length(redundant) > 0L && parsed$K1 == 0L) {
      .warn_iv_request_ignored("redundant", .K10_MODEL_DESC,
                               "there is nothing to test instrument redundancy against.")
    }
    if (noid && !is.null(redundant) && length(redundant) > 0L &&
        parsed$K1 > 0L) {
      .warn_iv_request_ignored("redundant", .NOID_MODEL_DESC,
                               "the redundancy test is not computed.")
    }
    if (!noid && !is.null(redundant) && length(redundant) > 0L &&
        parsed$K1 > 0L) {
      redundant_cols <- .expand_terms_to_colnames(
        redundant, parsed$excluded_names, parsed$excluded_colnames,
        parsed$excluded_assign
      )
      # Stata's ranktest never receives centering (ivreg2.ado:1753-1764), so
      # the redundancy test's KP omega is always uncentered, matching the
      # underid/weak-id treatment at .compute_id_tests() above -- `center`
      # is deliberately not forwarded here and the function default (FALSE)
      # governs.
      diagnostics$redundancy <- .compute_redundancy_test(
        X = parsed$X, Z = parsed$Z,
        weights = parsed$weights, cluster_vec = cluster_vec,
        vcov_type = effective_vcov_type,
        N = parsed$N, K1 = parsed$K1,
        endo_colnames = parsed$endo_colnames,
        excluded_colnames = parsed$excluded_colnames,
        redundant_vars = redundant_cols, dofminus = dofminus,
        weight_type = weight_type,
        kernel = kernel, bw = bw, time_index = time_index
      )
    }

    # --- Stata-quirk disclosure notes -------------------------------------
    # Five sites where an explicit user option is silently not honored inside
    # an automatically computed diagnostic get a "note-as-data" disclosure: a
    # note string attached to the affected slot, surfaced by diagnostics() as
    # a `note` column and printed by summary() as a footnote. Numbers are
    # unchanged; only the note field is added. Each note text is identical
    # across the slots it touches so summary() prints it once.
    #
    # Maintenance contract: these five triggers are re-derived here, post-hoc,
    # from options and flags that were already consumed further up the
    # pipeline, rather than being computed inside the diagnostic functions
    # themselves. That means ANY new diagnostic path that silently drops a
    # user option must (1) add its trigger to this block and (2) add a test
    # to test-diag-notes.R pinning both the fires-on-trigger and
    # absent-without-it cases -- the tests in that file are the tripwires
    # that catch a future diagnostic silently reintroducing an undisclosed
    # drop.
    non_bartlett_kernel <- !is.null(id_kernel) && id_kernel != "Bartlett"
    # Deliberately narrowed to robust-family VCE only, not every psd fit.
    # Under vcov = "iid", psd corrects the stored $S,
    # which is already positive semidefinite by construction (the iid moment
    # covariance is a scalar multiple of (Z'Z), an inner-product Gram matrix),
    # so the correction is provably inert on this path. Noting an
    # option that cannot possibly change the identification/redundancy
    # statistics would disclose the non-application of an inapplicable
    # option, not a real Stata-parity quirk, so the trigger fires only where
    # psd could actually matter: robust, HAC, AC, or cluster VCE.
    psd_on_robust <- !is.null(psd) &&
      (sw_flag ||
         effective_vcov_type %in% c("robust", "HAC", "AC", "CL"))
    kernel_note <- "Computed with the Bartlett kernel regardless of the kernel= option, matching Stata's ranktest."
    kiefer_note <- "Computed without the Kiefer VCE structure, which Stata's ranktest omits from the identification tests."
    psd_note <- "Computed without the psd correction, which Stata's ranktest never receives."
    center_endog_note <- "Computed without moment centering; Stata does not forward center to the endogeneity test."
    cue_endog_note <- "C-statistic from a recursive re-estimation that does not use CUE, matching Stata's ivreg2."
    orthog_recursive_note <- "C-statistic from a recursive re-estimation that does not use the CUE or LIML estimator, matching Stata's ivreg2."

    diagnostics$underid <- .attach_diag_note(
      diagnostics$underid, non_bartlett_kernel, kernel_note)
    diagnostics$underid <- .attach_diag_note(
      diagnostics$underid, isTRUE(kiefer), kiefer_note)
    diagnostics$underid <- .attach_diag_note(
      diagnostics$underid, psd_on_robust, psd_note)

    diagnostics$weak_id <- .attach_diag_note(
      diagnostics$weak_id, isTRUE(kiefer), kiefer_note)

    diagnostics$weak_id_robust <- .attach_diag_note(
      diagnostics$weak_id_robust, non_bartlett_kernel, kernel_note)
    diagnostics$weak_id_robust <- .attach_diag_note(
      diagnostics$weak_id_robust, isTRUE(kiefer), kiefer_note)
    diagnostics$weak_id_robust <- .attach_diag_note(
      diagnostics$weak_id_robust, psd_on_robust, psd_note)

    diagnostics$redundancy <- .attach_diag_note(
      diagnostics$redundancy, non_bartlett_kernel, kernel_note)
    diagnostics$redundancy <- .attach_diag_note(
      diagnostics$redundancy, psd_on_robust, psd_note)

    diagnostics$endogeneity <- .attach_diag_note(
      diagnostics$endogeneity, isTRUE(center), center_endog_note)
    diagnostics$endogeneity <- .attach_diag_note(
      diagnostics$endogeneity, method == "cue", cue_endog_note)

    diagnostics$orthog <- .attach_diag_note(
      diagnostics$orthog, method %in% c("cue", "liml"), orthog_recursive_note)

    } else {
      # b0 evaluates the CUE objective at a fixed parameter vector with no
      # optimization, so the entire identification-diagnostics family above
      # (underid, weak-id, first-stage, AR, Stock-Wright, endogeneity,
      # orthog, redundancy) is not computed (Stata line 3819). reduced_form
      # is NOT suppressed by b0 (verified empirically) so it is deliberately
      # absent from this list — see the reduced-form gate further below.
      # Explicit requests among these would otherwise vanish silently;
      # warn instead, following the package's explicit-request-ignored
      # convention (the same one used for the K1 = 0 and OLS cases).
      if (!is.null(endog) && length(endog) > 0L) {
        .warn_iv_request_ignored("endog", .B0_MODEL_DESC,
                                 "the endogeneity test is not computed.")
      }
      if (!is.null(orthog) && length(orthog) > 0L) {
        .warn_iv_request_ignored("orthog", .B0_MODEL_DESC,
                                 "the orthogonality test is not computed.")
      }
      if (!is.null(redundant) && length(redundant) > 0L) {
        .warn_iv_request_ignored("redundant", .B0_MODEL_DESC,
                                 "the redundancy test is not computed.")
      }
      if (isTRUE(first_stage_flag)) {
        .warn_iv_request_ignored("first_stage", .B0_MODEL_DESC,
                                 "the first-stage regressions are not stored.")
      }
    }  # end of if (is.null(b0)) — identification diagnostics block
  } else {
    # OLS path (one-part formula): the IV-only requests have no target here.
    # Warn rather than silently drop them (the explicit-request-ignored
    # convention; the K1 = 0 IV form gets the same treatment above and at
    # the reduced-form/first-stage gates below).
    if (length(endog) > 0L) {
      .warn_iv_request_ignored("endog", .OLS_MODEL_DESC,
                               "there is nothing to test.")
    }
    if (length(orthog) > 0L) {
      .warn_iv_request_ignored("orthog", .OLS_MODEL_DESC,
                               "there is nothing to test.")
    }
    if (length(redundant) > 0L) {
      .warn_iv_request_ignored("redundant", .OLS_MODEL_DESC,
                               "there is nothing to test.")
    }
    if (reduced_form != "none") {
      .warn_iv_request_ignored("reduced_form", .OLS_MODEL_DESC,
                               "there is no reduced form to store.")
    }
    if (isTRUE(first_stage_flag)) {
      .warn_iv_request_ignored("first_stage", .OLS_MODEL_DESC,
                               "there are no first-stage regressions to store.")
    }
  }
  if (length(diagnostics) == 0L) diagnostics <- NULL

  # Reduced-form regression. K1 = 0: no reduced form is computed or stored,
  # matching Stata — the entire RF/first-stage estimation block is gated on
  # `endo1_ct > 0` (ivreg2.ado:1256), and saverf on a `(=z)` model stores
  # nothing (e(rfeq) empty; verified against Stata ivreg2 4.1.12). Note the
  # y-on-Z
  # reduced form WOULD differ from the main y-on-X model (Z carries the
  # surplus instruments), but Stata never computes it, so there is no ground
  # truth to validate an extension against. Where Stata is silent about the
  # dropped request, we warn (the explicit-request-ignored convention).
  reduced_form_result <- NULL
  if (parsed$is_iv && parsed$K1 == 0L && reduced_form != "none") {
    .warn_iv_request_ignored("reduced_form", .K10_MODEL_DESC,
                             "there is no reduced form to store.")
  }
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

  # Build extractable first-stage model objects
  first_stage_models <- NULL
  if (isTRUE(first_stage_flag) && parsed$is_iv && parsed$K1 == 0L) {
    .warn_iv_request_ignored("first_stage", .K10_MODEL_DESC,
                             "there are no first-stage regressions to store.")
  }
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
