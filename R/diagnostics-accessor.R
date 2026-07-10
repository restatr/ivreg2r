# ============================================================================
# diagnostics() accessor
# ============================================================================
# Flattens the IV specification tests stored on `x$diagnostics` into a single
# tidy tibble, one row per test, so they can be filtered/joined like any
# other tidyverse output instead of navigated as a nested list.


# --------------------------------------------------------------------------
# diagnostics() generic + method
# --------------------------------------------------------------------------

#' Extract IV diagnostic tests as a tibble
#'
#' Flattens the IV specification tests stored on a fitted model into a single
#' tibble, one row per test, so they can be filtered, joined, or printed like
#' any other tidyverse output instead of navigated as a nested list
#' (`x$diagnostics$overid$stat` and friends).
#'
#' @param x A fitted model object.
#' @param ... Additional arguments (ignored).
#' @return A [tibble::tibble()] with one row per diagnostic test present on
#'   the model and the following columns:
#'   \describe{
#'     \item{`test`}{Character. A stable machine-readable key, meant for
#'       `filter(test == ...)`-style selection (e.g. `"overid"`,
#'       `"weak_id_robust"`, `"endogeneity"`).}
#'     \item{`test_name`}{Character. The display name printed by
#'       `summary()`.}
#'     \item{`statistic`}{Double. The test statistic, or — for the
#'       Stock-Yogo rows — the critical value rather than a statistic.}
#'     \item{`df`}{Integer. Degrees of freedom, `NA` when not applicable.}
#'     \item{`df2`}{Integer. Denominator degrees of freedom, `NA` unless the
#'       statistic is F-distributed.}
#'     \item{`p_value`}{Double. `NA` when not defined (e.g. Stock-Yogo
#'       critical values, or an exactly identified overidentification row).}
#'     \item{`tested_vars`}{Character. The variable(s) tested; populated only
#'       for the endogeneity, orthogonality, and redundancy rows.}
#'   }
#'
#'   Rows appear only for the tests actually computed on the model: an IV
#'   fit reports underidentification, weak identification, and
#'   overidentification by default, while `endog =`, `orthog =`, and
#'   `redundant =` add the corresponding endogeneity, orthogonality, and
#'   redundancy rows.
#'   The Stock-Yogo rows (keys of the form `"sy_iv_size_10"`) report critical
#'   values, not test statistics, so their `p_value` is always `NA`. An OLS
#'   fit (single-part formula) has no diagnostics at all and returns a
#'   zero-row tibble with the columns above.
#'
#' @examples
#' data(card)
#' fit <- ivreg2(
#'   lwage ~ exper + expersq + black + south | educ | nearc4 + nearc2,
#'   data = card, vcov = "robust"
#' )
#' diagnostics(fit)
#'
#' @family ivreg2 methods
#' @seealso [first_stage()], [glance.ivreg2()]
#' @export
diagnostics <- function(x, ...) UseMethod("diagnostics")

#' @rdname diagnostics
#' @export
diagnostics.ivreg2 <- function(x, ...) {
  diag <- x$diagnostics
  rows <- list()

  if (!is.null(diag$underid)) {
    uid <- diag$underid
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "underid", test_name = uid$test_name,
      statistic = uid$stat, df = uid$df, p_value = uid$p
    )
  }

  if (!is.null(diag$weak_id)) {
    wid <- diag$weak_id
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "weak_id", test_name = wid$test_name, statistic = wid$stat
    )
  }

  if (!is.null(diag$weak_id_robust)) {
    widr <- diag$weak_id_robust
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "weak_id_robust", test_name = widr$test_name,
      statistic = widr$stat
    )
  }

  if (!is.null(diag$weak_id_sy)) {
    sy <- diag$weak_id_sy
    for (i in seq_len(nrow(sy))) {
      key <- paste0(
        "sy_", .diag_slug(sy$type[i]), "_",
        sub("%$", "", sy$threshold[i])
      )
      rows[[length(rows) + 1L]] <- .diag_row(
        test = key,
        test_name = paste0(
          "Stock-Yogo critical value (", sy$type[i], ", ",
          sy$threshold[i], ")"
        ),
        statistic = sy$critical_value[i]
      )
    }
  }

  if (!is.null(diag$overid)) {
    oid <- diag$overid
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "overid", test_name = oid$test_name,
      statistic = oid$stat, df = oid$df, p_value = oid$p
    )
  }

  if (!is.null(diag$anderson_rubin_overid)) {
    aro <- diag$anderson_rubin_overid
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "anderson_rubin_overid_lr",
      test_name = "Anderson-Rubin overidentification (LR)",
      statistic = aro$lr_stat, df = aro$df, p_value = aro$lr_p
    )
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "anderson_rubin_overid_lin",
      test_name = "Anderson-Rubin overidentification (linearized)",
      statistic = aro$lin_stat, df = aro$df, p_value = aro$lin_p
    )
  }

  if (!is.null(diag$anderson_rubin)) {
    ar <- diag$anderson_rubin
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "anderson_rubin_f", test_name = "Anderson-Rubin Wald F",
      statistic = ar$f_stat, df = ar$f_df1, df2 = ar$f_df2, p_value = ar$f_p
    )
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "anderson_rubin_chi2", test_name = "Anderson-Rubin Wald chi-sq",
      statistic = ar$chi2_stat, df = ar$chi2_df, p_value = ar$chi2_p
    )
  }

  if (!is.null(diag$stock_wright)) {
    sw <- diag$stock_wright
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "stock_wright", test_name = "Stock-Wright LM S",
      statistic = sw$stat, df = sw$df, p_value = sw$p
    )
  }

  if (!is.null(diag$endogeneity)) {
    endog <- diag$endogeneity
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "endogeneity", test_name = endog$test_name,
      statistic = endog$stat, df = endog$df, p_value = endog$p,
      tested_vars = endog$tested_vars
    )
  }

  if (!is.null(diag$orthog)) {
    orth <- diag$orthog
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "orthog", test_name = orth$test_name,
      statistic = orth$stat, df = orth$df, p_value = orth$p,
      tested_vars = orth$tested_vars
    )
  }

  if (!is.null(diag$redundancy)) {
    red <- diag$redundancy
    rows[[length(rows) + 1L]] <- .diag_row(
      test = "redundancy", test_name = red$test_name,
      statistic = red$stat, df = red$df, p_value = red$p,
      tested_vars = red$tested_vars
    )
  }

  if (length(rows) == 0L) {
    return(tibble::tibble(
      test        = character(0),
      test_name   = character(0),
      statistic   = double(0),
      df          = integer(0),
      df2         = integer(0),
      p_value     = double(0),
      tested_vars = character(0)
    ))
  }

  do.call(rbind, rows)
}


# --------------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------------

#' Build one row of the diagnostics() tibble
#' @keywords internal
#' @noRd
.diag_row <- function(test, test_name, statistic, df = NA_integer_,
                       df2 = NA_integer_, p_value = NA_real_,
                       tested_vars = NA_character_) {
  tested_vars_out <- if (length(tested_vars) == 1L && is.na(tested_vars)) {
    NA_character_
  } else {
    paste(tested_vars, collapse = ", ")
  }
  tibble::tibble(
    test        = test,
    test_name   = test_name,
    statistic   = as.double(statistic),
    df          = as.integer(df),
    df2         = as.integer(df2),
    p_value     = as.double(p_value),
    tested_vars = tested_vars_out
  )
}

#' Slugify a Stock-Yogo test-type label for use in a machine key
#'
#' Lowercases and replaces any run of non-alphanumeric characters with a
#' single underscore, e.g. `"IV size"` -> `"iv_size"`.
#' @keywords internal
#' @noRd
.diag_slug <- function(x) {
  gsub("[^a-z0-9]+", "_", tolower(x))
}
