# --------------------------------------------------------------------------
# Stock-Yogo critical values lookup
# --------------------------------------------------------------------------

#' Look up Stock-Yogo (2005) critical values for weak identification tests
#'
#' Returns a data.frame of critical values for the Cragg-Donald / KP Wald F
#' statistic, covering different threshold tables depending on the estimation
#' method. Tables are available for K1 (endogenous regressors) <= 3 (IV bias)
#' or K1 <= 2 (all others), and L1 (excluded instruments) 1-100.
#'
#' Dispatch logic (matches Stata's `Disp_cdsy` in ivreg2.ado):
#' - 2SLS/OLS: IV relative bias + IV size distortion
#' - LIML (no Fuller): LIML size distortion only
#' - Fuller (LIML with fuller > 0): Fuller relative bias + Fuller maximum bias
#' - kclass: NULL (no Stock-Yogo tables)
#'
#' @param K1 Integer, number of endogenous regressors.
#' @param L1 Integer, number of excluded instruments.
#' @param method Character, estimation method: "ols", "2sls", "liml", "kclass".
#' @param fuller Numeric, Fuller parameter (0 = plain LIML).
#' @return A data.frame with columns `type`, `threshold`, `critical_value`,
#'   or NULL if no tables are available for the given K1/L1/method combination.
#' @keywords internal
#' @noRd
.stock_yogo_lookup <- function(K1, L1, method = "2sls", fuller = 0) {
  # kclass has no Stock-Yogo tables
  if (method == "kclass") return(NULL)

  rows <- list()

  .add_rows <- function(tables, type_label, max_K1) {
    if (K1 > max_K1) return()
    for (thresh in names(tables)) {
      cv <- if (L1 >= 1L && L1 <= 100L) tables[[thresh]][L1, K1] else NA_real_
      if (!is.na(cv)) {
        rows[[length(rows) + 1L]] <<- data.frame(
          type = type_label,
          threshold = thresh,
          critical_value = cv,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (method %in% c("ols", "2sls")) {
    # IV relative bias tables (K1 <= 3)
    .add_rows(list(
      "5%"  = sy_ivbias5,  "10%" = sy_ivbias10,
      "20%" = sy_ivbias20, "30%" = sy_ivbias30
    ), "IV relative bias", 3L)

    # IV size distortion tables (K1 <= 2)
    .add_rows(list(
      "10%" = sy_ivsize10, "15%" = sy_ivsize15,
      "20%" = sy_ivsize20, "25%" = sy_ivsize25
    ), "IV size", 2L)

  } else if (method == "liml" && fuller == 0) {
    # LIML size distortion tables (K1 <= 2)
    .add_rows(list(
      "10%" = sy_limlsize10, "15%" = sy_limlsize15,
      "20%" = sy_limlsize20, "25%" = sy_limlsize25
    ), "LIML size", 2L)

  } else if (method == "liml" && fuller > 0) {
    # Fuller relative bias tables (K1 <= 2)
    .add_rows(list(
      "5%"  = sy_fullrel5,  "10%" = sy_fullrel10,
      "20%" = sy_fullrel20, "30%" = sy_fullrel30
    ), "Fuller relative bias", 2L)

    # Fuller maximum bias tables (K1 <= 2)
    .add_rows(list(
      "5%"  = sy_fullmax5,  "10%" = sy_fullmax10,
      "20%" = sy_fullmax20, "30%" = sy_fullmax30
    ), "Fuller maximum bias", 2L)
  }

  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}
