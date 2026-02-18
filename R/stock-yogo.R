# --------------------------------------------------------------------------
# Stock-Yogo critical values lookup
# --------------------------------------------------------------------------

#' Look up Stock-Yogo (2005) critical values for weak identification tests
#'
#' Returns a data.frame of critical values for the Cragg-Donald / KP Wald F
#' statistic, covering IV relative bias and IV size distortion thresholds.
#' Tables are available for K1 (endogenous regressors) <= 3 (bias) or
#' K1 <= 2 (size), and L1 (excluded instruments) 1-100.
#'
#' @param K1 Integer, number of endogenous regressors.
#' @param L1 Integer, number of excluded instruments.
#' @return A data.frame with columns `type`, `threshold`, `critical_value`,
#'   or NULL if no tables are available for the given K1/L1 combination.
#' @keywords internal
#' @noRd
.stock_yogo_lookup <- function(K1, L1) {
  # No tables available for K1 > 3

  if (K1 > 3L) return(NULL)

  rows <- list()

  # IV relative bias tables (K1 <= 3)
  if (K1 <= 3L) {
    bias_tables <- list(
      "5%"  = sy_ivbias5,
      "10%" = sy_ivbias10,
      "20%" = sy_ivbias20,
      "30%" = sy_ivbias30
    )
    for (thresh in names(bias_tables)) {
      cv <- if (L1 >= 1L && L1 <= 100L) bias_tables[[thresh]][L1, K1] else NA_real_
      if (!is.na(cv)) {
        rows[[length(rows) + 1L]] <- data.frame(
          type = "IV relative bias",
          threshold = thresh,
          critical_value = cv,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # IV size distortion tables (K1 <= 2)
  if (K1 <= 2L) {
    size_tables <- list(
      "10%" = sy_ivsize10,
      "15%" = sy_ivsize15,
      "20%" = sy_ivsize20,
      "25%" = sy_ivsize25
    )
    for (thresh in names(size_tables)) {
      cv <- if (L1 >= 1L && L1 <= 100L) size_tables[[thresh]][L1, K1] else NA_real_
      if (!is.na(cv)) {
        rows[[length(rows) + 1L]] <- data.frame(
          type = "IV size",
          threshold = thresh,
          critical_value = cv,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}
