#' Grunfeld (1958) Corporate Investment Panel
#'
#' Corporate investment, market value, and capital stock for 10 large U.S.
#' firms observed annually over 1935--1954: a balanced panel of 10 firms
#' times 20 years (200 observations). Grunfeld, Y. (1958). \emph{The
#' Determinants of Corporate Investment}. PhD dissertation, University of
#' Chicago. See Kleiber, C. and Zeileis, A. (2010) for the definitive
#' account of the dataset's provenance and its many circulating variants.
#'
#' @format A data frame with 200 observations and 6 variables (10 firms
#'   times 20 years, 1935--1954):
#' \describe{
#'   \item{company}{Firm identifier (1--10). Use as \code{ivar}.}
#'   \item{year}{Calendar year (1935--1954). The time variable: use as
#'     \code{tvar}.}
#'   \item{invest}{Gross investment, current year.}
#'   \item{mvalue}{Market value of the firm, prior year.}
#'   \item{kstock}{Capital stock, prior year.}
#'   \item{time}{Time index (1--20), a linear trend.}
#' }
#'
#' @source
#' Grunfeld, Y. (1958). \emph{The Determinants of Corporate Investment}.
#' PhD dissertation, University of Chicago.
#'
#' Kleiber, C. and Zeileis, A. (2010). The Grunfeld data at 50.
#' \emph{German Economic Review}, 11(4), 404--417.
#'
#' Distributed via Stata's \code{webuse grunfeld}.
#'
#' Redistribution basis: \code{system.file("COPYRIGHTS", package = "ivreg2r")}.
#'
#' @examples
#' data(grunfeld)
#'
#' # ivreg2 help file line 1595: Driscoll-Kraay standard errors (bw = 2),
#' # small-sample correction
#' fit <- ivreg2(
#'   invest ~ mvalue + kstock, data = grunfeld,
#'   dkraay = 2, small = TRUE, tvar = "year", ivar = "company"
#' )
#' summary(fit)
#' @family ivreg2r datasets
"grunfeld"
