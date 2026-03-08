#' Phillips Curve Macroeconomic Dataset
#'
#' U.S. macroeconomic time series (1948--1996) with inflation and unemployment
#' rates plus lags and first differences. This is the same dataset used in
#' Stata's \code{ivreg2} help file for HAC and AC examples, originally from
#' Wooldridge (2019).
#'
#' @format A data frame with 49 observations and 11 variables:
#' \describe{
#'   \item{year}{Calendar year (1948--1996).}
#'   \item{inf}{Inflation rate (percent).}
#'   \item{unem}{Unemployment rate (percent).}
#'   \item{cinf}{Change in inflation (\code{inf - inf_1}).}
#'   \item{cunem}{Change in unemployment (\code{unem - unem_1}).}
#'   \item{unem_1}{Lagged unemployment (one year).}
#'   \item{inf_1}{Lagged inflation (one year).}
#'   \item{unem_2}{Lagged unemployment (two years).}
#'   \item{inf_2}{Lagged inflation (two years).}
#'   \item{cinf_1}{Lagged change in inflation.}
#'   \item{cunem_1}{Lagged change in unemployment.}
#' }
#'
#' @source
#' Wooldridge, J.M. (2019). \emph{Introductory Econometrics: A Modern
#' Approach}, 7th ed. Cengage Learning.
#'
#' @examples
#' data(phillips)
#' # OLS Phillips curve with HAC standard errors
#' fit <- ivreg2(cinf ~ unem, data = phillips,
#'               kernel = "bartlett", bw = 3, tvar = "year")
#' summary(fit)
"phillips"
