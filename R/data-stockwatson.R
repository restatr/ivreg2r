#' Stock and Watson U.S. Macroeconomic Dataset
#'
#' Quarterly U.S. macroeconomic time series (1959Q1--2000Q4) used in
#' Baum, Schaffer & Stillman (2007) to illustrate HAC standard errors,
#' two-step GMM, and CUE estimation. Original source: Stock and Watson (2003).
#'
#' The derived variables \code{inf}, \code{ggdp}, \code{dinf}, and the lagged
#' instrument columns are pre-computed following BSS 2007 (p. 474) so that
#' users can replicate their examples directly without time-series operators.
#'
#' @format A data frame with 168 observations and 17 variables:
#' \describe{
#'   \item{year}{Calendar year.}
#'   \item{quarter}{Quarter (1--4).}
#'   \item{date}{Stata quarterly date (0 = 1960Q1).}
#'   \item{UR}{Unemployment rate (percent).}
#'   \item{CPI}{Consumer Price Index.}
#'   \item{FFIR}{Federal funds interest rate.}
#'   \item{TBILL}{Treasury bill rate.}
#'   \item{TBON}{Treasury bond rate.}
#'   \item{ER}{Trade-weighted exchange rate.}
#'   \item{GDP}{Nominal GDP (billions of dollars).}
#'   \item{inf}{Annualized CPI inflation: \code{100 * log(CPI / L4.CPI)}.
#'     \code{NA} for the first 4 observations.}
#'   \item{ggdp}{Annualized GDP growth: \code{100 * log(GDP / L4.GDP)}.
#'     \code{NA} for the first 4 observations.}
#'   \item{dinf}{First difference of inflation (\code{inf - L.inf}).
#'     \code{NA} for the first 5 observations.}
#'   \item{ggdp_2}{Second lag of GDP growth (\code{L2.ggdp}).}
#'   \item{TBILL_1}{First lag of Treasury bill rate (\code{L.TBILL}).}
#'   \item{ER_1}{First lag of exchange rate (\code{L.ER}).}
#'   \item{TBON_1}{First lag of Treasury bond rate (\code{L.TBON}).}
#' }
#'
#' @source
#' Stock, J.H. and Watson, M.W. (2003). \emph{Introduction to Econometrics}.
#' Addison-Wesley.
#'
#' Downloaded from \url{http://fmwww.bc.edu/ec-p/data/stockwatson/macrodat.dta}.
#'
#' @references
#' Baum, C.F., Schaffer, M.E. and Stillman, S. (2007). Enhanced routines for
#' instrumental variables/generalized method of moments estimation and testing.
#' \emph{The Stata Journal}, 7(4), 465--506.
#'
#' @examples
#' data(stockwatson)
#'
#' # BSS 2007, p. 475: 2SLS with IID errors
#' fit_iv <- ivreg2(dinf ~ 1 | UR | ggdp_2 + TBILL_1 + ER_1 + TBON_1,
#'                  data = stockwatson)
#' summary(fit_iv)
#'
#' \donttest{
#' # BSS 2007, p. 476: two-step GMM with HAC (Bartlett, bw=5)
#' fit_gmm <- ivreg2(dinf ~ 1 | UR | ggdp_2 + TBILL_1 + ER_1 + TBON_1,
#'                   data = stockwatson, method = "gmm2s",
#'                   vcov = "robust", kernel = "bartlett", bw = 5,
#'                   tvar = "date")
#' summary(fit_gmm)
#' }
"stockwatson"
