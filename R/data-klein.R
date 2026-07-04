#' Klein (1950) Model I Macroeconomic Data
#'
#' Annual U.S. national-accounts aggregates for 1920--1941, from Klein's
#' Model I of the U.S. economy: consumption, profits, wages, investment,
#' capital stock, government spending, and taxes. Klein, L.R. (1950).
#' \emph{Economic Fluctuations in the United States, 1921--1941}. Cowles
#' Commission Monograph No. 11. Wiley.
#'
#' \strong{Two columns both encode "year" -- do not confuse them.}
#' \code{yr} is the calendar year (1920--1941) and is the time variable to
#' pass as \code{tvar} when using the time-series operators \code{l()} /
#' \code{d()}. \code{year} is calendar year minus 1931 (a linear trend
#' ranging from -11 to 10) and is used as an \emph{instrument} (a trend
#' term) in the help-file examples below -- it is not a second copy of the
#' time index.
#'
#' The upstream Stata dataset (\code{webuse klein}) also ships two
#' precomputed one-period lags, \code{profits1} (= \code{L.profits}) and
#' \code{totinc1} (= \code{L.totinc}). Those columns are deliberately
#' dropped here: this package's \code{l()} time-series operator computes
#' the same lags directly from \code{profits} and \code{totinc} given
#' \code{tvar = "yr"}, so shipping precomputed lag columns would be
#' redundant and could drift out of sync with \code{l()}. \code{capital1}
#' (the lagged capital stock) is kept because it is a primitive in this
#' dataset -- there is no contemporaneous \code{capital} column to lag it
#' from -- and it appears directly in the instrument lists below.
#'
#' @format A data frame with 22 observations and 12 variables:
#' \describe{
#'   \item{yr}{Calendar year (1920--1941). The time variable: use as
#'     \code{tvar} with \code{l()}/\code{d()}.}
#'   \item{consump}{Consumption.}
#'   \item{profits}{Private profits.}
#'   \item{wagepriv}{Private wage bill.}
#'   \item{invest}{Investment.}
#'   \item{capital1}{Lagged value of capital stock (a primitive; no
#'     contemporaneous capital column exists to derive it from).}
#'   \item{totinc}{Total income/demand.}
#'   \item{wagegovt}{Government wage bill.}
#'   \item{govt}{Government spending.}
#'   \item{taxnetx}{Indirect business taxes plus net exports.}
#'   \item{wagetot}{Total U.S. wage bill.}
#'   \item{year}{Calendar year minus 1931 (a linear trend, -11 to 10),
#'     used as an instrument -- not the time variable (see \code{yr}).}
#' }
#'
#' @source
#' Klein, L.R. (1950). \emph{Economic Fluctuations in the United States,
#' 1921--1941}. Cowles Commission Monograph No. 11. New York: Wiley.
#'
#' Distributed via Stata's \code{webuse klein}.
#'
#' @examples
#' data(klein)
#'
#' # ivreg2 help file line 1462: LIML, consumption on lagged profits, with
#' # profits and wagetot treated as endogenous
#' fit <- ivreg2(
#'   consump ~ l(profits, 1) | profits + wagetot |
#'     govt + taxnetx + year + wagegovt + capital1 + l(totinc, 1),
#'   data = klein, tvar = "yr", method = "liml"
#' )
#' summary(fit)
"klein"
