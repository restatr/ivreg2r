#' Arellano--Bond (1991) UK Employment Panel
#'
#' Employment, wages, and capital stock for UK companies: an unbalanced
#' annual panel over 1976--1984 (140 companies, 1,031 firm-year
#' observations). This is the companion dataset to Arellano, M. and Bond,
#' S. (1991), "Some Tests of Specification for Panel Data: Monte Carlo
#' Evidence and an Application to Employment Equations," \emph{Review of
#' Economic Studies}, 58(2), 277--297, shipped as distributed (levels plus
#' the log transforms and year dummies used throughout the Stata
#' literature).
#'
#' \code{n}, \code{w}, \code{k}, and \code{ys} are natural logs of
#' \code{emp}, \code{wage}, \code{cap}, and \code{indoutpt} respectively.
#' \code{yr1980}--\code{yr1984} are year-dummy indicators (1 in the named
#' year, 0 otherwise).
#'
#' @format A data frame with 1,031 observations and 16 variables:
#' \describe{
#'   \item{ind}{Industry code.}
#'   \item{year}{Calendar year (1976--1984). The time variable: use as
#'     \code{tvar}.}
#'   \item{emp}{Employment.}
#'   \item{wage}{Real wage.}
#'   \item{cap}{Gross capital stock.}
#'   \item{indoutpt}{Industry output.}
#'   \item{n}{Log employment, \code{log(emp)}.}
#'   \item{w}{Log real wage, \code{log(wage)}.}
#'   \item{k}{Log gross capital stock, \code{log(cap)}.}
#'   \item{ys}{Log industry output, \code{log(indoutpt)}.}
#'   \item{yr1980}{Year dummy: 1 if \code{year == 1980}.}
#'   \item{yr1981}{Year dummy: 1 if \code{year == 1981}.}
#'   \item{yr1982}{Year dummy: 1 if \code{year == 1982}.}
#'   \item{yr1983}{Year dummy: 1 if \code{year == 1983}.}
#'   \item{yr1984}{Year dummy: 1 if \code{year == 1984}.}
#'   \item{id}{Firm identifier (140 distinct companies). Use as \code{ivar}.}
#' }
#'
#' @source
#' Arellano, M. and Bond, S. (1991). Some tests of specification for panel
#' data: Monte Carlo evidence and an application to employment equations.
#' \emph{Review of Economic Studies}, 58(2), 277--297.
#'
#' Downloaded from \url{http://fmwww.bc.edu/ec-p/data/macro/abdata.dta}.
#'
#' @examples
#' data(abdata)
#'
#' # ivreg2 help file line 1558: one-way cluster on firm id
#' fit <- ivreg2(n ~ w + k, data = abdata, clusters = ~id)
#' summary(fit)
"abdata"
