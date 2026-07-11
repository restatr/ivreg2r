#' Cigarette Demand Panel (Baltagi--Levin / Baltagi--Griffin--Xiong)
#'
#' Annual per-capita cigarette sales for 46 U.S. states over 1963--1992
#' (a balanced panel, 1,380 state-year observations), along with the price
#' per pack, population, consumer price index, per-capita disposable
#' income, and the minimum price per pack among neighboring states used to
#' instrument for price endogeneity. Baltagi, B. H. and Levin, D. (1992),
#' "Cigarette taxation: raising revenues and reducing consumption,"
#' \emph{Structural Change and Economic Dynamics}, 3(2), 321--335. See
#' also Baltagi, B. H., Griffin, J. M. and Xiong, W. (2000), "To pool or
#' not to pool: homogeneous versus heterogeneous estimators applied to
#' cigarette demand," \emph{Review of Economics and Statistics}, 82(1),
#' 117--126.
#'
#' \code{price} and \code{ndi} are nominal; deflate by \code{cpi} for real
#' terms.
#'
#' @format A data frame with 1,380 observations and 9 variables (46 states
#'   times 30 years, 1963--1992):
#' \describe{
#'   \item{state}{State identifier code (46 distinct U.S. states). Use as
#'     \code{ivar}.}
#'   \item{year}{Year, coded 63--92 (i.e., 1963--1992). The time variable:
#'     use as \code{tvar}.}
#'   \item{price}{Price per pack of cigarettes (nominal).}
#'   \item{pop}{Population.}
#'   \item{pop16}{Population above the age of 16.}
#'   \item{cpi}{Consumer price index (1983 = 100).}
#'   \item{ndi}{Per-capita nominal disposable income.}
#'   \item{sales}{Cigarette sales, in packs per capita.}
#'   \item{pimin}{Minimum price per pack of cigarettes in adjoining
#'     states.}
#' }
#'
#' @source
#' Baltagi, B. H. and Levin, D. (1992). Cigarette taxation: raising
#' revenues and reducing consumption. \emph{Structural Change and
#' Economic Dynamics}, 3(2), 321--335.
#'
#' Baltagi, B. H., Griffin, J. M. and Xiong, W. (2000). To pool or not to
#' pool: homogeneous versus heterogeneous estimators applied to cigarette
#' demand. \emph{Review of Economics and Statistics}, 82(1), 117--126.
#'
#' Distributed via R package \code{plm}, \code{data(Cigar)}.
#'
#' Redistribution basis: \code{system.file("COPYRIGHTS", package = "ivreg2r")}.
#'
#' @examples
#' data(cigar)
#'
#' # Price-endogeneity IV specification in real (CPI-deflated) terms, with
#' # the neighboring-state minimum price as the excluded instrument, as in
#' # the GFIC empirical example.
#' cigar_real <- transform(
#'   cigar,
#'   lsales  = log(sales),
#'   lrprice = log(price / cpi),
#'   lrndi   = log(ndi / cpi),
#'   lrpimin = log(pimin / cpi)
#' )
#' fit <- ivreg2(lsales ~ lrndi | lrprice | lrpimin, data = cigar_real)
#' summary(fit)
#' @family ivreg2r datasets
"cigar"
