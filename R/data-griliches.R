#' Griliches (1976) Young Men's Wages Dataset
#'
#' Wages and characteristics of young men from the National Longitudinal
#' Survey (NLS), originally analyzed in Griliches (1976). This dataset is
#' a standard IV example used in Hayashi (2000, Ch. 3) and Baum, Schaffer
#' & Stillman (2007, pp. 492--494) to illustrate weak identification tests,
#' Anderson--Rubin inference, Stock--Wright S statistics, and instrument
#' redundancy.
#'
#' Two cross-sections are pooled (survey years 66--73 and 80). The Baum,
#' Schaffer & Stillman (2007) examples use year dummies via
#' \code{factor(year)}.
#'
#' @format A data frame with 758 observations and 20 variables:
#' \describe{
#'   \item{rns}{Residency in the South (1 = yes).}
#'   \item{rns80}{Residency in the South in 1980.}
#'   \item{mrt}{Marital status (1 = married).}
#'   \item{mrt80}{Marital status in 1980.}
#'   \item{smsa}{Resides in metropolitan area (1 = urban).}
#'   \item{smsa80}{Resides in metropolitan area in 1980.}
#'   \item{med}{Mother's education (years).}
#'   \item{iq}{IQ score.}
#'   \item{kww}{Score on Knowledge of the World of Work test.}
#'   \item{year}{Survey year (66--73 or 80).}
#'   \item{age}{Age at survey.}
#'   \item{age80}{Age in 1980.}
#'   \item{s}{Completed years of schooling.}
#'   \item{s80}{Completed years of schooling in 1980.}
#'   \item{expr}{Work experience (years).}
#'   \item{expr80}{Work experience in 1980 (years).}
#'   \item{tenure}{Job tenure (years).}
#'   \item{tenure80}{Job tenure in 1980 (years).}
#'   \item{lw}{Log wage.}
#'   \item{lw80}{Log wage in 1980.}
#' }
#'
#' @source
#' Griliches, Z. (1976). Wages of very young men. \emph{Journal of Political
#' Economy}, 84(4), S69--S85.
#'
#' Downloaded from \url{http://fmwww.bc.edu/ec-p/data/hayashi/griliches76.dta}.
#'
#' @references
#' Hayashi, F. (2000). \emph{Econometrics}. Princeton University Press.
#' (Chapter 3, p. 255.)
#'
#' Baum, C.F., Schaffer, M.E. and Stillman, S. (2007). Enhanced routines for
#' instrumental variables/generalized method of moments estimation and testing.
#' \emph{The Stata Journal}, 7(4), 465--506. (pp. 492--494.)
#'
#' @examples
#' data(griliches)
#'
#' # Baum, Schaffer & Stillman (2007), p. 493: IV with robust SEs, weak-ID
#' # diagnostics, redundancy test
#' fit <- ivreg2(lw ~ s + expr + tenure + rns + smsa + factor(year) |
#'                 iq | age + mrt,
#'               data = griliches, vcov = "robust",
#'               redundant = "mrt")
#' summary(fit)
"griliches"
