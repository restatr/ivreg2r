#' Card (1995) College Proximity Dataset
#'
#' Cross-sectional data from the National Longitudinal Survey of Young Men
#' (1966--1981) used by Card (1995) to estimate the return to schooling using
#' college proximity as an instrument for education.
#'
#' @format A data frame with 3,010 observations and 33 variables:
#' \describe{
#'   \item{nearc2}{Grew up near a 2-year college (binary).}
#'   \item{nearc4}{Grew up near a 4-year college (binary).}
#'   \item{educ}{Years of education.}
#'   \item{age}{Age in years (1976).}
#'   \item{fatheduc}{Father's years of education.}
#'   \item{motheduc}{Mother's years of education.}
#'   \item{weight}{NLS sampling weight.}
#'   \item{momdad14}{Lived with both parents at age 14 (binary).}
#'   \item{sinmom14}{Lived with single mother at age 14 (binary).}
#'   \item{step14}{Lived with stepparent at age 14 (binary).}
#'   \item{reg661}{Region dummy: New England (1966).}
#'   \item{reg662}{Region dummy: Middle Atlantic (1966).}
#'   \item{reg663}{Region dummy: East North Central (1966).}
#'   \item{reg664}{Region dummy: West North Central (1966).}
#'   \item{reg665}{Region dummy: South Atlantic (1966).}
#'   \item{reg666}{Region dummy: East South Central (1966).}
#'   \item{reg667}{Region dummy: West South Central (1966).}
#'   \item{reg668}{Region dummy: Mountain (1966).}
#'   \item{reg669}{Region dummy: Pacific (1966).}
#'   \item{south66}{Lived in the South in 1966 (binary).}
#'   \item{black}{Black (binary).}
#'   \item{smsa}{Lives in SMSA (binary, 1976).}
#'   \item{south}{Lives in the South (binary, 1976).}
#'   \item{smsa66}{Lived in SMSA in 1966 (binary).}
#'   \item{wage}{Hourly wage (cents, 1976).}
#'   \item{enroll}{Enrolled in school in 1976 (binary).}
#'   \item{KWW}{Knowledge of the World of Work test score.}
#'   \item{IQ}{IQ score.}
#'   \item{married}{Married (binary, 1976).}
#'   \item{libcrd14}{Had a library card at age 14 (binary).}
#'   \item{exper}{Years of labor market experience (\code{age - educ - 6}).}
#'   \item{lwage}{Log hourly wage.}
#'   \item{expersq}{Experience squared (\code{exper^2}).}
#' }
#'
#' @source
#' Card, D. (1995). "Using Geographic Variation in College Proximity to
#' Estimate the Return to Schooling." In L.N. Christofides, E.K. Grant, and
#' R. Swidinsky (Eds.), *Aspects of Labour Market Behaviour: Essays in Honour
#' of John Vanderkamp*. University of Toronto Press.
#'
#' @examples
#' data(card)
#' # IV regression: instrument education with college proximity
#' fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
#'               data = card)
#' summary(fit)
"card"
