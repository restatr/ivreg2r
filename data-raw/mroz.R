# Build script: Mroz (1987) female labor supply dataset
#
# Source: Mroz, T.A. (1987). "The Sensitivity of an Empirical Model of
# Married Women's Hours of Work to Economic and Statistical Assumptions."
# Econometrica, 55(4), 765-799.
#
# Data from the Panel Study of Income Dynamics (PSID), 1975. 753 married
# women, of whom 428 were working (inlf == 1) with observed wages.
# This extract matches the version distributed with Stata's bcuse package
# ("mroz"), cached at ../validation/data/mroz.dta (cached 2026-07-04).
# Building from the cached .dta (rather than a CSV round-trip) preserves
# the float32 bit patterns of the Stata distribution.

# Path is relative to the package root (pkg/). Run with setwd("pkg/") or
# from an RStudio project rooted at pkg/.
library(haven)

d <- read_dta("../validation/data/mroz.dta")

mroz <- data.frame(
  inlf      = as.integer(d$inlf),
  hours     = as.integer(d$hours),
  kidslt6   = as.integer(d$kidslt6),
  kidsge6   = as.integer(d$kidsge6),
  age       = as.integer(d$age),
  educ      = as.integer(d$educ),
  wage      = as.numeric(d$wage),
  repwage   = as.numeric(d$repwage),
  hushrs    = as.integer(d$hushrs),
  husage    = as.integer(d$husage),
  huseduc   = as.integer(d$huseduc),
  huswage   = as.numeric(d$huswage),
  faminc    = as.integer(d$faminc),
  mtr       = as.numeric(d$mtr),
  motheduc  = as.integer(d$motheduc),
  fatheduc  = as.integer(d$fatheduc),
  unem      = as.numeric(d$unem),
  city      = as.integer(d$city),
  exper     = as.integer(d$exper),
  nwifeinc  = as.numeric(d$nwifeinc),
  lwage     = as.numeric(d$lwage),
  expersq   = as.integer(d$expersq)
)

usethis::use_data(mroz, overwrite = TRUE)
