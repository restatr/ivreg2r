# Build script: Wooldridge wagepan panel dataset
#
# Source: Wooldridge, J.M. (2010). Econometric Analysis of Cross Section and
# Panel Data, 2nd ed. MIT Press.
#
# Panel data from the National Longitudinal Survey of Youth (NLSY), 1980-1987.
# 4,360 observations (545 persons x 8 years). This extract matches the version
# distributed with Stata's bcuse package ("wagepan"), cached at
# ../validation/data/wagepan.dta (cached 2026-07-04), trimmed to essential
# variables for vignette use. Building from the cached .dta (rather than a
# CSV round-trip) preserves the float32 bit patterns of the Stata
# distribution.

# Path is relative to the package root (pkg/). Run with setwd("pkg/") or
# from an RStudio project rooted at pkg/.
library(haven)

d <- read_dta("../validation/data/wagepan.dta")

# Trim to essential variables for vignette use
wagepan <- data.frame(
  nr       = as.integer(d$nr),
  year     = as.integer(d$year),
  lwage    = as.numeric(d$lwage),
  educ     = as.integer(d$educ),
  black    = as.integer(d$black),
  hisp     = as.integer(d$hisp),
  exper    = as.integer(d$exper),
  expersq  = as.integer(d$expersq),
  married  = as.integer(d$married),
  union    = as.integer(d$union),
  hours    = as.integer(d$hours)
)

usethis::use_data(wagepan, overwrite = TRUE)
