# Build script: Card (1995) dataset
#
# Source: Card, D. (1995). "Using Geographic Variation in College Proximity to
# Estimate the Return to Schooling." In L.N. Christofides, E.K. Grant, and
# R. Swidinsky (Eds.), Aspects of Labour Market Behaviour: Essays in Honour
# of John Vanderkamp. University of Toronto Press.
#
# Data originally from the NLS Young Men Cohort (1966-1981). This extract
# matches the version distributed with Stata's bcuse package ("card"),
# cached at ../validation/data/card.dta (cached 2026-07-04). Building from
# the cached .dta (rather than a CSV round-trip) preserves the float32 bit
# patterns of the Stata distribution, per planning/25 ground rule 1.

# Path is relative to the package root (pkg/). Run with setwd("pkg/") or
# from an RStudio project rooted at pkg/.
library(haven)

d <- read_dta("../validation/data/card.dta")

# Drop the row identifier used by Stata (id)
card <- data.frame(
  nearc2   = as.integer(d$nearc2),
  nearc4   = as.integer(d$nearc4),
  educ     = as.integer(d$educ),
  age      = as.integer(d$age),
  fatheduc = as.integer(d$fatheduc),
  motheduc = as.integer(d$motheduc),
  weight   = as.integer(d$weight),
  momdad14 = as.integer(d$momdad14),
  sinmom14 = as.integer(d$sinmom14),
  step14   = as.integer(d$step14),
  reg661   = as.integer(d$reg661),
  reg662   = as.integer(d$reg662),
  reg663   = as.integer(d$reg663),
  reg664   = as.integer(d$reg664),
  reg665   = as.integer(d$reg665),
  reg666   = as.integer(d$reg666),
  reg667   = as.integer(d$reg667),
  reg668   = as.integer(d$reg668),
  reg669   = as.integer(d$reg669),
  south66  = as.integer(d$south66),
  black    = as.integer(d$black),
  smsa     = as.integer(d$smsa),
  south    = as.integer(d$south),
  smsa66   = as.integer(d$smsa66),
  wage     = as.integer(d$wage),
  enroll   = as.integer(d$enroll),
  KWW      = as.integer(d$KWW),
  IQ       = as.integer(d$IQ),
  married  = as.integer(d$married),
  libcrd14 = as.integer(d$libcrd14),
  exper    = as.integer(d$exper),
  lwage    = as.numeric(d$lwage),
  expersq  = as.integer(d$expersq)
)

usethis::use_data(card, overwrite = TRUE)
