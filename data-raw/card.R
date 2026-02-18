# Build script: Card (1995) dataset
#
# Source: Card, D. (1995). "Using Geographic Variation in College Proximity to
# Estimate the Return to Schooling." In L.N. Christofides, E.K. Grant, and
# R. Swidinsky (Eds.), Aspects of Labour Market Behaviour: Essays in Honour
# of John Vanderkamp. University of Toronto Press.
#
# Data originally from the NLS Young Men Cohort (1966-1981). This extract
# matches the version distributed with Stata's bcuse package ("card").

# Path is relative to the package root (pkg/). Run with setwd("pkg/") or
# from an RStudio project rooted at pkg/.
card <- read.csv("tests/stata-benchmarks/fixtures/card_data.csv")

# Drop the row identifier used by Stata
card$id <- NULL

usethis::use_data(card, overwrite = TRUE)
