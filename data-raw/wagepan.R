# Build script: Wooldridge wagepan panel dataset
#
# Source: Wooldridge, J.M. (2010). Econometric Analysis of Cross Section and
# Panel Data, 2nd ed. MIT Press.
#
# Panel data from the National Longitudinal Survey of Youth (NLSY), 1980-1987.
# 4,360 observations (545 persons x 8 years). This extract matches the version
# distributed with Stata's bcuse package ("wagepan"), trimmed to essential
# variables for vignette use.

# Path is relative to the package root (pkg/). Run with setwd("pkg/") or
# from an RStudio project rooted at pkg/.
wagepan <- read.csv("tests/stata-benchmarks/fixtures/wagepan_data.csv")

# Trim to essential variables for vignette use
wagepan <- wagepan[, c("nr", "year", "lwage", "educ", "black", "hisp",
                        "exper", "expersq", "married", "union", "hours")]

usethis::use_data(wagepan, overwrite = TRUE)
