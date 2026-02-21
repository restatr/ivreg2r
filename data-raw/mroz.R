# Build script: Mroz (1987) female labor supply dataset
#
# Source: Mroz, T.A. (1987). "The Sensitivity of an Empirical Model of
# Married Women's Hours of Work to Economic and Statistical Assumptions."
# Econometrica, 55(4), 765-799.
#
# Data from the Panel Study of Income Dynamics (PSID), 1975. 753 married
# women, of whom 428 were working (inlf == 1) with observed wages.
# This extract matches the version distributed with Stata's bcuse package
# ("mroz").

# Path is relative to the package root (pkg/). Run with setwd("pkg/") or
# from an RStudio project rooted at pkg/.
mroz <- read.csv("tests/stata-benchmarks/fixtures/mroz_data.csv")

usethis::use_data(mroz, overwrite = TRUE)
