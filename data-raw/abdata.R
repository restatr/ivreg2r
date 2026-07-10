# Build script: Arellano-Bond (1991) UK employment panel
#
# Employment, wages, and capital for UK companies, unbalanced annual panel
# 1976-1984 (140 companies, 1,031 observations). Companion data to
# Arellano, M. and Bond, S. (1991). "Some Tests of Specification for Panel
# Data: Monte Carlo Evidence and an Application to Employment Equations."
# Review of Economic Studies, 58(2), 277-297.
#
# Source: http://fmwww.bc.edu/ec-p/data/macro/abdata.dta (cached at
# validation/data/abdata.dta), 1031 x 16.
# All 16 columns kept: n/w/k/ys are logs of emp/wage/cap/indoutpt, and
# yr1980-yr1984 are year dummies -- constructed but not lags/differences,
# so decision D8 does not apply, and the help-file examples reference them
# directly.

# Path is relative to the package root (pkg/). Run with setwd("pkg/") or
# from an RStudio project rooted at pkg/.
library(haven)

d <- read_dta("../validation/data/abdata.dta")

abdata <- data.frame(
  ind      = as.integer(d$ind),
  year     = as.integer(d$year),
  emp      = as.numeric(d$emp),
  wage     = as.numeric(d$wage),
  cap      = as.numeric(d$cap),
  indoutpt = as.numeric(d$indoutpt),
  n        = as.numeric(d$n),
  w        = as.numeric(d$w),
  k        = as.numeric(d$k),
  ys       = as.numeric(d$ys),
  yr1980   = as.integer(d$yr1980),
  yr1981   = as.integer(d$yr1981),
  yr1982   = as.integer(d$yr1982),
  yr1983   = as.integer(d$yr1983),
  yr1984   = as.integer(d$yr1984),
  id       = as.integer(d$id)
)

usethis::use_data(abdata, overwrite = TRUE)
