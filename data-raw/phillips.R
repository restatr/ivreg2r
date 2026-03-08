# Build script: Phillips (1958) curve dataset
#
# U.S. macroeconomic time series, 1948-1996 (49 observations).
# Variables: inflation, unemployment, and their lags/changes.
# This is the same dataset used in Stata's ivreg2 help file for
# HAC and AC examples. Source: Wooldridge (2019), originally from
# Economic Report of the President.

# Path is relative to the package root (pkg/).
library(haven)
d <- read_dta("../validation/data/phillips.dta")

phillips <- data.frame(
  year   = as.integer(d$year),
  inf    = as.numeric(d$inf),
  unem   = as.numeric(d$unem),
  cinf   = as.numeric(d$cinf),
  cunem  = as.numeric(d$cunem),
  unem_1 = as.numeric(d$unem_1),
  inf_1  = as.numeric(d$inf_1),
  unem_2 = as.numeric(d$unem_2),
  inf_2  = as.numeric(d$inf_2),
  cinf_1 = as.numeric(d$cinf_1),
  cunem_1 = as.numeric(d$cunem_1)
)

usethis::use_data(phillips, overwrite = TRUE)
