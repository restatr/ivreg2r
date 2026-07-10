# Build script: Grunfeld (1958) investment panel
#
# Corporate investment panel: 10 firms x 20 years (1935-1954), balanced.
# Grunfeld, Y. (1958). The Determinants of Corporate Investment. PhD
# dissertation, University of Chicago. See also Kleiber, C. and Zeileis, A.
# (2010), "The Grunfeld Data at 50," German Economic Review 11(4), 404-417,
# for provenance and the many circulating variants of this dataset.
#
# Source: webuse grunfeld (cached at validation/data/grunfeld.dta), 200 x 6.
# All six columns kept, including `time` (a 1-20 trend index -- constructed
# but not a lag/difference, so decision D8 does not apply).

# Path is relative to the package root (pkg/). Run with setwd("pkg/") or
# from an RStudio project rooted at pkg/.
library(haven)

d <- read_dta("../validation/data/grunfeld.dta")

grunfeld <- data.frame(
  company = as.integer(d$company),
  year    = as.integer(d$year),
  invest  = as.numeric(d$invest),
  mvalue  = as.numeric(d$mvalue),
  kstock  = as.numeric(d$kstock),
  time    = as.integer(d$time)
)

usethis::use_data(grunfeld, overwrite = TRUE)
