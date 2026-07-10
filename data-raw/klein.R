# Build script: Klein (1950) Model I macro data
#
# Interwar US national-accounts aggregates (consumption, profits, wages,
# investment, capital stock, government spending, taxes), annual, 1920-1941.
# Klein, L.R. (1950). Economic Fluctuations in the United States, 1921-1941.
# Cowles Commission Monograph No. 11. Wiley.
#
# Source: webuse klein (cached at validation/data/klein.dta), 22 x 14.
#
# Deliberate deviation (decision D8): the
# upstream file carries two precomputed lag columns, profits1 (= L.profits)
# and totinc1 (= L.totinc). These are dropped -- the help-file examples
# reference L.profits/L.totinc, expressible as l(profits, 1)/l(totinc, 1)
# since ticket F4, and D8 forbids shipping constructed lag columns.
# capital1 (lagged capital stock) is kept: it is a primitive with no
# contemporaneous "capital" column to lag it from, and it appears directly
# in the instrument lists of help-file examples H72-H76. Bundled shape:
# 22 x 12.

# Path is relative to the package root (pkg/). Run with setwd("pkg/") or
# from an RStudio project rooted at pkg/.
library(haven)

d <- read_dta("../validation/data/klein.dta")

klein <- data.frame(
  yr       = as.integer(d$yr),
  consump  = as.numeric(d$consump),
  profits  = as.numeric(d$profits),
  wagepriv = as.numeric(d$wagepriv),
  invest   = as.numeric(d$invest),
  capital1 = as.numeric(d$capital1),
  totinc   = as.numeric(d$totinc),
  wagegovt = as.numeric(d$wagegovt),
  govt     = as.numeric(d$govt),
  taxnetx  = as.numeric(d$taxnetx),
  wagetot  = as.numeric(d$wagetot),
  year     = as.integer(d$year)
)

usethis::use_data(klein, overwrite = TRUE)
