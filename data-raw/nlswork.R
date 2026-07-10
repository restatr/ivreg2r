# Build script: NLS Young Women panel (Stata [XT] manual extract)
#
# National Longitudinal Survey of Young Women, 14-24 years old in 1968,
# interview years 1968-1988; the extract distributed with Stata's [XT]
# manual (webuse nlswork). 28,534 observations, 21 variables.
#
# Source: webuse nlswork (cached at validation/data/nlswork.dta).
# `race` ships as haven_labelled (1 = White, 2 = Black, 3 = Other) and is
# coerced to plain integer here -- codes documented in the .Rd, no factor
# conversion (faithful to the upstream numeric coding).
# All columns kept. Integral columns coerced to integer (verified: no
# fractional values); ttl_exp, tenure, ln_wage stay numeric. NAs are
# common in this dataset -- as.integer()/as.numeric() preserve them.
# compress = "xz": this is by far the largest bundled dataset.

# Path is relative to the package root (pkg/). Run with setwd("pkg/") or
# from an RStudio project rooted at pkg/.
library(haven)

d <- read_dta("../validation/data/nlswork.dta")

nlswork <- data.frame(
  idcode   = as.integer(d$idcode),
  year     = as.integer(d$year),
  birth_yr = as.integer(d$birth_yr),
  age      = as.integer(d$age),
  race     = as.integer(d$race),
  msp      = as.integer(d$msp),
  nev_mar  = as.integer(d$nev_mar),
  grade    = as.integer(d$grade),
  collgrad = as.integer(d$collgrad),
  not_smsa = as.integer(d$not_smsa),
  c_city   = as.integer(d$c_city),
  south    = as.integer(d$south),
  ind_code = as.integer(d$ind_code),
  occ_code = as.integer(d$occ_code),
  union    = as.integer(d$union),
  wks_ue   = as.integer(d$wks_ue),
  ttl_exp  = as.numeric(d$ttl_exp),
  tenure   = as.numeric(d$tenure),
  hours    = as.integer(d$hours),
  wks_work = as.integer(d$wks_work),
  ln_wage  = as.numeric(d$ln_wage)
)

usethis::use_data(nlswork, overwrite = TRUE, compress = "xz")
