# Build script: Griliches (1976) wages dataset
#
# Young men's wages from the National Longitudinal Survey (NLS), used in
# Griliches (1976, J.Pol.Ec.) and extensively in Hayashi (2000, Ch. 3).
# Also used in Baum, Schaffer & Stillman (2007, Stata Journal) to
# illustrate weak identification, Anderson-Rubin, Stock-Wright, and
# redundancy tests (pp. 492-494).
#
# Source: http://fmwww.bc.edu/ec-p/data/hayashi/griliches76.dta
# (also available via Stata: use http://www.stata-press.com/data/imeus/griliches)
#
# 758 observations, 20 variables. Two cross-sections (year 66-73 and
# year 80). The BSS 2007 examples use the pooled sample with year
# dummies via factor(year).

library(haven)

d <- read_dta("../validation/data/griliches76.dta")

griliches <- data.frame(
  rns      = as.numeric(d$rns),
  rns80    = as.numeric(d$rns80),
  mrt      = as.numeric(d$mrt),
  mrt80    = as.numeric(d$mrt80),
  smsa     = as.numeric(d$smsa),
  smsa80   = as.numeric(d$smsa80),
  med      = as.numeric(d$med),
  iq       = as.numeric(d$iq),
  kww      = as.numeric(d$kww),
  year     = as.integer(d$year),
  age      = as.numeric(d$age),
  age80    = as.numeric(d$age80),
  s        = as.numeric(d$s),
  s80      = as.numeric(d$s80),
  expr     = as.numeric(d$expr),
  expr80   = as.numeric(d$expr80),
  tenure   = as.numeric(d$tenure),
  tenure80 = as.numeric(d$tenure80),
  lw       = as.numeric(d$lw),
  lw80     = as.numeric(d$lw80)
)

usethis::use_data(griliches, overwrite = TRUE)
