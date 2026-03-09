# Build script: Stock and Watson (2003) U.S. macroeconomic dataset
#
# Quarterly U.S. macroeconomic time series used in Baum, Schaffer &
# Stillman (2007, Stata Journal 7(4):465-506) to illustrate HAC, GMM,
# and CUE estimation. Original source: Stock and Watson (2003),
# "Introduction to Econometrics." Downloaded from
# http://fmwww.bc.edu/ec-p/data/stockwatson/macrodat.dta
#
# The raw dataset has 168 quarterly observations (1959Q1-2000Q4) with
# variables: UR, CPI, FFIR, TBILL, TBON, ER, GDP, date.
#
# BSS 2007 (p. 474) constructs two derived variables:
#   inf  = 100 * log(CPI / L4.CPI)     — annualized CPI inflation
#   ggdp = 100 * log(GDP / L4.GDP)     — annualized real GDP growth
#
# Their examples (pp. 475-480) use:
#   Dependent:    D.inf   (first difference of inflation)
#   Endogenous:   UR      (unemployment rate)
#   Instruments:  L2.ggdp, L.TBILL, L.ER, L.TBON
#
# We pre-compute inf, ggdp, and all needed lags/differences so users
# can run the BSS examples directly without Stata-style time-series
# operators.

library(haven)

d <- read_dta("../validation/data/macrodat.dta")

# Convert Stata quarterly date to year + quarter
# Stata quarterly: 0 = 1960q1, so year = 1960 + floor(date/4),
# quarter = (date %% 4) + 1
d$year <- 1960L + as.integer(floor(d$date / 4))
d$quarter <- as.integer((d$date %% 4) + 1L)

# Sort by time (should already be sorted, but be explicit)
d <- d[order(d$date), ]

# Derived variables from BSS 2007 p. 474
d$inf <- 100 * log(d$CPI / dplyr::lag(d$CPI, 4))
d$ggdp <- 100 * log(d$GDP / dplyr::lag(d$GDP, 4))

# First difference of inflation (D.inf in Stata)
d$dinf <- d$inf - dplyr::lag(d$inf, 1)

# Lags needed for the BSS instrument set:
#   L2.ggdp, L.TBILL, L.ER, L.TBON
d$ggdp_2 <- dplyr::lag(d$ggdp, 2)   # L2.ggdp
d$TBILL_1 <- dplyr::lag(d$TBILL, 1)  # L.TBILL
d$ER_1 <- dplyr::lag(d$ER, 1)        # L.ER
d$TBON_1 <- dplyr::lag(d$TBON, 1)    # L.TBON

stockwatson <- data.frame(
  year    = d$year,
  quarter = d$quarter,
  date    = as.numeric(d$date),
  UR      = as.numeric(d$UR),
  CPI     = as.numeric(d$CPI),
  FFIR    = as.numeric(d$FFIR),
  TBILL   = as.numeric(d$TBILL),
  TBON    = as.numeric(d$TBON),
  ER      = as.numeric(d$ER),
  GDP     = as.numeric(d$GDP),
  inf     = as.numeric(d$inf),
  ggdp    = as.numeric(d$ggdp),
  dinf    = as.numeric(d$dinf),
  ggdp_2  = as.numeric(d$ggdp_2),
  TBILL_1 = as.numeric(d$TBILL_1),
  ER_1    = as.numeric(d$ER_1),
  TBON_1  = as.numeric(d$TBON_1)
)

usethis::use_data(stockwatson, overwrite = TRUE)
