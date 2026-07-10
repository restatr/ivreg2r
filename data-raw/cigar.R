# Build script: Baltagi-Levin / Baltagi-Griffin-Xiong cigarette demand panel
#
# State-level cigarette demand panel: 46 U.S. states x 30 years (1963-1992),
# balanced. Baltagi, B. H. and Levin, D. (1992). "Cigarette taxation: raising
# revenues and reducing consumption." Structural Change and Economic
# Dynamics, 3(2), 321-335. See also Baltagi, B. H., Griffin, J. M. and
# Xiong, W. (2000). "To pool or not to pool: homogeneous versus
# heterogeneous estimators applied to cigarette demand." Review of
# Economics and Statistics, 82(1), 117-126.
#
# Source: plm::Cigar (cached at validation/data/cigar.dta), 1380 x 9.
# All nine columns kept as distributed -- no constructed columns here (the
# log transforms used in fixtures/vignette are computed on the fly), so the
# no-constructed-columns rule does not apply.

# Path is relative to the package root (pkg/). Run with setwd("pkg/") or
# from an RStudio project rooted at pkg/.
library(haven)

# The cached .dta is git-ignored (local cache convention for
# validation/data/); unlike the webuse/bcuse datasets, its upstream is the
# R plm package, so regenerate the cache from plm::Cigar when absent and
# re-verify bit-identity either way, as this dataset's provenance direction
# is reversed.
data(Cigar, package = "plm")
cache <- "../validation/data/cigar.dta"
if (!file.exists(cache)) {
  write_dta(Cigar, cache)
}
d <- read_dta(cache)

# tolerance = 0 makes this bit-exact for the numeric payloads (not just
# "close enough"); check.attributes = FALSE strips only haven's
# format/type-storage attributes picked up by the .dta round trip, which are
# not part of the data itself (verified necessary and sufficient: dropping
# it makes attribute differences fail the check, and no other diffs remain
# once it is set).
stopifnot(isTRUE(all.equal(as.data.frame(Cigar), as.data.frame(d),
                           check.attributes = FALSE, tolerance = 0)))

cigar <- data.frame(
  state = as.integer(d$state),
  year  = as.integer(d$year),
  price = as.numeric(d$price),
  pop   = as.numeric(d$pop),
  pop16 = as.numeric(d$pop16),
  cpi   = as.numeric(d$cpi),
  ndi   = as.numeric(d$ndi),
  sales = as.numeric(d$sales),
  pimin = as.numeric(d$pimin)
)

usethis::use_data(cigar, overwrite = TRUE)
