# ============================================================================
# Tests: bundled datasets (klein, grunfeld, abdata, nlswork, cigar) -- coherence
#
# Checks dims/column names for each of the datasets bundled alongside planning/25-data-provenance.md. The klein/abdata cross-checks against the ts-operator fixture CSVs (tsop_klein_data.csv / tsop_ab_data.csv, ticket F4) were retired 2026-07-06 along with those CSVs (data-CSV retirement): coherence is now structural, because the fixture generator loads the same ../validation/data/klein.dta and abdata.dta files from which pkg/data-raw builds the bundled .rda, following the M-14 precedent of doing source-of-record verification at build time in data-raw rather than in live tests. cigar's bit-identity cross-check against its cached source-of-record .dta file (validation/data/cigar.dta) runs at build time in data-raw/cigar.R for the same reason: a test gated on that git-ignored cache would SKIP on any checkout without it (public repo, fresh clones), which violates the 0-skips rule.
# ============================================================================

test_that("klein has the documented dimensions and column names", {
  data(klein)
  expect_equal(dim(klein), c(22L, 12L))
  expect_named(klein, c("yr", "consump", "profits", "wagepriv", "invest",
                         "capital1", "totinc", "wagegovt", "govt", "taxnetx",
                         "wagetot", "year"))
})

test_that("grunfeld has the documented dimensions and column names", {
  data(grunfeld)
  expect_equal(dim(grunfeld), c(200L, 6L))
  expect_named(grunfeld, c("company", "year", "invest", "mvalue", "kstock",
                           "time"))
  expect_equal(length(unique(grunfeld$company)), 10L)
  expect_equal(length(unique(grunfeld$year)), 20L)
})

test_that("abdata has the documented dimensions and column names", {
  data(abdata)
  expect_equal(dim(abdata), c(1031L, 16L))
  expect_named(abdata, c("ind", "year", "emp", "wage", "cap", "indoutpt",
                         "n", "w", "k", "ys", "yr1980", "yr1981", "yr1982",
                         "yr1983", "yr1984", "id"))
})

test_that("nlswork has the documented dimensions and column names", {
  data(nlswork)
  expect_equal(dim(nlswork), c(28534L, 21L))
  expect_named(nlswork, c("idcode", "year", "birth_yr", "age", "race", "msp",
                          "nev_mar", "grade", "collgrad", "not_smsa",
                          "c_city", "south", "ind_code", "occ_code", "union",
                          "wks_ue", "ttl_exp", "tenure", "hours", "wks_work",
                          "ln_wage"))
  expect_setequal(unique(nlswork$race), c(1L, 2L, 3L))
})

test_that("cigar has the documented dimensions and column names", {
  data(cigar)
  expect_equal(dim(cigar), c(1380L, 9L))
  expect_named(cigar, c("state", "year", "price", "pop", "pop16", "cpi",
                        "ndi", "sales", "pimin"))
  expect_equal(length(unique(cigar$state)), 46L)
  expect_equal(length(unique(cigar$year)), 30L)
})
