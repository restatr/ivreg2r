# ============================================================================
# Tests: bundled datasets (klein, grunfeld, abdata, nlswork, cigar) -- coherence
#
# Checks dims/column names for each of the datasets bundled alongside
# planning/25-data-provenance.md, and cross-checks klein/abdata against the
# full-precision fixture data exports from the ts-operator fixture generator
# (tsop_klein_data.csv / tsop_ab_data.csv, ticket F4) to prove the bundled
# .rda and the fixture CSVs carry identical numbers. cigar is cross-checked
# directly against its cached source-of-record .dta file (validation/data/
# cigar.dta), since it has no fixture CSV export.
# ============================================================================

expect_matches_fixture_csv <- function(data, csv_name, label) {
  path <- fixture_path(csv_name)
  skip_if(!file.exists(path))
  fx <- read.csv(path, na.strings = c("", "."))
  shared <- intersect(names(data), names(fx))
  expect_true(length(shared) > 0)
  expect_equal(nrow(data), nrow(fx))
  for (col in shared) {
    expect_equal(data[[col]], fx[[col]], tolerance = 1e-12,
                 info = paste(label, "column", col))
  }
}

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

test_that("bundled klein matches the tsop_klein_data.csv fixture export", {
  data(klein)
  expect_matches_fixture_csv(klein, "tsop_klein_data.csv", "klein")
})

test_that("bundled abdata matches the tsop_ab_data.csv fixture export", {
  data(abdata)
  expect_matches_fixture_csv(abdata, "tsop_ab_data.csv", "abdata")
})

test_that("bundled cigar matches the cached validation/data/cigar.dta source", {
  skip_if_not_installed("haven")
  dta_path <- file.path(testthat::test_path(), "..", "..", "..", "validation",
                         "data", "cigar.dta")
  skip_if(!file.exists(dta_path), "cigar.dta source-of-record not found")
  data(cigar)
  src <- as.data.frame(haven::read_dta(dta_path))
  expect_equal(nrow(cigar), nrow(src))
  for (col in names(cigar)) {
    expect_equal(cigar[[col]], as.numeric(src[[col]]), tolerance = 1e-12,
                 info = paste("cigar column", col))
  }
})
