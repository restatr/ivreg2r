# ============================================================================
# Tests: bundled datasets (klein, grunfeld, abdata, nlswork) -- coherence
#
# Checks dims/column names for each of the four datasets bundled alongside
# planning/25-data-provenance.md, and cross-checks klein/abdata against the
# full-precision fixture data exports from the ts-operator fixture generator
# (tsop_klein_data.csv / tsop_ab_data.csv, ticket F4) to prove the bundled
# .rda and the fixture CSVs carry identical numbers.
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

test_that("bundled klein matches the tsop_klein_data.csv fixture export", {
  data(klein)
  expect_matches_fixture_csv(klein, "tsop_klein_data.csv", "klein")
})

test_that("bundled abdata matches the tsop_ab_data.csv fixture export", {
  data(abdata)
  expect_matches_fixture_csv(abdata, "tsop_ab_data.csv", "abdata")
})
