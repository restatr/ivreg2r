# ============================================================================
# Tests: Stock-Yogo Critical Values (Ticket D3)
# ============================================================================
#
# Unit tests for .stock_yogo_lookup() and integration tests verifying
# diagnostics$weak_id_sy is populated correctly by ivreg2().

# --- Load datasets ---
fixture_dir <- file.path(
  testthat::test_path(), "..", "stata-benchmarks", "fixtures"
)

card_path <- file.path(fixture_dir, "card_data.csv")
if (file.exists(card_path)) {
  card <- read.csv(card_path)
}

sim_multi_path <- file.path(fixture_dir, "sim_multi_endo_data.csv")
if (file.exists(sim_multi_path)) {
  sim_multi <- read.csv(sim_multi_path)
}


# ============================================================================
# Unit tests: .stock_yogo_lookup()
# ============================================================================

test_that("K1=1, L1=1: ivsize10 = 16.38, no bias rows (all NA)", {
  result <- ivreg2r:::.stock_yogo_lookup(1L, 1L)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("type", "threshold", "critical_value"))

  # Only size rows should be present (bias rows 1-2 are NA for K1=1)
  size_rows <- result[result$type == "IV size", ]
  expect_true(nrow(size_rows) > 0L)
  expect_equal(size_rows$critical_value[size_rows$threshold == "10%"], 16.38)
  expect_equal(size_rows$critical_value[size_rows$threshold == "15%"], 8.96)
  expect_equal(size_rows$critical_value[size_rows$threshold == "20%"], 6.66)
  expect_equal(size_rows$critical_value[size_rows$threshold == "25%"], 5.53)

  # Bias rows: L1=1 is NA for all bias tables (K1=1 needs L1 >= 3)
  bias_rows <- result[result$type == "IV relative bias", ]
  expect_equal(nrow(bias_rows), 0L)
})

test_that("K1=1, L1=2: ivsize10 = 19.93, ivsize15 = 11.59", {
  result <- ivreg2r:::.stock_yogo_lookup(1L, 2L)
  size_rows <- result[result$type == "IV size", ]
  expect_equal(size_rows$critical_value[size_rows$threshold == "10%"], 19.93)
  expect_equal(size_rows$critical_value[size_rows$threshold == "15%"], 11.59)
})

test_that("K1=2, L1=4: ivbias10 = 7.56, ivsize10 = 16.87", {
  result <- ivreg2r:::.stock_yogo_lookup(2L, 4L)
  bias_rows <- result[result$type == "IV relative bias", ]
  size_rows <- result[result$type == "IV size", ]
  expect_equal(bias_rows$critical_value[bias_rows$threshold == "10%"], 7.56)
  expect_equal(size_rows$critical_value[size_rows$threshold == "10%"], 16.87)
})

test_that("K1=3, L1=5: ivbias5 = 9.53, no size tables for K1=3", {
  result <- ivreg2r:::.stock_yogo_lookup(3L, 5L)
  bias_rows <- result[result$type == "IV relative bias", ]
  size_rows <- result[result$type == "IV size", ]
  expect_equal(bias_rows$critical_value[bias_rows$threshold == "5%"], 9.53)
  expect_equal(nrow(size_rows), 0L)
})

test_that("K1 > 3 returns NULL", {
  expect_null(ivreg2r:::.stock_yogo_lookup(4L, 10L))
  expect_null(ivreg2r:::.stock_yogo_lookup(5L, 20L))
})

test_that("L1 > 100 returns NULL (all values out of range)", {
  expect_null(ivreg2r:::.stock_yogo_lookup(1L, 101L))
  expect_null(ivreg2r:::.stock_yogo_lookup(2L, 200L))
})

test_that("L1 < 1 returns NULL", {
  expect_null(ivreg2r:::.stock_yogo_lookup(1L, 0L))
  expect_null(ivreg2r:::.stock_yogo_lookup(2L, -1L))
})

test_that("Return is data.frame with correct columns", {
  result <- ivreg2r:::.stock_yogo_lookup(1L, 5L)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("type", "threshold", "critical_value"))
  expect_type(result$type, "character")
  expect_type(result$threshold, "character")
  expect_type(result$critical_value, "double")
})

test_that("K1=1, L1=3: bias tables start having values", {
  result <- ivreg2r:::.stock_yogo_lookup(1L, 3L)
  bias_rows <- result[result$type == "IV relative bias", ]
  # L1=3 is the first row with non-NA for K1=1 bias tables
  expect_equal(bias_rows$critical_value[bias_rows$threshold == "5%"], 13.91)
  expect_equal(bias_rows$critical_value[bias_rows$threshold == "10%"], 9.08)
  expect_equal(bias_rows$critical_value[bias_rows$threshold == "20%"], 6.46)
  expect_equal(bias_rows$critical_value[bias_rows$threshold == "30%"], 5.39)
})


# ============================================================================
# Integration tests: diagnostics$weak_id_sy in ivreg2()
# ============================================================================

test_that("card_just_id: weak_id_sy has IV size rows, first cv = 16.38", {
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  sy <- fit$diagnostics$weak_id_sy
  expect_s3_class(sy, "data.frame")

  size_rows <- sy[sy$type == "IV size", ]
  expect_true(nrow(size_rows) > 0L)
  expect_equal(size_rows$critical_value[size_rows$threshold == "10%"], 16.38)
})

test_that("card_overid: weak_id_sy contains ivsize10 = 19.93", {
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc2 + nearc4,
                data = card)
  sy <- fit$diagnostics$weak_id_sy
  expect_s3_class(sy, "data.frame")

  size_rows <- sy[sy$type == "IV size", ]
  expect_equal(size_rows$critical_value[size_rows$threshold == "10%"], 19.93)
  expect_equal(size_rows$critical_value[size_rows$threshold == "15%"], 11.59)
})

test_that("sim_multi_endo: weak_id_sy has both bias and size rows", {
  skip_if(!file.exists(sim_multi_path), "sim_multi_endo data not found")

  # K1=2 endogenous, L1=4 excluded instruments
  fit <- ivreg2(y ~ x1 + x2 | endo1 + endo2 | z1 + z2 + z3 + z4,
                data = sim_multi)
  sy <- fit$diagnostics$weak_id_sy
  expect_s3_class(sy, "data.frame")

  bias_rows <- sy[sy$type == "IV relative bias", ]
  size_rows <- sy[sy$type == "IV size", ]
  expect_true(nrow(bias_rows) > 0L)
  expect_true(nrow(size_rows) > 0L)

  expect_equal(bias_rows$critical_value[bias_rows$threshold == "10%"], 7.56)
  expect_equal(size_rows$critical_value[size_rows$threshold == "10%"], 16.87)
})

test_that("OLS model: weak_id_sy is NULL", {
  fit <- ivreg2(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$diagnostics$weak_id_sy)
})

test_that("weak_id_sy is VCE-invariant (same tables for iid vs robust)", {
  skip_if(!file.exists(card_path), "card data not found")

  fit_iid <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card)
  fit_hc1 <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                    data = card, vcov = "HC1")

  expect_equal(fit_iid$diagnostics$weak_id_sy, fit_hc1$diagnostics$weak_id_sy)
})


# ============================================================================
# Unit tests: LIML/Fuller/kclass Stock-Yogo tables (Ticket H4)
# ============================================================================

# --- Spot-check values for new tables ---

test_that("Fuller relative bias tables: spot-check values", {
  result <- ivreg2r:::.stock_yogo_lookup(1L, 1L, method = "liml", fuller = 1)
  full_rel <- result[result$type == "Fuller relative bias", ]
  expect_equal(full_rel$critical_value[full_rel$threshold == "5%"], 24.09)
  expect_equal(full_rel$critical_value[full_rel$threshold == "10%"], 19.36)
  expect_equal(full_rel$critical_value[full_rel$threshold == "20%"], 15.64)
  expect_equal(full_rel$critical_value[full_rel$threshold == "30%"], 12.71)
})

test_that("Fuller relative bias K1=2: spot-check", {
  result <- ivreg2r:::.stock_yogo_lookup(2L, 2L, method = "liml", fuller = 1)
  full_rel <- result[result$type == "Fuller relative bias", ]
  expect_equal(full_rel$critical_value[full_rel$threshold == "5%"], 15.50)
  expect_equal(full_rel$critical_value[full_rel$threshold == "10%"], 12.55)
})

test_that("Fuller maximum bias tables: spot-check values", {
  result <- ivreg2r:::.stock_yogo_lookup(1L, 1L, method = "liml", fuller = 1)
  full_max <- result[result$type == "Fuller maximum bias", ]
  expect_equal(full_max$critical_value[full_max$threshold == "5%"], 23.81)
  expect_equal(full_max$critical_value[full_max$threshold == "10%"], 19.40)
  expect_equal(full_max$critical_value[full_max$threshold == "20%"], 15.39)
  expect_equal(full_max$critical_value[full_max$threshold == "30%"], 12.76)
})

test_that("Fuller maximum bias K1=2: spot-check", {
  result <- ivreg2r:::.stock_yogo_lookup(2L, 2L, method = "liml", fuller = 4)
  full_max <- result[result$type == "Fuller maximum bias", ]
  expect_equal(full_max$critical_value[full_max$threshold == "5%"], 14.19)
  expect_equal(full_max$critical_value[full_max$threshold == "10%"], 11.92)
})

test_that("LIML size distortion tables: spot-check values", {
  result <- ivreg2r:::.stock_yogo_lookup(1L, 1L, method = "liml", fuller = 0)
  liml_size <- result[result$type == "LIML size", ]
  expect_equal(liml_size$critical_value[liml_size$threshold == "10%"], 16.38)
  expect_equal(liml_size$critical_value[liml_size$threshold == "15%"], 8.96)
  expect_equal(liml_size$critical_value[liml_size$threshold == "20%"], 6.66)
  expect_equal(liml_size$critical_value[liml_size$threshold == "25%"], 5.53)
})

test_that("LIML size K1=2, L1=2: spot-check", {
  result <- ivreg2r:::.stock_yogo_lookup(2L, 2L, method = "liml", fuller = 0)
  liml_size <- result[result$type == "LIML size", ]
  expect_equal(liml_size$critical_value[liml_size$threshold == "10%"], 7.03)
  expect_equal(liml_size$critical_value[liml_size$threshold == "15%"], 4.58)
  expect_equal(liml_size$critical_value[liml_size$threshold == "20%"], 3.95)
  expect_equal(liml_size$critical_value[liml_size$threshold == "25%"], 3.63)
})

test_that("Last row (L1=100) spot-checks for new tables", {
  # Fuller relative bias
  result <- ivreg2r:::.stock_yogo_lookup(1L, 100L, method = "liml", fuller = 1)
  full_rel <- result[result$type == "Fuller relative bias", ]
  expect_equal(full_rel$critical_value[full_rel$threshold == "5%"], 1.55)

  # Fuller maximum bias
  full_max <- result[result$type == "Fuller maximum bias", ]
  expect_equal(full_max$critical_value[full_max$threshold == "5%"], 1.46)

  # LIML size
  result2 <- ivreg2r:::.stock_yogo_lookup(1L, 100L, method = "liml", fuller = 0)
  liml_size <- result2[result2$type == "LIML size", ]
  expect_equal(liml_size$critical_value[liml_size$threshold == "10%"], 3.83)

  # LIML size K1=2
  result3 <- ivreg2r:::.stock_yogo_lookup(2L, 100L, method = "liml", fuller = 0)
  liml_size3 <- result3[result3$type == "LIML size", ]
  expect_equal(liml_size3$critical_value[liml_size3$threshold == "10%"], 4.59)
  expect_equal(liml_size3$critical_value[liml_size3$threshold == "25%"], 2.11)
})

# --- Method dispatch ---

test_that("method='2sls' returns IV tables (default behavior)", {
  result <- ivreg2r:::.stock_yogo_lookup(1L, 5L, method = "2sls")
  types <- unique(result$type)
  expect_true("IV relative bias" %in% types)
  expect_true("IV size" %in% types)
  expect_false("LIML size" %in% types)
  expect_false("Fuller relative bias" %in% types)
})

test_that("method='liml', fuller=0 returns LIML size only", {
  result <- ivreg2r:::.stock_yogo_lookup(1L, 5L, method = "liml", fuller = 0)
  types <- unique(result$type)
  expect_equal(types, "LIML size")
})

test_that("method='liml', fuller=1 returns Fuller rel bias + max bias", {
  result <- ivreg2r:::.stock_yogo_lookup(1L, 5L, method = "liml", fuller = 1)
  types <- unique(result$type)
  expect_true("Fuller relative bias" %in% types)
  expect_true("Fuller maximum bias" %in% types)
  expect_false("LIML size" %in% types)
  expect_false("IV size" %in% types)
})

test_that("method='kclass' returns NULL", {
  expect_null(ivreg2r:::.stock_yogo_lookup(1L, 5L, method = "kclass"))
  expect_null(ivreg2r:::.stock_yogo_lookup(2L, 10L, method = "kclass"))
})

# --- K1 boundary for LIML/Fuller tables ---

test_that("K1=2 returns LIML/Fuller results, K1=3 returns NULL", {
  # LIML size: K1=2 works, K1=3 returns NULL
  result <- ivreg2r:::.stock_yogo_lookup(2L, 5L, method = "liml", fuller = 0)
  expect_s3_class(result, "data.frame")

  expect_null(ivreg2r:::.stock_yogo_lookup(3L, 5L, method = "liml", fuller = 0))

  # Fuller: K1=2 works, K1=3 returns NULL
  result2 <- ivreg2r:::.stock_yogo_lookup(2L, 5L, method = "liml", fuller = 1)
  expect_s3_class(result2, "data.frame")

  expect_null(ivreg2r:::.stock_yogo_lookup(3L, 5L, method = "liml", fuller = 1))
})

# --- Backward compatibility ---

test_that("Default args (no method) match explicit method='2sls'", {
  default <- ivreg2r:::.stock_yogo_lookup(1L, 5L)
  explicit <- ivreg2r:::.stock_yogo_lookup(1L, 5L, method = "2sls")
  expect_equal(default, explicit)
})

test_that("Fuller with any positive value dispatches to Fuller tables", {
  result1 <- ivreg2r:::.stock_yogo_lookup(1L, 5L, method = "liml", fuller = 1)
  result4 <- ivreg2r:::.stock_yogo_lookup(1L, 5L, method = "liml", fuller = 4)
  # Same tables (Fuller parameter doesn't affect which tables are used)
  expect_equal(result1, result4)
})


# ============================================================================
# Integration tests: LIML/Fuller Stock-Yogo in ivreg2() (Ticket H4)
# ============================================================================

test_that("ivreg2 method='liml' stores LIML size tables", {
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, method = "liml")
  sy <- fit$diagnostics$weak_id_sy
  expect_s3_class(sy, "data.frame")

  types <- unique(sy$type)
  expect_equal(types, "LIML size")
  expect_equal(sy$critical_value[sy$threshold == "10%"], 16.38)
})

test_that("ivreg2 fuller=1 stores Fuller bias tables", {
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, fuller = 1)
  sy <- fit$diagnostics$weak_id_sy
  expect_s3_class(sy, "data.frame")

  types <- unique(sy$type)
  expect_true("Fuller relative bias" %in% types)
  expect_true("Fuller maximum bias" %in% types)
  expect_false("LIML size" %in% types)
})

test_that("ivreg2 kclass=0.5 has NULL Stock-Yogo", {
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card, kclass = 0.5)
  expect_null(fit$diagnostics$weak_id_sy)
})

test_that("2SLS Stock-Yogo unchanged after H4", {
  skip_if(!file.exists(card_path), "card data not found")

  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ | nearc4,
                data = card)
  sy <- fit$diagnostics$weak_id_sy
  expect_s3_class(sy, "data.frame")

  types <- unique(sy$type)
  expect_true("IV size" %in% types)
  expect_false("LIML size" %in% types)
})
