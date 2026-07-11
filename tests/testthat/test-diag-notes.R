# ============================================================================
# Tests: Stata-quirk disclosure notes on diagnostic slots (planning/31 R2)
# ============================================================================
# Five sites where an explicit user option is silently not honored inside an
# automatically computed diagnostic carry a "note-as-data" disclosure: a note
# on the affected slot, surfaced by diagnostics() as a `note` column and by
# summary() as a footnote. These tests assert each note fires on the trigger
# and is absent without it. No Stata fixtures are needed: the notes are
# metadata, not numbers, and the underlying statistics are validated
# elsewhere.

data(card)
data(phillips)

# Substrings that pin each note without hard-coding the full sentence.
KERNEL_SUB <- "Bartlett kernel regardless"
KIEFER_SUB <- "Kiefer VCE structure"
PSD_SUB    <- "psd correction"
CENTER_SUB <- "without moment centering"
CUE_SUB    <- "does not use CUE"
ORTHOG_SUB <- "does not use the CUE or LIML"

# Synthetic panel index for the kiefer fits (mirrors test-weight-types.R).
card_panel <- card
card_panel$id2 <- ceiling(seq_len(nrow(card)) / 10)
card_panel$tt <- (seq_len(nrow(card)) - 1) %% 10 + 1


# ============================================================================
# Trigger 1: non-Bartlett kernel -> KP id + redundancy notes
# ============================================================================

test_that("a non-Bartlett kernel notes the id and redundancy tests", {
  fit <- ivreg2(cinf ~ 1 | unem | l(unem, 1:3), data = phillips,
                tvar = "year", vcov = "HAC", kernel = "Parzen", bw = 3,
                redundant = "l(unem, 3)")

  expect_match(fit$diagnostics$underid$note, KERNEL_SUB)
  expect_match(fit$diagnostics$weak_id_robust$note, KERNEL_SUB)
  expect_match(fit$diagnostics$redundancy$note, KERNEL_SUB)

  d <- diagnostics(fit)
  expect_match(d$note[d$test == "underid"], KERNEL_SUB)
  expect_match(d$note[d$test == "weak_id_robust"], KERNEL_SUB)
  expect_match(d$note[d$test == "redundancy"], KERNEL_SUB)
})

test_that("a Bartlett-kernel fit carries no kernel note", {
  fit <- ivreg2(cinf ~ 1 | unem | l(unem, 1:3), data = phillips,
                tvar = "year", vcov = "HAC", kernel = "Bartlett", bw = 3,
                redundant = "l(unem, 3)")

  expect_null(fit$diagnostics$underid$note)
  expect_null(fit$diagnostics$weak_id_robust$note)
  expect_null(fit$diagnostics$redundancy$note)
  expect_true(all(is.na(diagnostics(fit)$note)))
})


# ============================================================================
# Trigger 2: center = TRUE -> endogeneity note
# ============================================================================

test_that("center = TRUE notes the endogeneity test", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ |
                  nearc4 + nearc2, data = card, vcov = "robust",
                center = TRUE)

  expect_match(fit$diagnostics$endogeneity$note, CENTER_SUB)
  expect_match(diagnostics(fit)$note[diagnostics(fit)$test == "endogeneity"],
               CENTER_SUB)
})

test_that("center = FALSE leaves the endogeneity test unannotated", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ |
                  nearc4 + nearc2, data = card, vcov = "robust")

  expect_null(fit$diagnostics$endogeneity$note)
})


# ============================================================================
# Trigger 3: kiefer -> id-test notes
# ============================================================================

test_that("a kiefer fit notes the identification tests", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4 + nearc2,
                data = card_panel, kiefer = TRUE, tvar = "tt", ivar = "id2")

  expect_match(fit$diagnostics$underid$note, KIEFER_SUB)
  expect_match(fit$diagnostics$weak_id$note, KIEFER_SUB)

  d <- diagnostics(fit)
  expect_match(d$note[d$test == "underid"], KIEFER_SUB)
  expect_match(d$note[d$test == "weak_id"], KIEFER_SUB)
})

test_that("a non-kiefer fit carries no kiefer note", {
  fit <- ivreg2(lwage ~ exper + expersq | educ | nearc4 + nearc2,
                data = card_panel)

  expect_null(fit$diagnostics$underid$note)
  expect_null(fit$diagnostics$weak_id$note)
})


# ============================================================================
# Trigger 4: psd on a robust fit -> KP id + redundancy notes
# ============================================================================

test_that("psd on a robust fit notes the id and redundancy tests", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ |
                  nearc4 + nearc2, data = card, vcov = "robust",
                psd = "psd0", redundant = "nearc2")

  expect_match(fit$diagnostics$underid$note, PSD_SUB)
  expect_match(fit$diagnostics$weak_id_robust$note, PSD_SUB)
  expect_match(fit$diagnostics$redundancy$note, PSD_SUB)

  d <- diagnostics(fit)
  expect_match(d$note[d$test == "underid"], PSD_SUB)
  expect_match(d$note[d$test == "weak_id_robust"], PSD_SUB)
  expect_match(d$note[d$test == "redundancy"], PSD_SUB)
})

test_that("a robust fit without psd carries no psd note", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ |
                  nearc4 + nearc2, data = card, vcov = "robust",
                redundant = "nearc2")

  expect_null(fit$diagnostics$underid$note)
  expect_null(fit$diagnostics$weak_id_robust$note)
  expect_null(fit$diagnostics$redundancy$note)
})


# ============================================================================
# Trigger 5: method = "cue"/"liml" -> endogeneity + orthogonality notes
# ============================================================================

test_that("a CUE fit notes the endogeneity C-statistic", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ |
                  nearc4 + nearc2, data = card, vcov = "robust",
                method = "cue")

  expect_match(fit$diagnostics$endogeneity$note, CUE_SUB)
  expect_match(diagnostics(fit)$note[diagnostics(fit)$test == "endogeneity"],
               CUE_SUB)
})

test_that("a 2SLS fit carries no CUE endogeneity note", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ |
                  nearc4 + nearc2, data = card, vcov = "robust")

  expect_false(any(grepl(CUE_SUB, fit$diagnostics$endogeneity$note)))
})

test_that("CUE and LIML fits note the orthogonality C-statistic", {
  fit_cue <- ivreg2(lwage ~ exper + expersq + black + south | educ |
                      nearc4 + nearc2, data = card, vcov = "robust",
                    method = "cue", orthog = "nearc2")
  fit_liml <- ivreg2(lwage ~ exper + expersq + black + south | educ |
                       nearc4 + nearc2, data = card, method = "liml",
                     orthog = "nearc2")

  expect_match(fit_cue$diagnostics$orthog$note, ORTHOG_SUB)
  expect_match(fit_liml$diagnostics$orthog$note, ORTHOG_SUB)
  expect_match(diagnostics(fit_cue)$note[diagnostics(fit_cue)$test == "orthog"],
               ORTHOG_SUB)
})

test_that("a 2SLS fit carries no orthogonality recursive-estimation note", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ |
                  nearc4 + nearc2, data = card, vcov = "robust",
                orthog = "nearc2")

  expect_null(fit$diagnostics$orthog$note)
})


# ============================================================================
# summary() prints each distinct note exactly once
# ============================================================================

test_that("summary() prints a shared note once as a footnote", {
  # The kernel note is attached to three slots (underid, weak_id_robust,
  # redundancy); the footnote block must collapse them to a single line.
  fit <- ivreg2(cinf ~ 1 | unem | l(unem, 1:3), data = phillips,
                tvar = "year", vcov = "HAC", kernel = "Parzen", bw = 3,
                redundant = "l(unem, 3)")
  out <- capture.output(print(summary(fit)))

  expect_equal(sum(grepl(KERNEL_SUB, out)), 1L)
  expect_true(any(grepl("^Note: ", out)))
})

test_that("summary() prints no footnote when no note is triggered", {
  fit <- ivreg2(lwage ~ exper + expersq + black + south | educ |
                  nearc4 + nearc2, data = card, vcov = "robust")
  out <- capture.output(print(summary(fit)))

  expect_false(any(grepl("^Note: ", out)))
})
