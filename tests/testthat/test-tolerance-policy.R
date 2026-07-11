# Machine check of the tolerance-override policy: every tolerance literal in
# the test sources that is looser than standard_tol_ceiling must be an
# allow-listed override in helper-tolerances.R (tolerance_overrides), with a
# root-cause comment at the call site. Parser-based so that comments and
# strings mentioning tolerances are ignored.

scan_tolerance_literals <- function(path) {
  pd <- utils::getParseData(parse(path, keep.source = TRUE))
  toks <- pd[pd$terminal, c("line1", "col1", "token", "text")]
  toks <- toks[order(toks$line1, toks$col1), ]
  arg_names <- c("tolerance", "tol", "tol_coef", "tol_se", "tol_stat",
                 "tol_pval")
  i <- which(toks$token %in% c("SYMBOL_SUB", "SYMBOL_FORMALS") &
               toks$text %in% arg_names)
  i <- i[i + 2L <= nrow(toks) &
           toks$token[i + 1L] %in% c("EQ_SUB", "EQ_FORMALS") &
           toks$token[i + 2L] == "NUM_CONST"]
  if (length(i) == 0L) return(NULL)
  data.frame(file = basename(path),
             line = toks$line1[i + 2L],
             value = as.numeric(toks$text[i + 2L]))
}

test_that("every looser-than-standard tolerance literal is allow-listed", {
  files <- list.files(test_path(), pattern = "^(test-|helper-).*\\.[rR]$",
                      full.names = TRUE)
  files <- files[basename(files) != "test-tolerance-policy.R"]
  skip_if(length(files) == 0L, "test sources not visible to the scanner")

  hits <- do.call(rbind, lapply(files, scan_tolerance_literals))
  loose <- hits[hits$value > standard_tol_ceiling, , drop = FALSE]

  if (nrow(loose) == 0L) {
    observed <- data.frame(file = character(), value = numeric(),
                           count = integer())
  } else {
    observed <- aggregate(line ~ file + value, data = loose, FUN = length)
    names(observed)[names(observed) == "line"] <- "count"
  }
  # (a) every loose literal must match an allow-list row with the exact
  # occurrence count; (b) every allow-list row must still match real
  # occurrences (stale rows fail too).
  merged <- merge(observed, tolerance_overrides[c("file", "value", "count")],
                  by = c("file", "value"), all = TRUE,
                  suffixes = c("_found", "_allowed"))
  bad <- merged[is.na(merged$count_found) | is.na(merged$count_allowed) |
                  merged$count_found != merged$count_allowed, , drop = FALSE]
  expect_true(nrow(bad) == 0L, info = paste0(
    "Tolerance policy violations (see tolerance_overrides in ",
    "helper-tolerances.R):\n",
    paste(capture.output(print(bad)), collapse = "\n"),
    "\nAll looser-than-standard sites found:\n",
    paste(sprintf("  %s:%d value=%g", loose$file, loose$line, loose$value),
          collapse = "\n")))
})
