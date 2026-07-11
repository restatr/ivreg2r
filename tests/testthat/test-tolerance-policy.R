# Machine check of the tolerance-override policy: every tolerance literal in
# the test sources that is looser than standard_tol_ceiling must be an
# allow-listed override in helper-tolerances.R (tolerance_overrides), with a
# root-cause comment at the call site. Parser-based so that comments and
# strings mentioning tolerances are ignored. Deliberate consequence: a
# literal 1e-4 on a stat/pval comparison is flagged even though it equals
# stata_tol$stat — standard tolerances must be spelled stata_tol$..., which
# the scanner skips (symbolic, not a literal). If a comparison helper gains
# a new tolerance-like parameter, add its name to arg_names below.
#
# Known limitations (accepted; per-commit review is the backstop): the
# (file, value, count) allow-list key prevents inventory drift but is not a
# site identity — a same-file, same-value site swap in a single edit would
# balance the count; and the root-cause-comment requirement is convention,
# not machine-enforced (the scanner strips comments). A second scanner below
# closes the positional bypass: expect_equal(x, y, 1e-3) really does apply
# 1e-3 as tolerance through `...`, invisibly to the named-argument scan, so
# bare positional numeric literals beyond the first two arguments of
# expect_equal are banned outright.

scan_tolerance_literals <- function(path) {
  pd <- utils::getParseData(parse(path, keep.source = TRUE))
  # drop COMMENT tokens so an inline comment between "=" and the number
  # cannot break the name/EQ/NUM adjacency
  toks <- pd[pd$terminal & pd$token != "COMMENT",
             c("line1", "col1", "token", "text")]
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
             # strip the suffix of integer literals like 0L before coercing
             value = as.numeric(sub("L$", "", toks$text[i + 2L])))
}

# Bare positional numeric literals beyond the first two arguments of
# expect_equal(): each is an undeclared tolerance (testthat forwards it via
# `...` to waldo::compare, where it lands in `tolerance`).
scan_positional_numeric <- function(path) {
  pd <- utils::getParseData(parse(path, keep.source = TRUE))
  fns <- pd[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "expect_equal", ]
  if (nrow(fns) == 0L) return(NULL)
  hits <- list()
  for (k in seq_len(nrow(fns))) {
    call_id <- pd$parent[match(fns$parent[k], pd$id)]
    kids <- pd[pd$parent == call_id, ]
    kids <- kids[order(kids$line1, kids$col1), ]
    # collect positional argument exprs (skip named name/EQ_SUB/value runs;
    # the first expr child is the function itself)
    pos_args <- integer(0)
    i <- 1L
    seen_fn <- FALSE
    while (i <= nrow(kids)) {
      tok <- kids$token[i]
      if (tok %in% c("SYMBOL_SUB", "STR_CONST", "EQ_SUB")) {
        # named argument: skip the name, the =, and the value expr
        i <- i + if (tok == "EQ_SUB") 2L else 3L
        next
      }
      if (tok == "expr") {
        if (!seen_fn) seen_fn <- TRUE else pos_args <- c(pos_args, kids$id[i])
      }
      i <- i + 1L
    }
    if (length(pos_args) < 3L) next
    for (arg_id in pos_args[-(1:2)]) {
      # terminal descendants of this argument expression
      ids <- arg_id
      repeat {
        grown <- unique(c(ids, pd$id[pd$parent %in% ids]))
        if (length(grown) == length(ids)) break
        ids <- grown
      }
      terms <- pd[pd$id %in% ids & pd$terminal, ]
      is_num <- nrow(terms) == 1L && terms$token == "NUM_CONST"
      is_neg_num <- nrow(terms) == 2L && identical(sort(terms$token),
                                                   sort(c("'-'", "NUM_CONST")))
      if (is_num || is_neg_num) {
        hits[[length(hits) + 1L]] <- data.frame(
          file = basename(path),
          line = pd$line1[pd$id == arg_id]
        )
      }
    }
  }
  if (length(hits) == 0L) return(NULL)
  do.call(rbind, hits)
}

test_that("no bare positional numeric is passed to expect_equal beyond arg 2", {
  files <- list.files(test_path(), pattern = "^(test-|helper-).*\\.[rR]$",
                      full.names = TRUE)
  files <- files[basename(files) != "test-tolerance-policy.R"]
  skip_if(length(files) == 0L, "test sources not visible to the scanner")

  hits <- do.call(rbind, lapply(files, scan_positional_numeric))
  expect_true(is.null(hits) || nrow(hits) == 0L, info = paste0(
    "Positional numeric literals act as an undeclared tolerance in ",
    "expect_equal (forwarded via `...` to waldo::compare). Pass tolerance ",
    "by name instead. Offending sites:\n",
    if (!is.null(hits)) paste(sprintf("  %s:%d", hits$file, hits$line),
                              collapse = "\n") else ""))
})

test_that("every looser-than-standard tolerance literal is allow-listed", {
  files <- list.files(test_path(), pattern = "^(test-|helper-).*\\.[rR]$",
                      full.names = TRUE)
  files <- files[basename(files) != "test-tolerance-policy.R"]
  skip_if(length(files) == 0L, "test sources not visible to the scanner")

  hits <- do.call(rbind, lapply(files, scan_tolerance_literals))
  if (is.null(hits)) {
    hits <- data.frame(file = character(), line = integer(),
                       value = numeric())
  }
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
