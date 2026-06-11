# --------------------------------------------------------------------------
# Time-series operators l() / d()
# --------------------------------------------------------------------------
# Formula-level lag and difference operators resolved against `tvar`/`ivar`,
# mirroring Stata's tsset-aware L./D. operators (ivreg2.ado relies on
# fvrevar/fvexpand for these, ivreg2.ado:490-493, 4021-4034).
#
# Mechanism (fixest's panel operators are the design precedent, with one
# improvement): `l()`/`d()` calls in the formula are evaluated during
# model.frame() against closures placed in an environment that wraps the
# formula's own environment (.make_ts_env()), rather than located by call-
# stack inspection as fixest does. Lag *ranges* (l(x, 1:3)) are expanded to
# scalar-k terms before model.frame() so each term yields exactly one column
# whose name is its canonical term label, e.g. "l(unem, 1)" — option matching
# (orthog/endog/redundant/partial), collinearity reporting, and broom output
# all flow through the existing term-label machinery unchanged.
#
# Lags are computed by time VALUE within panel unit (gap-aware): a gap at
# t - k, or the start of a unit's series, yields NA, and the row is dropped
# from the estimation sample by na.action — matching Stata's markout-based
# sample shrinkage.

.TS_OPERATORS <- c("l", "d")


# --------------------------------------------------------------------------
# Exported stubs
# --------------------------------------------------------------------------
#' Time-series operators for ivreg2 formulas
#'
#' `l()` (lag) and `d()` (difference) may be used inside an [ivreg2()] formula
#' to create lags and differences of a variable on the fly, resolved against
#' the time variable `tvar` (and panel variable `ivar` for panel data). They
#' mirror Stata's `tsset`-aware `L.` and `D.` operators:
#'
#' | R syntax | Meaning | Stata equivalent |
#' |---|---|---|
#' | `l(x)`, `l(x, 1)` | lag 1 | `L.x` |
#' | `l(x, 1:3)` | lags 1, 2, 3 (three terms) | `L(1/3).x` |
#' | `l(x, 0)` | `x` itself | `L0.x` |
#' | `d(x)`, `d(x, 1)` | first difference `x - l(x, 1)` | `D.x` |
#' | `d(x, 2)` | second-order difference `x - 2 l(x, 1) + l(x, 2)` | `D2.x` |
#'
#' @section Lag semantics:
#' Lags are computed by time **value**, not row position: `l(x, k)` at time
#' `t` is the value of `x` at time `t - k` within the same panel unit. If no
#' observation exists at `t - k` (start of the series, or a gap in the time
#' variable), the lag is `NA` and the row is dropped from the estimation
#' sample — exactly Stata's behavior. Lags never cross panel-unit boundaries.
#' The data does not need to be sorted.
#'
#' @section Differences from fixest:
#' Following Stata (not fixest), `d(x, k)` is the k-th **order** (iterated)
#' difference `(1-L)^k x`, matching Stata's `Dk.x`; fixest's `d(x, k)` is the
#' distance-k difference `x - l(x, k)`. The two agree for `k = 1`.
#'
#' @section Usage constraints:
#' Operators require `tvar` (and `ivar` for panel data) to be passed to
#' [ivreg2()]; `l()` and `d()` are not standalone lag utilities and signal an
#' error when called outside an `ivreg2()` formula (use, e.g.,
#' `dplyr::lag()` with a consecutive-periods check, or `fixest`/`collapse`/
#' `plm`, for general-purpose panel lagging). The lag/difference order `k`
#' must be a constant non-negative integer (or integer vector for ranges);
#' leads (negative `k`), nested operators (`l(d(x))`), and ranges on the
#' left-hand side (the model has a single dependent variable) are not
#' supported. When operators appear among the *regressors*,
#' `predict(fit, newdata)` recomputes them within `newdata`, so `newdata`
#' must contain the history rows needed for the lags; rows whose lags are
#' missing get `NA` predictions (again matching Stata). Operators confined
#' to the instrument part impose no requirement on `newdata`.
#'
#' @param x A numeric variable in the model data.
#' @param k Lag order(s) for `l()`; difference order(s) for `d()`.
#'   Non-negative integer (vector allowed, expanded to one term per element).
#' @return When used inside an [ivreg2()] formula: a numeric vector aligned
#'   to the model data. Calling these functions anywhere else is an error.
#' @examples
#' data(phillips)
#' # Stata: ivreg2 cinf (unem = l(1/3).unem), bw(3)
#' fit <- ivreg2(cinf ~ 1 | unem | l(unem, 1:3), data = phillips,
#'               tvar = "year", vcov = "AC", kernel = "bartlett", bw = 3)
#' coef(fit)
#' @name ts-operators
#' @aliases l d
NULL

#' @rdname ts-operators
#' @export
l <- function(x, k = 1) {
  stop("l() can only be used inside an ivreg2() formula, with `tvar` ",
       "(and `ivar` for panel data) specified. See ?`ts-operators`.",
       call. = FALSE)
}

#' @rdname ts-operators
#' @export
d <- function(x, k = 1) {
  stop("d() can only be used inside an ivreg2() formula, with `tvar` ",
       "(and `ivar` for panel data) specified. See ?`ts-operators`.",
       call. = FALSE)
}


# --------------------------------------------------------------------------
# Detection
# --------------------------------------------------------------------------
#' Does an expression contain a time-series operator call?
#'
#' Recursively scans a language object (typically a formula) for calls to
#' `l()` or `d()`. Symbols named `l`/`d` (e.g. a variable called `l`) do not
#' count — only calls.
#'
#' @param expr A language object (formula, call, symbol, or atomic).
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.has_ts_operators <- function(expr) {
  if (!is.call(expr)) {
    return(FALSE)
  }
  fn <- expr[[1L]]
  if (is.symbol(fn) && as.character(fn) %in% .TS_OPERATORS) {
    return(TRUE)
  }
  for (i in seq_along(expr)[-1L]) {
    arg <- expr[[i]]
    if (is.call(arg) && .has_ts_operators(arg)) {
      return(TRUE)
    }
  }
  FALSE
}


# --------------------------------------------------------------------------
# Normalization / range expansion
# --------------------------------------------------------------------------
#' Normalize ts-operator calls in a formula
#'
#' Expands lag ranges (`l(x, 1:3)` becomes `l(x, 1) + l(x, 2) + l(x, 3)`,
#' fixest's expand_lags() precedent) and canonicalizes every operator call to
#' the explicit two-argument form (`d(w)` becomes `d(w, 1)`), so that term
#' labels, column names, and user option strings all share one canonical
#' spelling. Validates `k` (constant, non-negative integer) and rejects
#' nested operators and leads.
#'
#' @param formula A formula (possibly multi-part, pre-`as.Formula`).
#' @return The normalized formula with its original environment.
#' @keywords internal
#' @noRd
.normalize_ts_formula <- function(formula) {
  env <- environment(formula)
  expr <- .normalize_ts_expr(formula, env, top = TRUE)
  stats::as.formula(expr, env = env)
}

#' Recursive worker for .normalize_ts_formula
#'
#' @param e A language object.
#' @param env Environment in which to evaluate `k` arguments.
#' @param top Logical: is `e` still in top-level term-sum position (under
#'   `~`, `+`, `-`, `|`, or `(`)? Lag ranges may only appear there — a range
#'   inside, say, `log()` has no well-defined expansion.
#' @param lhs Logical: is `e` on the left-hand side of the `~`? Lag ranges
#'   are rejected there — the model has a single dependent variable (Stata
#'   parity: "multiple dependent variables specified", ivreg2.ado:502-507).
#' @return The transformed language object.
#' @keywords internal
#' @noRd
.normalize_ts_expr <- function(e, env, top = TRUE, lhs = FALSE) {
  if (!is.call(e)) {
    return(e)
  }
  fn <- e[[1L]]
  fname <- if (is.symbol(fn)) as.character(fn) else ""
  if (fname %in% .TS_OPERATORS) {
    return(.normalize_ts_call(e, env, top, lhs))
  }
  top_child <- top && fname %in% c("~", "+", "-", "|", "(")
  args <- as.list(e)
  for (i in seq_along(args)[-1L]) {
    if (is.call(args[[i]])) {
      # In a two-sided `y ~ rhs` call, args[[2]] is the LHS
      lhs_child <- lhs || (fname == "~" && length(args) == 3L && i == 2L)
      args[[i]] <- .normalize_ts_expr(args[[i]], env, top = top_child,
                                      lhs = lhs_child)
    }
  }
  as.call(args)
}

#' Canonicalize a single l()/d() call, expanding ranges
#' @keywords internal
#' @noRd
.normalize_ts_call <- function(e, env, top, lhs = FALSE) {
  op <- as.character(e[[1L]])
  m <- tryCatch(
    match.call(function(x, k = 1) NULL, e),
    error = function(err) {
      stop("Invalid ", op, "() term in formula: ", deparse1(e), ". ",
           "Usage: ", op, "(x, k).", call. = FALSE)
    }
  )
  if (!"x" %in% names(as.list(m))) {
    stop(op, "() requires a variable: ", deparse1(e), ".", call. = FALSE)
  }
  x_expr <- m$x
  if (is.call(x_expr) && .has_ts_operators(x_expr)) {
    stop("Nested time-series operators are not supported: ", deparse1(e), ".",
         call. = FALSE)
  }
  k <- if (is.null(m$k)) {
    1
  } else {
    tryCatch(
      eval(m$k, env),
      error = function(err) {
        stop("Could not evaluate the lag/difference order in ", deparse1(e),
             ": k must be a constant (e.g. 1 or 1:3).", call. = FALSE)
      }
    )
  }
  if (!is.numeric(k) || length(k) < 1L || anyNA(k) || any(!is.finite(k)) ||
      any(k != trunc(k))) {
    stop("The lag/difference order in ", deparse1(e),
         " must be integer-valued.", call. = FALSE)
  }
  if (any(k < 0)) {
    stop("Negative orders in ", op, "() are not supported ",
         "(leads are out of scope): ", deparse1(e), ".", call. = FALSE)
  }
  calls <- lapply(as.numeric(k), function(ki) {
    as.call(list(as.name(op), x_expr, ki))
  })
  if (length(calls) == 1L) {
    return(calls[[1L]])
  }
  if (lhs) {
    stop("Lag ranges are not allowed on the left-hand side of the formula ",
         "(the model has a single dependent variable): ", deparse1(e), ".",
         call. = FALSE)
  }
  if (!top) {
    stop("Lag ranges like ", op, "(x, 1:3) can only be used as top-level ",
         "formula terms, not inside other functions: ", deparse1(e), ".",
         call. = FALSE)
  }
  Reduce(function(a, b) call("+", a, b), calls)
}


# --------------------------------------------------------------------------
# Evaluation context
# --------------------------------------------------------------------------
#' Validate the time/panel variables for ts-operator resolution
#'
#' Stata analogue: tsset's own checks, which ivreg2 relies on
#' (ivreg2.ado:315-320). Stata's tsset rejects missing time values and
#' repeated time values (within panel); lags are computed on the full
#' dataset before any if/in restriction, so the checks run on the full data.
#'
#' @param data A data frame (the estimation data, or `newdata` in predict).
#' @param tvar,ivar Column names (ivar may be NULL).
#' @param what Label for error messages ("data" or "newdata").
#' @return List with `tvar_vec` and `ivar_vec` (full-length vectors).
#' @keywords internal
#' @noRd
.ts_validate_context <- function(data, tvar, ivar = NULL, what = "data") {
  if (!tvar %in% names(data)) {
    stop("Time variable '", tvar, "' not found in ", what, ".", call. = FALSE)
  }
  tvar_vec <- data[[tvar]]
  if (anyNA(tvar_vec)) {
    stop("Time variable '", tvar, "' contains NA values.", call. = FALSE)
  }
  if (!is.numeric(tvar_vec)) {
    stop("Time variable '", tvar, "' must be numeric.", call. = FALSE)
  }
  # Inf - k == Inf, so an infinite time value would match() as its own lag
  # (silent self-lag). Stata cannot represent infinity, so tsset's effective
  # guarantee includes finiteness.
  if (any(!is.finite(tvar_vec))) {
    stop("Time variable '", tvar, "' must contain only finite values.",
         call. = FALSE)
  }
  ivar_vec <- NULL
  if (!is.null(ivar)) {
    if (!ivar %in% names(data)) {
      stop("Panel variable '", ivar, "' not found in ", what, ".",
           call. = FALSE)
    }
    ivar_vec <- data[[ivar]]
    if (anyNA(ivar_vec)) {
      stop("Panel variable '", ivar, "' contains NA values.", call. = FALSE)
    }
  }
  if (is.null(ivar_vec)) {
    if (anyDuplicated(tvar_vec)) {
      stop("Time variable '", tvar, "' has repeated time values. ",
           "Time-series operators require unique time values ",
           "(specify `ivar` for panel data).", call. = FALSE)
    }
  } else {
    if (anyDuplicated(data.frame(ivar_vec, tvar_vec))) {
      stop("Repeated time values within panel: (", ivar, ", ", tvar,
           ") pairs must be unique for time-series operators.",
           call. = FALSE)
    }
  }
  list(tvar_vec = tvar_vec, ivar_vec = ivar_vec)
}

#' Build the evaluation environment for in-formula l()/d()
#'
#' Returns a new environment whose parent is `parent_env` (normally the
#' formula's own environment) containing working `l`/`d` closures bound to
#' the supplied time/panel vectors. model.frame() evaluates formula terms
#' with this environment as the enclosure, so the closures shadow the
#' exported error stubs (and any user/fixest function of the same name)
#' during estimation only.
#'
#' @param parent_env Environment to wrap.
#' @param tvar_vec,ivar_vec Full-length time/panel vectors aligned to the
#'   rows of the data the formula will be evaluated on.
#' @return An environment containing `l` and `d`.
#' @keywords internal
#' @noRd
.make_ts_env <- function(parent_env, tvar_vec, ivar_vec = NULL) {
  env <- new.env(parent = parent_env)
  env$l <- function(x, k = 1) {
    .check_ts_k_scalar(k, "l")
    .ts_lag_values(x, k, tvar_vec, ivar_vec)
  }
  env$d <- function(x, k = 1) {
    .check_ts_k_scalar(k, "d")
    out <- x
    for (i in seq_len(k)) {
      out <- out - .ts_lag_values(out, 1, tvar_vec, ivar_vec)
    }
    out
  }
  env
}

#' Defensive scalar-k check for the in-formula closures
#'
#' Ranges are expanded before model.frame(), so by the time these closures
#' run k is scalar — unless the user evaluates the formula through a path
#' that skipped normalization.
#' @keywords internal
#' @noRd
.check_ts_k_scalar <- function(k, op) {
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || !is.finite(k) ||
      k != trunc(k) || k < 0) {
    stop(op, "() requires a single non-negative integer order.",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Gap-aware lag by time value within panel unit
#'
#' For each row i, returns the value of `x` at time `tvar_vec[i] - k` within
#' the same panel unit, or NA if no such observation exists (series start or
#' gap). Pure match()-based — no sorting, vectorized, by time value.
#'
#' @param x Numeric vector (full data length).
#' @param k Scalar non-negative integer lag.
#' @param tvar_vec Numeric time vector, same length as `x`.
#' @param ivar_vec Panel id vector or NULL.
#' @return Numeric vector of lagged values with NA where undefined.
#' @keywords internal
#' @noRd
.ts_lag_values <- function(x, k, tvar_vec, ivar_vec = NULL) {
  if (!is.numeric(x)) {
    stop("Time-series operators require numeric variables.", call. = FALSE)
  }
  if (length(x) != length(tvar_vec)) {
    stop("Time-series operator applied to a vector whose length (",
         length(x), ") differs from the time variable (", length(tvar_vec),
         ").", call. = FALSE)
  }
  if (k == 0) {
    return(x)
  }
  if (is.null(ivar_vec)) {
    idx <- match(tvar_vec - k, tvar_vec)
  } else {
    idx <- match(paste(ivar_vec, tvar_vec - k, sep = "\r"),
                 paste(ivar_vec, tvar_vec, sep = "\r"))
  }
  unname(x[idx])
}


# --------------------------------------------------------------------------
# Option-string canonicalization
# --------------------------------------------------------------------------
#' Canonicalize user-supplied term labels containing ts operators
#'
#' Applies the same normalization as the formula pass to character vectors
#' from `endog`/`orthog`/`redundant`/`partial`, so that `"l(unem,1)"` and
#' `"l(unem, 1)"` both match the canonical term label, and ranges like
#' `"l(unem, 1:2)"` expand to their component terms (Stata analogue:
#' fvexpand canonicalizes `l1.unem` before option-varlist matching,
#' ivreg2.ado:511-530). Strings that do not parse, or contain no ts
#' operator, pass through unchanged.
#'
#' @param x Character vector (or NULL).
#' @return Character vector with operator entries canonicalized/expanded.
#' @keywords internal
#' @noRd
.canonicalize_ts_labels <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  out <- character(0L)
  for (s in x) {
    e <- tryCatch(str2lang(s), error = function(err) NULL)
    if (is.null(e) || !is.call(e) || !.has_ts_operators(e)) {
      out <- c(out, s)
      next
    }
    norm <- tryCatch(.normalize_ts_expr(e, baseenv(), top = TRUE),
                     error = function(err) NULL)
    if (is.null(norm)) {
      out <- c(out, s)
      next
    }
    out <- c(out, .flatten_term_sum(norm))
  }
  out
}

#' Flatten a `+` chain of calls into deparsed labels
#' @keywords internal
#' @noRd
.flatten_term_sum <- function(e) {
  if (is.call(e) && identical(e[[1L]], as.name("+"))) {
    c(.flatten_term_sum(e[[2L]]), .flatten_term_sum(e[[3L]]))
  } else {
    deparse1(e)
  }
}
