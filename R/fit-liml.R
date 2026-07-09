# --------------------------------------------------------------------------
# .fit_kclass
# --------------------------------------------------------------------------
#' Fit k-class estimator (LIML, Fuller, user-supplied k)
#'
#' Computes k-class IV estimates. When `method = "liml"`, the LIML eigenvalue
#' lambda is computed from the concentration matrices. When `fuller > 0`,
#' applies the Fuller (1977) bias correction `k = lambda - fuller / (N - L)`.
#' When `method = "kclass"`, uses the user-supplied `kclass` value directly.
#'
#' The 2SLS-style bread `(X_hat'WX_hat)^{-1}` is returned for diagnostics
#' and downstream HC/CL VCV computation (H2), alongside the k-class VCV
#' `sigma^2 * solve(XhXh)`.
#'
#' @importFrom stats lm.fit lm.wfit
#' @param parsed A `parsed_formula` object from `.parse_formula()`.
#' @param method Character: `"liml"` or `"kclass"`.
#' @param kclass Numeric scalar: user-supplied k value (used only when
#'   `method = "kclass"`).
#' @param fuller Numeric scalar: Fuller modification parameter (default 0).
#' @param small Logical: if `TRUE`, use `N-K` denominator for sigma.
#' @param dofminus Integer: large-sample DoF adjustment (default 0).
#' @param sdofminus Integer: small-sample DoF adjustment (default 0).
#' @return A named list with `coefficients`, `residuals`, `fitted.values`,
#'   `vcov`, `sigma`, `df.residual`, `rank`, `r.squared`, `adj.r.squared`,
#'   `rss`, `r2u`, `r2c`, `mss`, `bread`, `bread_kclass`, `X_hat`, `proj_coef`,
#'   `lambda`,
#'   `kclass_value`, `method`, `fuller_param`.
#' @keywords internal
.fit_kclass <- function(parsed, method = "liml", kclass = NULL, fuller = 0,
                        small = FALSE, dofminus = 0L, sdofminus = 0L) {
  y <- parsed$y
  X <- parsed$X
  Z <- parsed$Z
  N <- parsed$N
  K <- parsed$K
  K1 <- parsed$K1
  L <- parsed$L
  L1 <- parsed$L1
  has_intercept <- parsed$has_intercept
  w <- parsed$weights

  # Plain just-identified LIML is mathematically equal to 2SLS (k = lambda = 1).
  # Short-circuit so the result is bit-identical to .fit_2sls() rather than
  # drifting at ulp scale through a different solve path.
  if (method == "liml" && L == K &&
      (is.null(fuller) || fuller == 0) && is.null(kclass)) {
    fit_2sls <- .fit_2sls(parsed, small = small,
                          dofminus = dofminus, sdofminus = sdofminus)
    fit_2sls$bread_kclass <- fit_2sls$bread
    fit_2sls$lambda <- 1.0
    fit_2sls$kclass_value <- 1.0
    fit_2sls$method <- "liml"
    fit_2sls$fuller_param <- 0
    return(fit_2sls)
  }

  # --- Stage 1: project X onto column space of Z (same as 2SLS) ---
  if (is.null(w)) {
    lm_Z <- lm.fit(Z, X)
  } else {
    lm_Z <- lm.wfit(Z, X, w)
  }
  X_hat <- as.matrix(lm_Z$fitted.values)
  colnames(X_hat) <- colnames(X)

  # --- First-stage map A (L x K): X_hat = Z %*% A ---
  # NA (aliased) coefficients zero-filled; see .fit_2sls().
  proj_coef <- as.matrix(lm_Z$coefficients)
  proj_coef[is.na(proj_coef)] <- 0
  dimnames(proj_coef) <- list(colnames(Z), colnames(X))

  # --- Stage 2: regress y on X_hat for rank check and 2SLS bread ---
  if (is.null(w)) {
    lm_XZ <- lm.fit(X_hat, y)
  } else {
    lm_XZ <- lm.wfit(X_hat, y, w)
  }

  if (lm_XZ$rank < K) {
    stop("The projected regressor matrix is rank-deficient (rank ",
         lm_XZ$rank, " < ", K, " regressors), so the instruments may not ",
         "identify the endogenous regressors.",
         call. = FALSE)
  }

  # 2SLS-style bread: (X_hat'WX_hat)^{-1} via stored QR
  p <- lm_XZ$rank
  R_qr <- lm_XZ$qr$qr[seq_len(p), seq_len(p), drop = FALSE]
  XtPX_inv <- chol2inv(R_qr)
  piv <- lm_XZ$qr$pivot
  piv_order <- order(piv[seq_len(p)])
  XtPX_inv <- XtPX_inv[piv_order, piv_order, drop = FALSE]
  colnames(XtPX_inv) <- rownames(XtPX_inv) <-
    colnames(X)[sort(piv[seq_len(p)])]

  # --- Compute LIML eigenvalue (when method == "liml") ---
  lambda <- NA_real_

  if (method == "liml") {
    if (L == K) {
      # Exactly identified: LIML = 2SLS, lambda = 1
      lambda <- 1.0
    } else {
      lambda <- .compute_liml_lambda(y, X, Z, parsed, w)
    }

    # Guard: spurious lambda near zero
    if (!is.na(lambda) && lambda < 1e-8) {
      warning("The LIML eigenvalue (lambda = ", format(lambda, digits = 4),
              ") is near zero, so the estimates may be unreliable.",
              call. = FALSE)
    }
  }

  # --- Determine k ---
  if (!is.null(kclass)) {
    k <- kclass
  } else if (fuller > 0) {
    k <- lambda - fuller / (N - L)
  } else {
    k <- lambda
  }

  # --- k-class cross-products ---
  if (is.null(w)) {
    XtX   <- crossprod(X)
    XtPzX <- crossprod(X_hat)
    Xty   <- crossprod(X, y)
    XtPzy <- crossprod(X_hat, y)
  } else {
    sw    <- sqrt(w)
    X_w   <- sw * X
    Xh_w  <- sw * X_hat
    y_w   <- sw * y
    XtX   <- crossprod(X_w)
    XtPzX <- crossprod(Xh_w)
    Xty   <- crossprod(X_w, y_w)
    XtPzy <- crossprod(Xh_w, y_w)
  }

  # k-class normal equation matrix: (1-k)*X'WX + k*X_hat'WX_hat
  XhXh <- (1 - k) * XtX + k * XtPzX
  XhXh <- (XhXh + t(XhXh)) / 2  # enforce symmetry
  Xhy  <- (1 - k) * Xty + k * XtPzy

  coef <- tryCatch(
    drop(.chol_solve(XhXh, Xhy)),
    error = function(e) {
      stop("The k-class normal-equation matrix is singular or nearly singular. ",
           "Check for collinearity among the regressors and instruments.",
           call. = FALSE)
    }
  )
  names(coef) <- colnames(X)

  # --- Residuals use original X ---
  fitted <- drop(X %*% coef)
  names(fitted) <- names(y)
  resid <- y - fitted
  rss <- if (is.null(w)) drop(crossprod(resid)) else sum(w * resid^2)

  # --- Sigma ---
  denom <- if (small) N - K - dofminus - sdofminus else N - dofminus
  sigma <- sqrt(rss / denom)

  # --- k-class IID VCV: sigma^2 * solve(XhXh) ---
  XhXh_inv <- tryCatch(
    .chol_solve(XhXh, diag(ncol(XhXh))),
    error = function(e) {
      stop("The k-class normal-equation matrix is singular or nearly singular. ",
           "Check for collinearity among the regressors and instruments.",
           call. = FALSE)
    }
  )
  XhXh_inv <- (XhXh_inv + t(XhXh_inv)) / 2
  colnames(XhXh_inv) <- rownames(XhXh_inv) <- colnames(X)
  V <- sigma^2 * XhXh_inv
  colnames(V) <- rownames(V) <- colnames(X)

  # --- R-squared (centered and uncentered) ---
  if (is.null(w)) {
    tss_c <- sum((y - mean(y))^2)
    tss_u <- sum(y^2)
  } else {
    wmean <- sum(w * y) / sum(w)
    tss_c <- sum(w * (y - wmean)^2)
    tss_u <- sum(w * y^2)
  }
  r2c <- 1 - rss / tss_c
  r2u <- 1 - rss / tss_u
  r2  <- if (has_intercept) r2c else r2u
  tss <- if (has_intercept) tss_c else tss_u
  mss <- tss - rss

  # --- Adjusted R-squared ---
  adj_r2 <- if (has_intercept) {
    1 - (1 - r2) * (N - 1) / (N - K - dofminus - sdofminus)
  } else {
    1 - (1 - r2) * N / (N - K - dofminus - sdofminus)
  }

  list(
    coefficients  = coef,
    residuals     = resid,
    fitted.values = fitted,
    vcov          = V,
    sigma         = sigma,
    df.residual   = as.integer(N - K - dofminus - sdofminus),
    rank          = as.integer(lm_XZ$rank),
    r.squared     = r2,
    adj.r.squared = adj_r2,
    rss           = rss,
    r2u           = r2u,
    r2c           = r2c,
    mss           = mss,
    bread         = XtPX_inv,
    bread_kclass  = XhXh_inv,
    X_hat         = X_hat,
    proj_coef     = proj_coef,
    lambda        = lambda,
    kclass_value  = k,
    method        = method,
    fuller_param  = fuller,
    tss_u         = tss_u,
    tss_c         = tss_c
  )
}


# --------------------------------------------------------------------------
# .compute_liml_lambda
# --------------------------------------------------------------------------
#' Compute the LIML eigenvalue from concentration matrices
#'
#' Constructs Y = (y, X1) (response + endogenous regressors), computes the
#' residual concentration matrix QWW = Y'M_Z Y and the restricted
#' concentration matrix QWW1 (where Z2 = included instruments only), then
#' returns the minimum eigenvalue of the symmetrized generalized eigenproblem.
#'
#' @param y Response vector.
#' @param X Full regressor matrix.
#' @param Z Full instrument matrix.
#' @param parsed Parsed formula object (for endo_names, excluded_names).
#' @param w Weight vector (NULL for unweighted).
#' @return Numeric scalar: LIML eigenvalue lambda.
#' @keywords internal
.compute_liml_lambda <- function(y, X, Z, parsed, w) {
  # Y = [y, X1] where X1 = endogenous regressor columns
  endo_cols <- match(parsed$endo_colnames, colnames(X))
  X1 <- X[, endo_cols, drop = FALSE]
  Y_mat <- cbind(y, X1)

  # Z2 = included instruments (exogenous regressors in Z)
  incl_cols <- which(!colnames(Z) %in% parsed$excluded_colnames)

  if (is.null(w)) {
    YtY <- crossprod(Y_mat)
    ZtY <- crossprod(Z, Y_mat)
    ZtZ <- crossprod(Z)
    QWW <- tryCatch(
      YtY - crossprod(ZtY, .chol_solve(ZtZ, ZtY)),
      error = function(e) {
        stop("The instrument cross-product matrix is singular. ",
             "Check for collinear instruments.", call. = FALSE)
      }
    )

    if (length(incl_cols) > 0L) {
      Z2  <- Z[, incl_cols, drop = FALSE]
      Z2tY  <- crossprod(Z2, Y_mat)
      Z2tZ2 <- crossprod(Z2)
      QWW1 <- tryCatch(
        YtY - crossprod(Z2tY, .chol_solve(Z2tZ2, Z2tY)),
        error = function(e) {
          stop("The included-instrument cross-product matrix is singular. ",
               "Check for collinear instruments.", call. = FALSE)
        }
      )
    } else {
      QWW1 <- YtY
    }
  } else {
    sw   <- sqrt(w)
    Y_w  <- sw * Y_mat
    Z_w  <- sw * Z

    YtWY <- crossprod(Y_w)
    ZtWY <- crossprod(Z_w, Y_w)
    ZtWZ <- crossprod(Z_w)
    QWW  <- tryCatch(
      YtWY - crossprod(ZtWY, .chol_solve(ZtWZ, ZtWY)),
      error = function(e) {
        stop("The instrument cross-product matrix is singular. ",
             "Check for collinear instruments.", call. = FALSE)
      }
    )

    if (length(incl_cols) > 0L) {
      Z2_w    <- Z_w[, incl_cols, drop = FALSE]
      Z2tWY   <- crossprod(Z2_w, Y_w)
      Z2tWZ2  <- crossprod(Z2_w)
      QWW1 <- tryCatch(
        YtWY - crossprod(Z2tWY, .chol_solve(Z2tWZ2, Z2tWY)),
        error = function(e) {
          stop("The included-instrument cross-product matrix is singular. ",
               "Check for collinear instruments.", call. = FALSE)
        }
      )
    } else {
      QWW1 <- YtWY
    }
  }

  # Enforce symmetry
  QWW  <- (QWW  + t(QWW))  / 2
  QWW1 <- (QWW1 + t(QWW1)) / 2

  # QWW^{-1/2} via eigendecomposition
  eig <- eigen(QWW, symmetric = TRUE)
  # Tiny negative eigenvalues are round-off on the (PSD) residual concentration
  # matrix Y'M_Z Y; clamp them silently. Warn only when a negative eigenvalue is
  # too large to be noise, which points to severe collinearity among the
  # instruments (an ill-conditioned Z'Z in the projection).
  cl <- .clamp_psd_eigenvalues(eig$values)
  d <- cl$values
  if (cl$non_psd) {
    warning("The LIML residual concentration matrix (Y'M_Z Y) was not ",
            "positive semidefinite (most negative eigenvalue ",
            formatC(cl$most_negative, format = "g", digits = 3),
            ") in the LIML eigenvalue computation; this usually indicates ",
            "severe collinearity among the instruments.", call. = FALSE)
  }
  inv_sqrt_d <- ifelse(d > 0, 1 / sqrt(d), 0)
  QWW_inv_sqrt <- eig$vectors %*% (inv_sqrt_d * t(eig$vectors))

  # Symmetrized eigenvalue problem: min eigenvalue of QWW^{-1/2} QWW1 QWW^{-1/2}
  M <- QWW_inv_sqrt %*% QWW1 %*% QWW_inv_sqrt
  M <- (M + t(M)) / 2

  min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)
}
