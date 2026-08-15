#' Ledoit-Wolf nonlinear shrinkage with optional light mixing
#'
#' Apply QIS nonlinear shrinkage to individual-level data and optionally blend
#' it with the sample correlation matrix.
#'
#' @param X Individual-level data matrix with observations in rows.
#' @param shrinkage Mixing weight for the nonlinear shrinkage estimator.
#'   Default is 0, which returns the positive-definite sample correlation
#'   matrix. Use 1 for the full nonlinear shrinkage estimator.
#' @param eigenmin Minimum eigenvalue target. Default is 0.001.
#' @param scale If `TRUE`, center and standardize `X` before estimation.
#'   Default is `FALSE`, assuming `X` has already been scaled.
#'
#' @return A positive-definite nonlinear-shrinkage correlation matrix.
#' @export
nonlinear_shrinkage <- function(X, shrinkage = 0, eigenmin = 1e-3, scale = FALSE) {
  if (missing(X) || is.null(X)) {
    stop("nonlinear_shrinkage requires individual-level X.", call. = FALSE)
  }
  prepared <- .ld_prepare_matrix(X, "X", scale = scale)
  S_sample <- .ld_matrix_cor(prepared$X, scale = FALSE)
  S_nl <- .ld_nonlinear_qis(prepared$X, eigenmin = eigenmin)
  S_nl <- stats::cov2cor(S_nl)
  S_nl[is.na(S_nl)] <- 0
  S_nl <- .ld_symmetrize(S_nl)
  diag(S_nl) <- 1

  if (length(shrinkage) != 1L || !is.finite(shrinkage)) {
    stop("shrinkage must be a finite scalar.", call. = FALSE)
  }
  shrinkage <- min(max(shrinkage[1L], 0), 1)
  out <- (1 - shrinkage) * S_sample + shrinkage * S_nl
  out <- .ld_fspd(.ld_symmetrize(out), eigenmin = eigenmin)
  out <- stats::cov2cor(out)
  out[is.na(out)] <- 0
  out <- .ld_symmetrize(out)
  diag(out) <- 1
  out
}

.ld_eigen_ascending <- function(x) {
  eig <- .ld_eigen(x)
  ord <- order(eig$values)
  list(values = eig$values[ord], vectors = eig$vectors[, ord, drop = FALSE])
}

.ld_reconstruct_from_eigen <- function(vectors, values) {
  out <- CppMatrix::matrixMultiply(t(t(vectors) * values), vectors, transB = TRUE)
  .ld_symmetrize(out)
}

# Adapted from covShrinkage QIS by Ledoit and Wolf (MIT License, 2022).
.ld_nonlinear_qis <- function(X, eigenmin = 1e-3) {
  X <- .ld_as_matrix(X, "X")
  N <- nrow(X)
  p <- ncol(X)
  if (N < 2L || p < 1L) {
    stop("X must have at least two rows and one column.", call. = FALSE)
  }
  n <- N - 1L
  if (n <= 1L) {
    stop("The effective sample size must be greater than 1.", call. = FALSE)
  }

  c_ratio <- p / n
  sample_cov <- CppMatrix::matrixMultiply(X, X, transA = TRUE) / n
  sample_cov <- .ld_symmetrize(sample_cov)
  spectral <- .ld_eigen_ascending(sample_cov)
  lambda <- pmax(spectral$values, 0)
  u <- spectral$vectors

  h <- min(c_ratio^2, 1 / c_ratio^2)^0.35 / p^0.35
  nonzero_index <- max(1L, p - n + 1L):p
  lambda_nonzero <- pmax(lambda[nonzero_index], max(eigenmin, .Machine$double.eps))
  invlambda <- 1 / lambda_nonzero

  m <- min(p, n)
  Lj <- matrix(rep(invlambda, each = m), nrow = m)
  Lj_i <- Lj - t(Lj)
  theta <- rowMeans(Lj * Lj_i / (Lj_i^2 + h^2 * Lj^2))
  Htheta <- rowMeans(Lj * (h * Lj) / (Lj_i^2 + h^2 * Lj^2))
  Atheta2 <- theta^2 + Htheta^2

  if (p <= n) {
    delta <- 1 / (
      (1 - c_ratio)^2 * invlambda +
        2 * c_ratio * (1 - c_ratio) * invlambda * theta +
        c_ratio^2 * invlambda * Atheta2
    )
  } else {
    delta0 <- 1 / ((c_ratio - 1) * mean(invlambda))
    delta <- c(rep(delta0, p - n), 1 / (invlambda * Atheta2))
  }
  delta <- delta * (sum(lambda) / sum(delta))
  delta <- pmax(as.numeric(delta), eigenmin)
  .ld_reconstruct_from_eigen(u, delta)
}
.ld_nonlinear_residual <- function(X, U, shrinkage = 0, eigenmin = 1e-3, scale = FALSE) {
  prepared <- .ld_prepare_matrix(X, "X", scale = scale)
  X_scaled <- prepared$X
  p <- ncol(X_scaled)
  q <- p - ncol(U)
  if (q <= 0L) {
    return(matrix(0, p, p))
  }
  B <- qr.Q(qr(U), complete = TRUE)[, (ncol(U) + 1L):p, drop = FALSE]
  Z <- CppMatrix::matrixMultiply(X_scaled, B)
  S_sample <- CppMatrix::matrixMultiply(Z, Z, transA = TRUE) / prepared$n_eff
  S_sample <- .ld_symmetrize(S_sample)
  S_nl <- .ld_nonlinear_qis(Z, eigenmin = eigenmin)
  shrinkage <- min(max(shrinkage[1L], 0), 1)
  S_B <- (1 - shrinkage) * S_sample + shrinkage * S_nl
  S_B <- .ld_fspd(.ld_symmetrize(S_B), eigenmin = eigenmin)
  .ld_symmetrize(CppMatrix::matrixMultiply(CppMatrix::matrixMultiply(B, S_B), B, transB = TRUE))
}
