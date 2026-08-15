.ld_as_square_matrix <- function(x, name = "x") {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (nrow(x) != ncol(x)) {
    stop(name, " must be a square matrix.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

.ld_as_matrix <- function(x, name = "x") {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  storage.mode(x) <- "double"
  if (anyNA(x)) {
    stop(name, " must not contain missing values.", call. = FALSE)
  }
  x
}

.ld_symmetrize <- function(x) {
  (x + t(x)) / 2
}

.ld_matrix_multiply <- function(A, B, transA = FALSE, transB = FALSE) {
  tryCatch(
    CppMatrix::matrixMultiply(
      A,
      B,
      transA = transA,
      transB = transB
    ),
    error = function(e) {
      if (isTRUE(transA)) {
        A <- t(A)
      }
      if (isTRUE(transB)) {
        B <- t(B)
      }
      A %*% B
    }
  )
}

.ld_center_columns <- function(x) {
  tryCatch(
    CppMatrix::matrixScale(
      x,
      center = TRUE,
      standardized = FALSE,
      robust = FALSE
    ),
    error = function(e) sweep(x, 2L, colMeans(x), "-")
  )
}

.ld_eigen <- function(x) {
  x <- .ld_symmetrize(x)
  eig <- tryCatch(
    CppMatrix::matrixEigen(x),
    error = function(e) eigen(x, symmetric = TRUE)
  )
  values <- as.numeric(eig$values)
  vectors <- eig$vectors
  ord <- order(values, decreasing = TRUE)
  list(values = values[ord], vectors = vectors[, ord, drop = FALSE])
}

.ld_min_eigen <- function(x) {
  x <- .ld_symmetrize(x)
  values <- tryCatch(
    CppMatrix::matrixEigen(x)$values,
    error = function(e) eigen(x, symmetric = TRUE, only.values = TRUE)$values
  )
  values <- as.numeric(values)
  if (!length(values) || all(is.na(values))) {
    return(NA_real_)
  }
  min(values, na.rm = TRUE)
}

.ld_clean_correlation <- function(x, name = "x") {
  x <- .ld_as_square_matrix(x, name)
  x <- .ld_symmetrize(x)
  out <- stats::cov2cor(x)
  out[is.na(out)] <- 0
  out <- .ld_symmetrize(out)
  diag(out) <- 1
  out
}

.ld_matrix_cor <- function(x, scale = FALSE, name = "X") {
  if (isTRUE(scale)) {
    x <- .ld_as_matrix(x, name)
    out <- tryCatch(
      CppMatrix::matrixCor(x),
      error = function(e) stats::cor(x, use = "pairwise.complete.obs")
    )
    return(.ld_clean_correlation(out, "correlation matrix"))
  }

  x <- .ld_as_matrix(x, name)
  if (nrow(x) < 2L || ncol(x) < 1L) {
    stop(name, " must have at least 2 rows and 1 column.", call. = FALSE)
  }
  out <- tryCatch(
    .ld_matrix_multiply(x, x, transA = TRUE) / (nrow(x) - 1L),
    error = function(e) crossprod(x) / (nrow(x) - 1L)
  )
  .ld_clean_correlation(out, "correlation matrix")
}

.ld_scale_matrix <- function(x, name = "X") {
  x <- .ld_as_matrix(x, name)
  if (nrow(x) < 3L || ncol(x) < 2L) {
    stop(name, " must have at least 3 rows and 2 columns.", call. = FALSE)
  }
  center <- colMeans(x)
  x <- sweep(x, 2L, center, "-")
  scale <- sqrt(colSums(x^2) / (nrow(x) - 1L))
  scale[!is.finite(scale) | scale <= 0] <- 1
  x <- sweep(x, 2L, scale, "/")
  list(X = x, n = nrow(x), n_eff = nrow(x) - 1L)
}

.ld_prepare_matrix <- function(x, name = "X", scale = FALSE) {
  if (isTRUE(scale)) {
    return(.ld_scale_matrix(x, name))
  }
  x <- .ld_as_matrix(x, name)
  if (nrow(x) < 2L || ncol(x) < 1L) {
    stop(name, " must have at least 2 rows and 1 column.", call. = FALSE)
  }
  list(X = x, n = nrow(x), n_eff = nrow(x) - 1L)
}

.ld_resolve_input <- function(S = NULL, X = NULL, n = NULL, name = "S", scale = FALSE) {
  if (is.null(S) && is.null(X)) {
    stop("Provide either ", name, " or X.", call. = FALSE)
  }
  if (!is.null(S)) {
    S <- .ld_as_square_matrix(S, name)
    if (!is.null(X)) {
      prepared <- .ld_prepare_matrix(X, "X", scale = scale)
      X <- prepared$X
      if (is.null(n)) {
        n <- prepared$n
      }
    }
    return(list(S = S, X = X, n = n))
  }
  prepared <- .ld_prepare_matrix(X, "X", scale = scale)
  list(
    S = .ld_matrix_cor(prepared$X, scale = FALSE),
    X = prepared$X,
    n = prepared$n
  )
}

.ld_validate_n <- function(n) {
  if (is.null(n) || length(n) != 1L || !is.finite(n) || n <= 1) {
    stop("n must be a finite sample size greater than 1.", call. = FALSE)
  }
  as.numeric(n)
}

.ld_theory_threshold <- function(p, n) {
  n <- .ld_validate_n(n)
  2 * sqrt(log(max(p, 2L)) / n)
}

.ld_resolve_threshold <- function(p, n, lambda) {
  if (!is.null(lambda)) {
    return(lambda)
  }
  if (!is.null(n)) {
    return(.ld_theory_threshold(p, n))
  }
  stop("Provide n to compute the default lambda or supply lambda.", call. = FALSE)
}

.ld_theory_poet_threshold <- function(p, n) {
  n <- .ld_validate_n(n)
  2 * max(sqrt(log(max(p, 2L)) / n), 1 / sqrt(p))
}

.ld_resolve_poet_threshold <- function(p, n, lambda) {
  if (!is.null(lambda)) {
    return(lambda)
  }
  .ld_theory_poet_threshold(p, n)
}

.ld_theory_bandwidth <- function(p, n, alpha) {
  n <- .ld_validate_n(n)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("alpha must be a positive finite scalar.", call. = FALSE)
  }
  max_k <- .ld_max_bandwidth(p)
  k <- ceiling(n^(1 / (2 * alpha + 1)))
  min(max(as.integer(k), 0L), max_k)
}

.ld_max_bandwidth <- function(p) {
  max(as.integer(floor(p / 2)), 0L)
}

.ld_resolve_bandwidth <- function(
    p,
    n,
    K,
    alpha
) {
  max_k <- .ld_max_bandwidth(p)
  if (!is.null(K)) {
    return(min(max(as.integer(round(K[1L])), 0L), max_k))
  }
  if (is.null(n)) {
    stop("Provide n to compute the default K or supply K.", call. = FALSE)
  }
  .ld_theory_bandwidth(p, n, alpha)
}

.ld_band_weight <- function(p, k) {
  d <- 0:(p - 1L)
  stats::toeplitz(as.numeric(d <= k))
}

.ld_taper_weight <- function(p, k, k_gap) {
  if (length(k_gap) != 1L || !is.finite(k_gap) || k_gap < 1) {
    stop("k_gap must be a finite scalar greater than or equal to 1.", call. = FALSE)
  }
  d <- 0:(p - 1L)
  upper_raw <- k_gap * k
  w <- numeric(p)
  w[d <= k] <- 1
  upper <- min(ceiling(upper_raw), p - 1L)
  if (upper > k) {
    idx <- which(d > k & d <= upper)
    denom <- upper_raw - k
    if (denom > 0) {
      w[idx] <- pmax((upper_raw - d[idx]) / denom, 0)
    }
  }
  stats::toeplitz(w)
}

.ld_fspd <- function(x, eig_min = 0, method = "FSopt", verbose = FALSE) {
  x <- .ld_symmetrize(x)
  p <- nrow(x)
  method <- match.arg(method, c("FSopt", "Sopt", "Fopt", "Max", "Infty"))

  eig <- .ld_eigen(x)
  values <- eig$values
  lmin <- min(values)
  lmax <- max(values)
  lmean <- mean(values)

  info <- list(
    method = method,
    operated = FALSE,
    alpha = 1,
    target = NA_real_,
    min_before = lmin,
    min_after = lmin
  )

  if (!is.finite(lmin) || lmin >= eig_min) {
    attr(x, "fspd") <- info
    return(x)
  }

  if (identical(method, "Infty")) {
    out <- x + diag(eig_min - lmin, p)
    out <- .ld_symmetrize(out)
    info$operated <- TRUE
    info$alpha <- NA_real_
    info$target <- Inf
    info$min_after <- .ld_min_eigen(out)
    attr(out, "fspd") <- info
    return(out)
  }

  target_s <- (lmin + lmax) / 2
  target_f <- if (abs(lmean - lmin) > .Machine$double.eps) {
    lmean + stats::var(values) * (p - 1) / p / (lmean - lmin)
  } else {
    NA_real_
  }

  target <- switch(
    method,
    Sopt = target_s,
    Fopt = target_f,
    Max = lmax,
    FSopt = max(c(target_s, target_f), na.rm = TRUE)
  )

  if (!is.finite(target) || target <= lmin + .Machine$double.eps) {
    out <- x + diag(eig_min - lmin, p)
    out <- .ld_symmetrize(out)
    info$operated <- TRUE
    info$alpha <- NA_real_
    info$target <- target
    info$min_after <- .ld_min_eigen(out)
    attr(out, "fspd") <- info
    return(out)
  }

  alpha <- 1 - (eig_min - lmin) / (target - lmin)
  alpha <- min(max(alpha, 0), 1)
  out <- alpha * x
  diag(out) <- diag(out) + (1 - alpha) * target
  out <- .ld_symmetrize(out)

  min_after <- .ld_min_eigen(out)
  if (is.finite(min_after) && min_after < eig_min) {
    diag(out) <- diag(out) + (eig_min - min_after)
    out <- .ld_symmetrize(out)
    min_after <- .ld_min_eigen(out)
  }

  if (isTRUE(verbose)) {
    message(
      "FSPD ", method, ": alpha = ", signif(alpha, 4),
      ", target = ", signif(target, 4),
      ", min eigenvalue ", signif(lmin, 4), " -> ", signif(min_after, 4)
    )
  }

  info$operated <- TRUE
  info$alpha <- alpha
  info$target <- target
  info$min_after <- min_after
  attr(out, "fspd") <- info
  out
}

.ld_fix_residual_diag <- function(x) {
  x <- .ld_symmetrize(x)
  d <- diag(x)
  pos_min <- suppressWarnings(min(d[d > 0]))
  if (is.finite(pos_min)) {
    d[d < 0] <- pos_min
  } else {
    d[d < 0] <- 0
  }
  diag(x) <- d
  x
}

.ld_factor_count_from_values <- function(d, k_min, k_max, cutoff_method = "D.ratio") {
  p <- length(d)
  if (p <= 2L) {
    return(1L)
  }

  k_min <- max(as.integer(k_min), 2L)
  k_max <- min(as.integer(k_max), p - 2L)
  if (k_min > k_max) {
    return(min(max(1L, k_max), p - 1L))
  }

  z <- rep(NA_real_, p)
  if (identical(cutoff_method, "D.ratio")) {
    j_lo <- max(k_min, 2L)
    j_hi <- min(k_max, p - 1L)
    if (j_lo <= j_hi) {
      j <- seq.int(j_lo, j_hi)
      den <- d[j] - d[j + 1L]
      num <- d[j - 1L] - d[j]
      z[j - 1L] <- ifelse(den != 0, num / den, NA_real_)
    }
  } else if (identical(cutoff_method, "ratio")) {
    j_lo <- max(k_min, 2L)
    j_hi <- min(k_max, p)
    if (j_lo <= j_hi) {
      j <- seq.int(j_lo, j_hi)
      z[j - 1L] <- ifelse(d[j] != 0, d[j - 1L] / d[j], NA_real_)
    }
  } else {
    stop("Unknown cutoff_method.", call. = FALSE)
  }

  idx <- suppressWarnings(which.max(z))
  pck <- if (length(idx) == 0L || is.na(idx)) max(k_min, 1L) else idx + 1L
  min(max(pck, 1L), p - 1L)
}

.ld_poet_components <- function(S, n, cutoff_method, k_min, k_max) {
  S <- .ld_as_square_matrix(S, "S")
  p <- nrow(S)
  S[is.na(S)] <- 0
  diag(S) <- 1
  S[abs(S) < 1e-4] <- 0
  S <- .ld_symmetrize(S)

  eig <- .ld_eigen(S)
  d <- eig$values
  U <- eig$vectors

  pck <- .ld_factor_count_from_values(d, k_min, k_max, cutoff_method)

  Uk <- U[, seq_len(pck), drop = FALSE]
  dk <- d[seq_len(pck)]
  denom <- p - pck - pck * p / n
  hatc <- if (is.finite(denom) && denom > 0) {
    (sum(diag(S)) - sum(dk)) / denom
  } else {
    0
  }
  dk <- pmax(.ld_mcp(dk, hatc * p / n, a = 3), 0)
  P <- tcrossprod(t(t(Uk) * dk), Uk)
  P <- .ld_symmetrize(P)
  E <- .ld_fix_residual_diag(S - P)

  list(S = S, P = P, E = E, U = Uk, factors = pck, eigenvalues = d)
}
