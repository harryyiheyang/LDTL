#' Tuning-free covariance transfer learning
#'
#' Pool a target sample covariance with a source covariance using the closed-form
#' unbiased-risk-estimation (URE) weight. The source and target samples need not
#' overlap, and source individual-level data are not required.
#'
#' The estimator is
#' \deqn{\widehat\Sigma=(1-\widehat\lambda)S_1+
#' \widehat\lambda S_0,}
#' where
#' \deqn{\widehat\lambda=\left[\widehat V_1 /
#' \|S_0-S_1\|_F^2\right]_{[0,1]}.}
#' Here `S1 = crossprod(X_target) / n` and `V1` is the matrix sample-variance
#' estimator. Thus no cross-validation or tuning grid is used.
#' Large matrix products and covariance construction use `CppMatrix` routines,
#' with base R retained only as a fallback and for unsupported scalar
#' reductions.
#'
#' @param X_target Target individual-level data matrix, with observations in
#'   rows. For the exact finite-sample URE identity, its rows should be
#'   independent and already centered at the known population mean.
#' @param S_source Source covariance matrix on the same variables, ordering,
#'   scale, and covariance normalization as the target covariance.
#' @param center If `TRUE`, subtract target sample column means. This is useful
#'   in practice but turns the exact known-mean URE into a plug-in rule. The
#'   default is `TRUE`.
#'
#' @return An object of class `cov_tl`, represented by a list with the pooled
#'   `covariance`, analytic `lambda` and `alpha`, target covariance, target
#'   variance estimate, observed source-target squared distance, and URE
#'   quadratic coefficients.
#' @export
cov_tl <- function(X_target, S_source, center = TRUE) {
  target <- .ld_tl_target_moments(X_target, center = center)
  S_source <- .ld_tl_source_covariance(S_source, ncol(target$X))

  distance_squared <- .ld_tl_frobenius_distance_squared(
    S_source,
    target$covariance
  )
  weight <- .ld_tl_closed_form_weight(
    numerator = target$variance,
    denominator = distance_squared
  )

  pooled <- .ld_symmetrize(
    (1 - weight$lambda) * target$covariance +
      weight$lambda * S_source
  )

  structure(
    list(
      covariance = pooled,
      lambda = weight$lambda,
      alpha = weight$alpha,
      target_covariance = target$covariance,
      source_covariance = S_source,
      target_variance = target$variance,
      distance_squared = distance_squared,
      ure = c(
        constant = target$variance,
        linear = -2 * target$variance,
        quadratic = distance_squared
      ),
      n_target = target$n,
      center = isTRUE(center),
      exact_known_mean_ure = !isTRUE(center)
    ),
    class = "cov_tl"
  )
}

#' Tuning-free eigenspace transfer learning
#'
#' Choose the source pooling rate by minimizing the first-order target
#' eigenspace (spectral-projector tangent) risk. The analytic rule is
#' \deqn{\widehat\lambda_{ES}=\left[\widehat T_1 /
#' \widehat Q_{ES}\right]_{[0,1]},}
#' where `T1` is the target tangent-noise estimate and `Q_ES` is the observed
#' source-target distance after applying the target spectral-projector
#' derivative. The final eigenspace is computed from the pooled covariance.
#'
#' If `S_pilot` is supplied from data independent of `X_target`, the empirical
#' tangent-risk quadratic is conditionally unbiased for that fixed derivative.
#' If it is omitted, the target sample eigensystem is reused, giving the usual
#' full-sample first-order plug-in estimator. Neither version uses
#' cross-validation or a tuning grid.
#' The covariance, score, cross-block, projector, and eigendecomposition
#' computations use `CppMatrix` routines whenever available.
#'
#' @param X_target Target individual-level data matrix, with observations in
#'   rows. For the fixed-derivative URE identity, its rows should be independent
#'   and centered at the known population mean.
#' @param S_source Source covariance matrix on the same variables, ordering,
#'   scale, and covariance normalization as the target covariance.
#' @param rank Number of leading target eigenvectors to estimate. Must be
#'   between 1 and `ncol(X_target) - 1`.
#' @param center If `TRUE`, subtract target sample column means. This gives a
#'   practical plug-in rule rather than the exact known-mean identity. The
#'   default is `TRUE`.
#' @param S_pilot Optional covariance matrix used only to estimate the
#'   spectral-projector derivative. It should be independent of `X_target` for
#'   the conditional unbiased-risk interpretation.
#' @param eigengap_tol Optional positive tolerance below which the pilot
#'   eigengap is treated as unidentified. By default a scale-adjusted numerical
#'   tolerance is used.
#'
#' @return An object of class `eigspac_tl`, represented by a list containing the
#'   pooled covariance, its leading eigenvectors and projector, analytic
#'   `lambda` and `alpha`, tangent-risk quantities, and pilot eigengap.
#' @export
eigspac_tl <- function(
    X_target,
    S_source,
    rank,
    center = TRUE,
    S_pilot = NULL,
    eigengap_tol = NULL
) {
  target <- .ld_tl_target_moments(X_target, center = center)
  p <- ncol(target$X)
  S_source <- .ld_tl_source_covariance(S_source, p)
  rank <- .ld_tl_rank(rank, p)

  if (is.null(S_pilot)) {
    S_derivative <- target$covariance
    pilot_independent <- FALSE
  } else {
    S_derivative <- .ld_tl_source_covariance(S_pilot, p, "S_pilot")
    pilot_independent <- TRUE
  }

  pilot_eig <- .ld_eigen(S_derivative)
  scale <- max(1, abs(pilot_eig$values))
  if (is.null(eigengap_tol)) {
    eigengap_tol <- sqrt(.Machine$double.eps) * scale
  }
  if (length(eigengap_tol) != 1L || !is.finite(eigengap_tol) ||
      eigengap_tol <= 0) {
    stop("eigengap_tol must be a positive finite scalar.", call. = FALSE)
  }

  eigengap <- pilot_eig$values[rank] - pilot_eig$values[rank + 1L]
  if (!is.finite(eigengap) || eigengap <= eigengap_tol) {
    stop(
      "The pilot eigengap at rank is too small for the tangent-risk rule.",
      call. = FALSE
    )
  }

  leading <- seq_len(rank)
  trailing <- seq.int(rank + 1L, p)
  U_leading <- pilot_eig$vectors[, leading, drop = FALSE]
  U_trailing <- pilot_eig$vectors[, trailing, drop = FALSE]
  gaps <- outer(
    pilot_eig$values[leading],
    pilot_eig$values[trailing],
    "-"
  )

  scores_leading <- .ld_matrix_multiply(target$X, U_leading)
  scores_trailing <- .ld_matrix_multiply(target$X, U_trailing)
  target_cross_block <- .ld_matrix_multiply(
    scores_leading,
    scores_trailing,
    transA = TRUE
  ) / target$n

  cross_score_squares <- .ld_matrix_multiply(
    scores_leading^2,
    scores_trailing^2,
    transA = TRUE
  )
  centered_cross_squares <-
    cross_score_squares - target$n * target_cross_block^2
  centered_cross_squares <- pmax(centered_cross_squares, 0)
  tangent_sum <- sum(centered_cross_squares / gaps^2)
  tangent_variance <- 2 * tangent_sum / (target$n * (target$n - 1L))
  tangent_variance <- max(tangent_variance, 0)

  source_trailing <- .ld_matrix_multiply(S_source, U_trailing)
  source_cross_block <- .ld_matrix_multiply(
    U_leading,
    source_trailing,
    transA = TRUE
  )
  delta_cross_block <- source_cross_block - target_cross_block
  tangent_distance_squared <- 2 * sum((delta_cross_block / gaps)^2)
  weight <- .ld_tl_closed_form_weight(
    numerator = tangent_variance,
    denominator = tangent_distance_squared
  )

  pooled <- .ld_symmetrize(
    (1 - weight$lambda) * target$covariance +
      weight$lambda * S_source
  )
  pooled_eig <- .ld_eigen(pooled)
  vectors <- pooled_eig$vectors[, leading, drop = FALSE]
  projector <- .ld_matrix_multiply(vectors, vectors, transB = TRUE)

  structure(
    list(
      covariance = pooled,
      vectors = vectors,
      projector = .ld_symmetrize(projector),
      eigenvalues = pooled_eig$values[leading],
      rank = rank,
      lambda = weight$lambda,
      alpha = weight$alpha,
      target_covariance = target$covariance,
      source_covariance = S_source,
      tangent_variance = tangent_variance,
      tangent_distance_squared = tangent_distance_squared,
      tangent_ure = c(
        constant = tangent_variance,
        linear = -2 * tangent_variance,
        quadratic = tangent_distance_squared
      ),
      pilot_covariance = S_derivative,
      pilot_eigengap = eigengap,
      pilot_supplied = pilot_independent,
      n_target = target$n,
      center = isTRUE(center),
      conditional_ure_if_pilot_independent =
        pilot_independent && !isTRUE(center)
    ),
    class = "eigspac_tl"
  )
}

.ld_tl_target_moments <- function(X_target, center = FALSE) {
  X_target <- .ld_as_matrix(X_target, "X_target")
  if (nrow(X_target) < 2L || ncol(X_target) < 1L) {
    stop(
      "X_target must have at least 2 rows and 1 column.",
      call. = FALSE
    )
  }
  if (any(!is.finite(X_target))) {
    stop("X_target must contain only finite values.", call. = FALSE)
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("center must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(center)) {
    X_target <- .ld_center_columns(X_target)
  }

  n <- nrow(X_target)
  S_target <- .ld_symmetrize(
    .ld_matrix_multiply(X_target, X_target, transA = TRUE) / n
  )
  row_norm_squared <- rowSums(X_target^2)
  variance <-
    (sum(row_norm_squared^2) - n * sum(S_target^2)) /
    (n * (n - 1L))

  list(
    X = X_target,
    n = n,
    covariance = S_target,
    variance = max(variance, 0)
  )
}

.ld_tl_frobenius_distance_squared <- function(A, B) {
  if (!all(dim(A) == dim(B))) {
    stop("A and B must have the same dimensions.", call. = FALSE)
  }
  distance <- norm(A - B, type = "F")
  distance^2
}

.ld_tl_source_covariance <- function(S_source, p, name = "S_source") {
  S_source <- .ld_as_square_matrix(S_source, name)
  if (!all(dim(S_source) == c(p, p))) {
    stop(name, " must have the same number of variables as X_target.", call. = FALSE)
  }
  if (any(!is.finite(S_source))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  .ld_symmetrize(S_source)
}

.ld_tl_closed_form_weight <- function(numerator, denominator) {
  tolerance <- .Machine$double.eps * max(1, numerator, denominator)
  if (denominator <= tolerance) {
    lambda <- if (numerator > tolerance) 1 else 0
  } else {
    lambda <- min(max(numerator / denominator, 0), 1)
  }
  alpha <- if (lambda >= 1) Inf else lambda / (1 - lambda)
  list(lambda = lambda, alpha = alpha)
}

.ld_tl_rank <- function(rank, p) {
  if (length(rank) != 1L || !is.finite(rank) || rank != round(rank)) {
    stop("rank must be a single integer.", call. = FALSE)
  }
  rank <- as.integer(rank)
  if (rank < 1L || rank >= p) {
    stop("rank must be between 1 and ncol(X_target) - 1.", call. = FALSE)
  }
  rank
}
