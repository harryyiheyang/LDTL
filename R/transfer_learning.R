#' Compute reusable source covariance and fourth-moment summaries
#'
#' Scan source individuals in C++ blocks and return the source empirical second
#' moment together with the scalar fourth-moment variance estimate needed by
#' finite-source stabilized transfer. The returned object is reusable across
#' target cohorts and does not retain `X_source`.
#'
#' If `S_source` is supplied, it is reused and the C++ scan only computes
#' \eqn{\sum_i \|X_{0i}\|_2^4}. This reduces the additional source pass from
#' covariance cost to \eqn{O(n_0p)}. The supplied matrix must have been formed
#' from exactly the same individuals and preprocessing as `X_source`, with
#' normalization \eqn{1/n_0}.
#'
#' @param X_source Source individual-level numeric or integer matrix with
#'   individuals in rows. Missing genotypes must be imputed before calling.
#' @param S_source Optional precomputed source second-moment matrix. If omitted,
#'   it is computed in the same C++ block scan.
#' @param center If `TRUE`, subtract source sample column means. This is a
#'   practical plug-in moment rule; the exact displayed fourth-moment identity
#'   assumes the population mean is known before the scan.
#' @param block_size Number of individuals per C++ block. `NULL` chooses an
#'   automatic value using a 512 MiB working-set target, bounded between 2,048
#'   and 65,536 individuals.
#' @param n_threads Number of OpenMP threads. Zero uses the C++ runtime default.
#'
#' @return An object of class `ld_source_moments` containing `covariance`,
#'   `variance`, `fourth_sum`, sample size, preprocessing metadata, and native
#'   execution diagnostics.
#' @export
source_moments_method <- function(
    X_source,
    S_source = NULL,
    center = FALSE,
    block_size = NULL,
    n_threads = 0L
) {
  if (!is.matrix(X_source) || !(is.numeric(X_source) || is.integer(X_source))) {
    stop(
      "X_source must be a numeric or integer matrix with individuals in rows.",
      call. = FALSE
    )
  }
  n <- nrow(X_source)
  p <- ncol(X_source)
  if (n < 2L || p < 1L) {
    stop("X_source must have at least 2 rows and 1 column.", call. = FALSE)
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("center must be TRUE or FALSE.", call. = FALSE)
  }
  block_size <- .ld_tl_source_block_size(n, p, block_size)
  if (length(n_threads) != 1L || !is.finite(n_threads) ||
      n_threads < 0 || n_threads != round(n_threads)) {
    stop("n_threads must be zero or a positive integer.", call. = FALSE)
  }
  n_threads <- as.integer(n_threads)

  compute_covariance <- is.null(S_source)
  if (!compute_covariance) {
    S_source <- .ld_tl_source_covariance(S_source, p)
  }
  native <- cpp_source_moments(
    X_source,
    center = isTRUE(center),
    block_size = block_size,
    n_threads = n_threads,
    compute_covariance = compute_covariance
  )
  if (compute_covariance) {
    S_source <- .ld_symmetrize(native$covariance)
    variance_raw <- native$variance_raw
  } else {
    n_double <- as.double(n)
    variance_raw <-
      (native$fourth_sum - n_double * sum(S_source^2)) /
      (n_double * (n_double - 1))
  }

  structure(
    list(
      covariance = S_source,
      variance = max(variance_raw, 0),
      variance_raw = variance_raw,
      fourth_sum = native$fourth_sum,
      n = n,
      p = p,
      mean = native$mean,
      center = isTRUE(center),
      exact_known_mean_moments = !isTRUE(center),
      normalization = "1/n",
      block_size = native$block_size,
      threads_used = native$threads_used,
      openmp = native$openmp,
      covariance_computed = compute_covariance
    ),
    class = "ld_source_moments"
  )
}

#' Tuning-free covariance transfer learning
#'
#' Pool a target sample covariance with a source covariance using a closed-form
#' URE weight. With a `ld_source_moments` object, source sampling noise is
#' included through
#' \deqn{\widehat\lambda = \widehat V_1 /
#' \max\{\widehat V_1 + \widehat V_0,
#'        \|S_0-S_1\|_F^2\}.}
#' If `source` is only a covariance matrix, source fourth-moment information is
#' unavailable and the summary-only denominator \eqn{\|S_0-S_1\|_F^2} is used.
#'
#' @param X_target Target individual-level matrix with observations in rows.
#' @param source Either an object returned by [source_moments_method()] or a
#'   source covariance matrix.
#' @param center If `TRUE`, subtract target sample column means.
#'
#' @return An object of class `cov_tl` with the pooled covariance, selected
#'   weight, source and target noise estimates, and mismatch diagnostics.
#' @export
cov_tl <- function(X_target, source, center = TRUE) {
  target <- .ld_tl_target_moments(X_target, center = center)
  resolved <- .ld_tl_resolve_source(source, ncol(target$X))
  S_source <- resolved$covariance
  distance_squared <- .ld_tl_frobenius_distance_squared(
    S_source,
    target$covariance
  )

  if (resolved$has_variance) {
    noise_floor <- target$variance + resolved$variance
    denominator <- max(noise_floor, distance_squared)
    mismatch_squared_raw <- distance_squared - noise_floor
    mismatch_squared <- max(mismatch_squared_raw, 0)
    weight_method <- "finite_source"
  } else {
    noise_floor <- NA_real_
    denominator <- distance_squared
    mismatch_squared_raw <- NA_real_
    mismatch_squared <- NA_real_
    weight_method <- "summary_only_ure"
  }
  weight <- .ld_tl_closed_form_weight(
    numerator = target$variance,
    denominator = denominator
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
      weight_method = weight_method,
      target_covariance = target$covariance,
      source_covariance = S_source,
      target_variance = target$variance,
      source_variance = resolved$variance,
      distance_squared = distance_squared,
      noise_floor = noise_floor,
      mismatch_squared_raw = mismatch_squared_raw,
      mismatch_squared = mismatch_squared,
      denominator = denominator,
      n_target = target$n,
      n_source = resolved$n,
      center = isTRUE(center),
      source_moments_supplied = resolved$has_variance,
      exact_known_mean_ure = !isTRUE(center),
      exact_known_mean_moments =
        !isTRUE(center) && isTRUE(resolved$exact_known_mean_moments)
    ),
    class = "cov_tl"
  )
}

#' Legacy summary-only covariance transfer learning
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
#' @return An object of class `cov_tl1`, represented by a list with the pooled
#'   `covariance`, analytic `lambda` and `alpha`, target covariance, target
#'   variance estimate, observed source-target squared distance, and URE
#'   quadratic coefficients.
#' @keywords internal
cov_tl1 <- function(X_target, S_source, center = TRUE) {
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
    class = "cov_tl1"
  )
}

#' Legacy summary-only eigenspace transfer learning
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
#' @return An object of class `eigspac_tl1`, represented by a list containing the
#'   pooled covariance, its leading eigenvectors and projector, analytic
#'   `lambda` and `alpha`, tangent-risk quantities, and pilot eigengap.
#' @keywords internal
eigspac_tl1 <- function(
    X_target,
    S_source,
    rank,
    center = TRUE,
    S_pilot = NULL,
    eigengap_tol = NULL
) {
  components <- .ld_tl_eigenspace_components(
    X_target = X_target,
    S_source = S_source,
    rank = rank,
    center = center,
    S_pilot = S_pilot,
    eigengap_tol = eigengap_tol
  )
  weight <- .ld_tl_closed_form_weight(
    numerator = components$tangent_variance,
    denominator = components$tangent_distance_squared
  )
  .ld_tl_eigenspace_result(
    components,
    weight,
    extra = list(
      weight_method = "tangent_ure",
      tangent_ure = c(
        constant = components$tangent_variance,
        linear = -2 * components$tangent_variance,
        quadratic = components$tangent_distance_squared
      )
    ),
    class_name = "eigspac_tl1"
  )
}

#' Tuning-free eigenspace transfer learning
#'
#' Choose the source pooling rate by minimizing the first-order target
#' eigenspace risk. When source noise information is available, the denominator
#' includes a source tangent-noise term before truncating the mismatch at zero.
#'
#' A reusable scalar covariance fourth moment does not identify the exact
#' source tangent variance for an arbitrary target eigenspace. When
#' `source_tangent_variance` is omitted, this function therefore uses the
#' explicit conservative plug-in bound
#' \deqn{\widehat T_0^{bound}=\widehat V_0/\widehat\gamma^2,}
#' where `gamma` is the pilot eigengap. An exact target-specific source tangent
#' estimate may instead be supplied by the caller. With a covariance-matrix
#' source and no supplied tangent variance, the summary-only tangent URE rule is
#' used.
#'
#' @param X_target Target individual-level data matrix.
#' @param source An `ld_source_moments` object from [source_moments_method()] or
#'   a source covariance matrix.
#' @param rank Number of leading target eigenvectors.
#' @param center Whether to center target columns.
#' @param S_pilot Optional independent covariance used for the projector
#'   derivative.
#' @param eigengap_tol Optional positive eigengap tolerance.
#' @param source_tangent_variance Optional nonnegative target-specific estimate
#'   of source tangent noise. If omitted, the Frobenius-noise/eigengap bound is
#'   used.
#'
#' @return An object of class `eigspac_tl`.
#' @export
eigspac_tl <- function(
    X_target,
    source,
    rank,
    center = TRUE,
    S_pilot = NULL,
    eigengap_tol = NULL,
    source_tangent_variance = NULL
) {
  resolved <- .ld_tl_resolve_source(source, ncol(X_target))
  components <- .ld_tl_eigenspace_components(
    X_target = X_target,
    S_source = resolved$covariance,
    rank = rank,
    center = center,
    S_pilot = S_pilot,
    eigengap_tol = eigengap_tol
  )

  if (is.null(source_tangent_variance) && resolved$has_variance) {
    source_tangent_variance <-
      resolved$variance / components$pilot_eigengap^2
    source_tangent_variance_type <- "frobenius_eigengap_upper_proxy"
  } else if (!is.null(source_tangent_variance)) {
    if (length(source_tangent_variance) != 1L ||
        !is.finite(source_tangent_variance) ||
        source_tangent_variance < 0) {
      stop(
        "source_tangent_variance must be a nonnegative finite scalar.",
        call. = FALSE
      )
    }
    source_tangent_variance_type <- "target_specific_supplied"
  } else {
    weight <- .ld_tl_closed_form_weight(
      numerator = components$tangent_variance,
      denominator = components$tangent_distance_squared
    )
    return(.ld_tl_eigenspace_result(
      components,
      weight,
      extra = list(
        weight_method = "summary_only_tangent_ure",
        source_variance = NA_real_,
        source_tangent_variance = NA_real_,
        source_tangent_variance_type = "unavailable",
        tangent_noise_floor = NA_real_,
        tangent_mismatch_squared_raw = NA_real_,
        tangent_mismatch_squared = NA_real_,
        denominator = components$tangent_distance_squared,
        n_source = NA_real_
      ),
      class_name = "eigspac_tl"
    ))
  }

  noise_floor <- components$tangent_variance + source_tangent_variance
  denominator <- max(noise_floor, components$tangent_distance_squared)
  mismatch_squared_raw <- components$tangent_distance_squared - noise_floor
  mismatch_squared <- max(mismatch_squared_raw, 0)
  weight <- .ld_tl_closed_form_weight(
    numerator = components$tangent_variance,
    denominator = denominator
  )

  .ld_tl_eigenspace_result(
    components,
    weight,
    extra = list(
      weight_method = "finite_source_tangent",
      source_variance = resolved$variance,
      source_tangent_variance = source_tangent_variance,
      source_tangent_variance_type = source_tangent_variance_type,
      tangent_noise_floor = noise_floor,
      tangent_mismatch_squared_raw = mismatch_squared_raw,
      tangent_mismatch_squared = mismatch_squared,
      denominator = denominator,
      n_source = resolved$n
    ),
    class_name = "eigspac_tl"
  )
}

.ld_tl_eigenspace_components <- function(
    X_target,
    S_source,
    rank,
    center,
    S_pilot,
    eigengap_tol
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

  list(
    target = target,
    source_covariance = S_source,
    rank = rank,
    leading = leading,
    tangent_variance = tangent_variance,
    tangent_distance_squared = tangent_distance_squared,
    pilot_covariance = S_derivative,
    pilot_eigengap = eigengap,
    pilot_supplied = pilot_independent,
    center = isTRUE(center)
  )
}

.ld_tl_eigenspace_result <- function(
    components,
    weight,
    extra,
    class_name
) {
  pooled <- .ld_symmetrize(
    (1 - weight$lambda) * components$target$covariance +
      weight$lambda * components$source_covariance
  )
  pooled_eig <- .ld_eigen(pooled)
  vectors <- pooled_eig$vectors[, components$leading, drop = FALSE]
  projector <- .ld_matrix_multiply(vectors, vectors, transB = TRUE)
  common <- list(
    covariance = pooled,
    vectors = vectors,
    projector = .ld_symmetrize(projector),
    eigenvalues = pooled_eig$values[components$leading],
    rank = components$rank,
    lambda = weight$lambda,
    alpha = weight$alpha,
    target_covariance = components$target$covariance,
    source_covariance = components$source_covariance,
    tangent_variance = components$tangent_variance,
    tangent_distance_squared = components$tangent_distance_squared,
    pilot_covariance = components$pilot_covariance,
    pilot_eigengap = components$pilot_eigengap,
    pilot_supplied = components$pilot_supplied,
    n_target = components$target$n,
    center = components$center,
    conditional_ure_if_pilot_independent =
      components$pilot_supplied && !components$center
  )
  structure(c(common, extra), class = class_name)
}

#' Max-score fold-adjusted covariance-path transfer PCA
#'
#' Teacher-B covariance-path estimator using ordinary target cross-validation.
#' A grid point is an effective fraction of the available source sample, so the
#' final grid point is exact sample-size pooling in every training fold. The
#' empirical target-score maximizer is then refitted on all observations.
#'
#' This is a separate estimator from the analytic CovTL and EigenTL families.
#' It does not use source fourth moments, although an `ld_source_moments` object
#' may be supplied as a convenient source covariance/sample-size container.
#'
#' @param X_target Target individual-level data matrix.
#' @param source An `ld_source_moments` object or a source covariance matrix.
#' @param rank Number of leading target eigenvectors.
#' @param n_source Source sample size. Required when `source` is a matrix and
#'   checked against the stored sample size when source moments are supplied.
#' @param folds Number of target validation folds.
#' @param source_fraction_grid Increasing effective-source fractions in
#'   `[0, 1]`, including zero and one.
#' @param center If `TRUE`, center each training and validation fold at the
#'   corresponding training-fold target mean. The exact known-mean theory uses
#'   `center = FALSE`.
#' @param fold_id Optional integer fold labels. Supplying these is the preferred
#'   way to make comparisons exactly reproducible. Otherwise balanced labels
#'   are randomly permuted using R's current RNG state.
#'
#' @return An object of class `path_tpca_max_score` containing the selected
#'   full-data covariance/projector, complete fold scores, path weights, and
#'   fold assignment.
#' @export
path_tpca_max_score <- function(
    X_target,
    source,
    rank,
    n_source = NULL,
    folds = 5L,
    source_fraction_grid = c(
      0, .05, .10, .20, .35, .50, .65, .80, .90, .95, 1
    ),
    center = TRUE,
    fold_id = NULL
) {
  .ld_path_tpca_fit(
    X_target = X_target,
    source = source,
    rank = rank,
    n_source = n_source,
    folds = folds,
    source_fraction_grid = source_fraction_grid,
    center = center,
    fold_id = fold_id,
    selector = "max_score",
    standard_error_multiplier = 0
  )
}

#' Paired one-standard-error fold-adjusted covariance-path transfer PCA
#'
#' Teacher-B path-adaptive estimator. It builds the same fold-adjusted path as
#' [path_tpca_max_score()], but returns the most transferred candidate whose
#' paired held-out target-score deficit is within one estimated standard error
#' of the empirical winner. The selected effective source fraction is then
#' refitted using all target observations and the source covariance.
#'
#' @inheritParams path_tpca_max_score
#' @param standard_error_multiplier Nonnegative multiplier on the paired
#'   standard error. The proposed default is one.
#'
#' @return An object of class `path_tpca_one_se` containing the selected
#'   full-data covariance/projector, paired competitive set, fold scores, and
#'   refit diagnostics.
#' @export
path_tpca_one_se <- function(
    X_target,
    source,
    rank,
    n_source = NULL,
    folds = 5L,
    source_fraction_grid = c(
      0, .05, .10, .20, .35, .50, .65, .80, .90, .95, 1
    ),
    center = TRUE,
    fold_id = NULL,
    standard_error_multiplier = 1
) {
  .ld_path_tpca_fit(
    X_target = X_target,
    source = source,
    rank = rank,
    n_source = n_source,
    folds = folds,
    source_fraction_grid = source_fraction_grid,
    center = center,
    fold_id = fold_id,
    selector = "paired_one_se",
    standard_error_multiplier = standard_error_multiplier
  )
}

.ld_path_tpca_fit <- function(
    X_target,
    source,
    rank,
    n_source,
    folds,
    source_fraction_grid,
    center,
    fold_id,
    selector,
    standard_error_multiplier
) {
  X_target <- .ld_as_matrix(X_target, "X_target")
  if (nrow(X_target) < 2L || ncol(X_target) < 2L ||
      any(!is.finite(X_target))) {
    stop(
      "X_target must have at least 2 rows, 2 columns, and only finite values.",
      call. = FALSE
    )
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("center must be TRUE or FALSE.", call. = FALSE)
  }
  n_target <- nrow(X_target)
  p <- ncol(X_target)
  rank <- .ld_tl_rank(rank, p)
  source_resolved <- .ld_tl_resolve_source(source, p)
  n_source <- .ld_path_source_size(source_resolved, n_source)
  zeta <- .ld_path_source_fraction_grid(source_fraction_grid)
  fold_id <- .ld_path_fold_id(n_target, folds, fold_id)
  fold_levels <- sort(unique(fold_id))
  n_candidates <- length(zeta)
  observation_scores <- matrix(NA_real_, n_target, n_candidates)
  fold_weights <- matrix(
    NA_real_,
    nrow = length(fold_levels),
    ncol = n_candidates
  )

  for (fold_index in seq_along(fold_levels)) {
    validation_indices <- which(fold_id == fold_levels[fold_index])
    training_indices <- which(fold_id != fold_levels[fold_index])
    n_fit <- length(training_indices)
    if (n_fit < 1L || !length(validation_indices)) {
      stop("Every target fold must have nonempty training and validation sets.", call. = FALSE)
    }
    X_fit <- X_target[training_indices, , drop = FALSE]
    X_validation <- X_target[validation_indices, , drop = FALSE]
    if (isTRUE(center)) {
      training_mean <- colMeans(X_fit)
      X_fit <- sweep(X_fit, 2L, training_mean, "-")
      X_validation <- sweep(X_validation, 2L, training_mean, "-")
    }
    S_fit <- .ld_symmetrize(
      .ld_matrix_multiply(X_fit, X_fit, transA = TRUE) / n_fit
    )
    weights <- .ld_path_weights(zeta, n_fit, n_source)
    fold_weights[fold_index, ] <- weights

    for (candidate_index in seq_along(zeta)) {
      candidate_covariance <- .ld_symmetrize(
        (1 - weights[candidate_index]) * S_fit +
          weights[candidate_index] * source_resolved$covariance
      )
      candidate_eig <- .ld_eigen(candidate_covariance)
      candidate_vectors <-
        candidate_eig$vectors[, seq_len(rank), drop = FALSE]
      validation_projection <- .ld_matrix_multiply(
        X_validation,
        candidate_vectors
      )
      observation_scores[validation_indices, candidate_index] <-
        rowSums(validation_projection^2)
    }
  }

  if (any(!is.finite(observation_scores))) {
    stop("Nonfinite held-out target scores were produced.", call. = FALSE)
  }
  mean_scores <- colMeans(observation_scores)
  best_index <- which.max(mean_scores)
  score_gaps <- mean_scores[best_index] - mean_scores
  paired_differences <-
    observation_scores[, best_index] - observation_scores
  paired_standard_errors <- apply(
    paired_differences,
    2L,
    stats::sd
  ) / sqrt(n_target)
  paired_standard_errors[best_index] <- 0

  if (selector == "max_score") {
    selected_index <- best_index
    competitive <- seq_along(zeta) == best_index
  } else {
    if (length(standard_error_multiplier) != 1L ||
        !is.finite(standard_error_multiplier) ||
        standard_error_multiplier < 0) {
      stop(
        "standard_error_multiplier must be a nonnegative finite scalar.",
        call. = FALSE
      )
    }
    tolerance <- .Machine$double.eps * max(1, abs(mean_scores))
    competitive <- score_gaps <=
      standard_error_multiplier * paired_standard_errors + tolerance
    selected_index <- max(which(competitive))
  }

  if (isTRUE(center)) {
    full_mean <- colMeans(X_target)
    X_full <- sweep(X_target, 2L, full_mean, "-")
  } else {
    full_mean <- rep(0, p)
    X_full <- X_target
  }
  target_covariance <- .ld_symmetrize(
    .ld_matrix_multiply(X_full, X_full, transA = TRUE) / n_target
  )
  full_weights <- .ld_path_weights(zeta, n_target, n_source)
  selected_weight <- full_weights[selected_index]
  covariance <- .ld_symmetrize(
    (1 - selected_weight) * target_covariance +
      selected_weight * source_resolved$covariance
  )
  fit_eig <- .ld_eigen(covariance)
  vectors <- fit_eig$vectors[, seq_len(rank), drop = FALSE]
  projector <- .ld_symmetrize(
    .ld_matrix_multiply(vectors, vectors, transB = TRUE)
  )
  class_name <- if (selector == "max_score") {
    "path_tpca_max_score"
  } else {
    "path_tpca_one_se"
  }

  structure(
    list(
      covariance = covariance,
      vectors = vectors,
      projector = projector,
      eigenvalues = fit_eig$values[seq_len(rank)],
      rank = rank,
      selector = selector,
      selected_index = selected_index,
      selected_zeta = zeta[selected_index],
      selected_weight = selected_weight,
      pooled_endpoint_selected = selected_index == length(zeta),
      max_score_index = best_index,
      max_score_zeta = zeta[best_index],
      max_score_weight = full_weights[best_index],
      source_fraction_grid = zeta,
      full_data_weights = full_weights,
      fold_weights = fold_weights,
      mean_scores = mean_scores,
      score_gaps = score_gaps,
      paired_standard_errors = paired_standard_errors,
      competitive_set = which(competitive),
      observation_scores = observation_scores,
      fold_id = fold_id,
      target_covariance = target_covariance,
      source_covariance = source_resolved$covariance,
      target_mean = full_mean,
      n_target = n_target,
      n_source = n_source,
      center = isTRUE(center),
      standard_error_multiplier = if (selector == "max_score") {
        NA_real_
      } else {
        standard_error_multiplier
      },
      exact_fold_pooling_endpoint = TRUE,
      source_fourth_moment_used = FALSE
    ),
    class = class_name
  )
}

.ld_path_source_size <- function(source_resolved, n_source) {
  stored <- source_resolved$n
  if (is.null(n_source)) {
    if (!is.finite(stored)) {
      stop(
        "n_source is required when source is supplied as a covariance matrix.",
        call. = FALSE
      )
    }
    return(as.double(stored))
  }
  if (length(n_source) != 1L || !is.finite(n_source) ||
      n_source < 1 || n_source != round(n_source)) {
    stop("n_source must be a positive integer.", call. = FALSE)
  }
  if (is.finite(stored) && as.double(n_source) != as.double(stored)) {
    stop("n_source does not match the sample size stored in source.", call. = FALSE)
  }
  as.double(n_source)
}

.ld_path_source_fraction_grid <- function(source_fraction_grid) {
  zeta <- as.numeric(source_fraction_grid)
  if (length(zeta) < 2L || any(!is.finite(zeta)) ||
      any(zeta < 0) || any(zeta > 1) || is.unsorted(zeta, strictly = TRUE) ||
      zeta[1L] != 0 || zeta[length(zeta)] != 1) {
    stop(
      "source_fraction_grid must be strictly increasing in [0, 1] and include endpoints 0 and 1.",
      call. = FALSE
    )
  }
  zeta
}

.ld_path_fold_id <- function(n_target, folds, fold_id) {
  if (is.null(fold_id)) {
    if (length(folds) != 1L || !is.finite(folds) ||
        folds < 2 || folds > n_target || folds != round(folds)) {
      stop("folds must be an integer between 2 and nrow(X_target).", call. = FALSE)
    }
    folds <- as.integer(folds)
    return(sample(rep(seq_len(folds), length.out = n_target)))
  }
  if (length(fold_id) != n_target || anyNA(fold_id) ||
      any(!is.finite(fold_id))) {
    stop("fold_id must contain one finite label per target observation.", call. = FALSE)
  }
  labels <- match(fold_id, unique(fold_id))
  if (length(unique(labels)) < 2L) {
    stop("fold_id must define at least two nonempty folds.", call. = FALSE)
  }
  as.integer(labels)
}

.ld_path_weights <- function(zeta, n_target_fit, n_source) {
  effective_source <- zeta * n_source
  effective_source / (n_target_fit + effective_source)
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

.ld_tl_source_block_size <- function(n, p, block_size) {
  if (is.null(block_size)) {
    target_bytes <- 512 * 1024^2
    block_size <- floor(target_bytes / (8 * p))
    block_size <- min(65536, max(2048, block_size))
    return(as.integer(min(n, block_size)))
  }
  if (length(block_size) != 1L || !is.finite(block_size) ||
      block_size < 1 || block_size != round(block_size)) {
    stop("block_size must be NULL or a positive integer.", call. = FALSE)
  }
  as.integer(min(n, block_size))
}

.ld_tl_resolve_source <- function(source, p) {
  if (inherits(source, "ld_source_moments")) {
    required <- c("covariance", "variance", "n", "p")
    if (!all(required %in% names(source))) {
      stop("The source moment object is incomplete.", call. = FALSE)
    }
    covariance <- .ld_tl_source_covariance(
      source$covariance,
      p,
      "source$covariance"
    )
    if (length(source$variance) != 1L || !is.finite(source$variance) ||
        source$variance < 0) {
      stop("source$variance must be a nonnegative finite scalar.", call. = FALSE)
    }
    if (length(source$n) != 1L || !is.finite(source$n) || source$n < 2) {
      stop("source$n must be at least 2.", call. = FALSE)
    }
    return(list(
      covariance = covariance,
      variance = as.double(source$variance),
      n = as.double(source$n),
      has_variance = TRUE,
      exact_known_mean_moments = isTRUE(source$exact_known_mean_moments)
    ))
  }

  list(
    covariance = .ld_tl_source_covariance(source, p, "source"),
    variance = NA_real_,
    n = NA_real_,
    has_variance = FALSE,
    exact_known_mean_moments = FALSE
  )
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
