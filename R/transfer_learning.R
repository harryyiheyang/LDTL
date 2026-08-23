#' Compute reusable source moments from a PLINK additive raw file
#'
#' Read a PLINK 2 `--export A` file in C++ blocks without materializing the
#' individual-by-variant matrix in R. The first pass estimates variant means
#' and root-mean-square scales after mean imputation; the second pass computes
#' the scalar fourth-moment sum needed by finite-source transfer learning.
#' When `R_source` is omitted, the second pass also accumulates the standardized
#' covariance with BLAS. When `R_source` is supplied, it is reused and the
#' additional scan remains \eqn{O(n_0p)}.
#'
#' Missing dosages encoded as `NA` or `.` are imputed to their variant mean.
#' Standardized second moments use normalization `1/n`, so the returned
#' covariance has unit diagonal up to numerical precision.
#'
#' @param raw_file Path to a PLINK 2 additive raw file produced by
#'   `--export A`.
#' @param R_source Optional precomputed source correlation matrix formed from
#'   the same individuals, variants, allele coding, and independent scaling.
#' @param block_size Number of samples parsed per C++/BLAS block. The default
#'   uses about `2048 * p * 8` bytes for the genotype block.
#' @param n_threads Number of OpenMP threads used by supported preprocessing
#'   steps. Zero uses the C++ runtime default. BLAS threading is controlled by
#'   the linked BLAS implementation.
#'
#' @return An object of class `ldtl_source_moments`, with the source covariance,
#'   variance, fourth-moment, sample-size, and dimension fields consumed by
#'   [cov_tl()] and [eigen_tl()]. File-scan diagnostics, variant IDs, means,
#'   scales,
#'   missing counts, and estimated C++ working bytes are also returned.
#' @export
source_moments_plink <- function(
    raw_file,
    R_source = NULL,
    block_size = 2048L,
    n_threads = 0L
) {
  if (!is.character(raw_file) || length(raw_file) != 1L ||
      is.na(raw_file) || !nzchar(raw_file)) {
    stop("raw_file must be one nonempty file path.", call. = FALSE)
  }
  raw_file <- normalizePath(
    path.expand(raw_file),
    winslash = "/",
    mustWork = TRUE
  )
  if (length(block_size) != 1L || !is.finite(block_size) ||
      block_size < 1 || block_size != round(block_size)) {
    stop("block_size must be a positive integer.", call. = FALSE)
  }
  if (length(n_threads) != 1L || !is.finite(n_threads) ||
      n_threads < 0 || n_threads != round(n_threads)) {
    stop("n_threads must be zero or a positive integer.", call. = FALSE)
  }

  compute_covariance <- is.null(R_source)
  native <- cpp_source_moments_raw(
    raw_file,
    block_size = as.integer(block_size),
    n_threads = as.integer(n_threads),
    compute_covariance = compute_covariance
  )
  if (compute_covariance) {
    R_source <- native$covariance
    variance_raw <- native$variance_raw
  } else {
    R_source <- .ld_tl_source_covariance(R_source, native$p, "R_source")
    n_double <- as.double(native$n)
    variance_raw <-
      (native$fourth_sum - n_double * sum(R_source^2)) /
      (n_double * (n_double - 1))
  }

  structure(
    list(
      R_source = R_source,
      variance = max(variance_raw, 0),
      variance_raw = variance_raw,
      fourth_sum = native$fourth_sum,
      n_source = native$n,
      p = native$p,
      mean = native$mean,
      scale = native$scale,
      observed = native$observed,
      missing_calls = native$missing_calls,
      variants = native$variants,
      exact_known_mean_moments = FALSE,
      source_format = "PLINK2 --export A",
      passes = native$passes,
      block_size = native$block_size,
      threads_used = native$threads_used,
      working_bytes = native$working_bytes,
      covariance_computed = compute_covariance
    ),
    class = "ldtl_source_moments"
  )
}

#' Compute reusable source covariance and fourth-moment summaries
#'
#' Scan source individuals in C++ blocks and return the source empirical second
#' moment together with the scalar fourth-moment variance estimate needed by
#' finite-source stabilized transfer. The returned object is reusable across
#' target cohorts and does not retain `X_source`.
#'
#' If `R_source` is supplied, it is reused and the C++ scan only computes
#' \eqn{\sum_i \|X_{0i}\|_2^4}. This reduces the additional source pass from
#' covariance cost to \eqn{O(n_0p)}. The supplied matrix must have been formed
#' from exactly the same individuals and preprocessing as `X_source`, with
#' normalization \eqn{1/n_0}.
#'
#' @param X_source Source individual-level numeric or integer matrix with
#'   individuals in rows. Missing genotypes must be imputed before calling.
#' @param R_source Optional precomputed source second-moment matrix. If omitted,
#'   it is computed in the same C++ block scan.
#' @param center If `TRUE`, subtract source sample column means. This is a
#'   practical plug-in moment rule; the exact displayed fourth-moment identity
#'   assumes the population mean is known before the scan.
#' @param block_size Number of individuals per C++ block. `NULL` chooses an
#'   automatic value using a 512 MiB working-set target, bounded between 2,048
#'   and 65,536 individuals.
#' @param n_threads Number of OpenMP threads. Zero uses the C++ runtime default.
#'
#' @return An object of class `ldtl_source_moments` containing `R_source`,
#'   `variance`, `fourth_sum`, sample size, preprocessing metadata, and native
#'   execution diagnostics.
#' @export
source_moments <- function(
    X_source,
    R_source = NULL,
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

  compute_covariance <- is.null(R_source)
  if (!compute_covariance) {
    R_source <- .ld_tl_source_covariance(R_source, p, "R_source")
  }
  native <- cpp_source_moments(
    X_source,
    center = isTRUE(center),
    block_size = block_size,
    n_threads = n_threads,
    compute_covariance = compute_covariance
  )
  if (compute_covariance) {
    R_source <- .ld_symmetrize(native$covariance)
    variance_raw <- native$variance_raw
  } else {
    n_double <- as.double(n)
    variance_raw <-
      (native$fourth_sum - n_double * sum(R_source^2)) /
      (n_double * (n_double - 1))
  }

  structure(
    list(
      R_source = R_source,
      variance = max(variance_raw, 0),
      variance_raw = variance_raw,
      fourth_sum = native$fourth_sum,
      n_source = n,
      p = p,
      mean = native$mean,
      exact_known_mean_moments = !isTRUE(center),
      block_size = native$block_size,
      threads_used = native$threads_used,
      openmp = native$openmp,
      covariance_computed = compute_covariance
    ),
    class = "ldtl_source_moments"
  )
}

#' Tuning-free covariance transfer learning
#'
#' Pool a target sample covariance with a source covariance using a closed-form
#' URE weight. With an `ldtl_source_moments` object, source sampling noise is
#' included through
#' \deqn{\widehat\lambda = \widehat V_1 /
#' \max\{\widehat V_1 + \widehat V_0,
#'        \|S_0-S_1\|_F^2\}.}
#' If `source` is only a covariance matrix, source fourth-moment information is
#' unavailable and the summary-only denominator \eqn{\|S_0-S_1\|_F^2} is used.
#'
#' @param X_target Target individual-level matrix with observations in rows.
#' @param source Either an object returned by [source_moments()] or a
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
      target_variance = target$variance,
      distance_squared = distance_squared,
      noise_floor = noise_floor,
      mismatch_squared_raw = mismatch_squared_raw,
      mismatch_squared = mismatch_squared,
      denominator = denominator
    ),
    class = "cov_tl"
  )
}

#' Fold-adjusted covariance-path eigenspace transfer learning
#'
#' Fit the source-transfer path proposed by Teacher B. Each grid point is an
#' effective fraction of the available source sample. Its covariance pooling
#' weight is recomputed for every target training fold, and the selected path
#' point is refitted using all target observations.
#'
#' `method = "one_se"` returns the most transferred candidate whose paired
#' held-out target-score deficit is within one standard error of the empirical
#' winner. `method = "min"` returns the candidate with minimum held-out
#' reconstruction risk, equivalently maximum captured target score.
#'
#' @param X_target Target individual-level data matrix.
#' @param source An `ldtl_source_moments` object or a source covariance matrix.
#' @param rank Number of leading target eigenvectors.
#' @param n_source Source sample size. Required when `source` is a matrix and
#'   checked against the stored sample size when source moments are supplied.
#' @param folds Number of target validation folds.
#' @param source_fraction_grid Increasing effective-source fractions in
#'   `[0, 1]`, including zero and one.
#' @param center If `TRUE`, center each training and validation fold at the
#'   corresponding training-fold target mean.
#' @param fold_id Optional fold labels. Supply a shared vector when comparing
#'   sources so all source candidates use identical target folds.
#' @param method Path selector. The default paired one-standard-error rule
#'   favors more transfer among competitive candidates; `"min"` selects the
#'   minimum reconstruction-risk candidate.
#' @param standard_error_multiplier Nonnegative multiplier for `method =
#'   "one_se"`. The default is one.
#'
#' @return An object of class `eigen_tl` containing the selected full-data
#'   covariance, eigenspace, complete held-out path scores, paired uncertainty,
#'   and source-transfer diagnostics.
#' @export
eigen_tl <- function(
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
    method = c("one_se", "min"),
    standard_error_multiplier = 1
) {
  method <- match.arg(method)
  .ld_eigen_tl_fit(
    X_target = X_target,
    source = source,
    rank = rank,
    n_source = n_source,
    folds = folds,
    source_fraction_grid = source_fraction_grid,
    center = center,
    fold_id = fold_id,
    selector = method,
    standard_error_multiplier = standard_error_multiplier
  )
}

.ld_eigen_tl_fit <- function(
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

  if (selector == "min") {
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
  selected_score <- mean_scores[selected_index]
  target_score <- mean_scores[1L]
  reconstruction_gain <- selected_score - target_score

  structure(
    list(
      covariance = covariance,
      vectors = vectors,
      projector = projector,
      eigenvalues = fit_eig$values[seq_len(rank)],
      selected_index = selected_index,
      selected_zeta = zeta[selected_index],
      selected_weight = selected_weight,
      pooled_endpoint_selected = selected_index == length(zeta),
      best_index = best_index,
      best_zeta = zeta[best_index],
      best_weight = full_weights[best_index],
      selected_score = selected_score,
      target_score = target_score,
      reconstruction_gain = reconstruction_gain,
      positive_reconstruction_gain = max(reconstruction_gain, 0),
      full_data_weights = full_weights,
      fold_weights = fold_weights,
      mean_scores = mean_scores,
      score_gaps = score_gaps,
      paired_standard_errors = paired_standard_errors,
      competitive_set = which(competitive),
      observation_scores = observation_scores
    ),
    class = "eigen_tl"
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
  if (inherits(source, "ldtl_source_moments")) {
    required <- c("R_source", "variance", "n_source", "p")
    if (!all(required %in% names(source))) {
      stop("The source moment object is incomplete.", call. = FALSE)
    }
    covariance <- .ld_tl_source_covariance(
      source$R_source,
      p,
      "source$R_source"
    )
    if (length(source$variance) != 1L || !is.finite(source$variance) ||
        source$variance < 0) {
      stop("source$variance must be a nonnegative finite scalar.", call. = FALSE)
    }
    if (length(source$n_source) != 1L || !is.finite(source$n_source) ||
        source$n_source < 2) {
      stop("source$n_source must be at least 2.", call. = FALSE)
    }
    return(list(
      covariance = covariance,
      variance = as.double(source$variance),
      n = as.double(source$n_source),
      has_variance = TRUE,
      exact_known_mean_moments = isTRUE(source$exact_known_mean_moments)
    ))
  }

  if (inherits(source, "ld_source_moments")) {
    required <- c("covariance", "variance", "n", "p")
    if (!all(required %in% names(source))) {
      stop("The legacy source moment object is incomplete.", call. = FALSE)
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
