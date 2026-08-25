#' Joint multi-source covariance transfer learning
#'
#' Compute target moments once and solve one joint multi-target shrinkage
#' problem using a list of source covariance summaries. Let `T` be the target
#' covariance and `D_s = S_s - T` for source covariance `S_s`. The method computes
#' \deqn{A_{sr}=\langle D_s,D_r\rangle_F,\qquad b_s=\widehat V_T}
#' and minimizes
#' \deqn{\alpha^T A\alpha-2b^T\alpha}
#' subject to `alpha >= 0` and `sum(alpha) <= 1`. The result is
#' \deqn{(1-\sum_s\alpha_s)T+\sum_s\alpha_sS_s.}
#'
#' The large-matrix pass is implemented in C++ without constructing the
#' `D_s` matrices. Only the small source-by-source Gram matrix is returned to
#' R, where all faces of the simplex are solved with base linear algebra. The
#' method therefore has no tuning grid or target-data split. Face enumeration
#' is intended for a small source library.
#'
#' Source finite-sample noise and all source-source interactions enter the
#' observed Gram matrix directly. Exact URE statements require independent
#' target and source samples with known population centering; sample centering
#' is a plug-in rule. URE is unbiased for every fixed coefficient vector under
#' those assumptions. `selected_ure` is the criterion at its data-selected
#' minimizer and is not an unbiased post-selection risk estimate.
#'
#' Source fourth moments are not used by the joint QP: the observed Gram already
#' contains finite-source sampling noise.
#'
#' @param X_target Target individual-level numeric matrix with observations in
#'   rows.
#' @param sources Named list of at least two source summaries. Each element is
#'   either an `ldtl_source_moments` object returned by [source_moments()] or
#'   [source_moments_plink()], or a source covariance
#'   matrix.
#' @param center If `TRUE`, subtract target sample column means. Exact
#'   finite-sample URE statements require a known population mean, corresponding
#'   to `center = FALSE` after external centering.
#'
#' @return An object of class `multisource_cov_tl` containing the joint covariance,
#'   source coefficients, target coefficient, Gram matrix, the URE linear term,
#'   selected URE value, active sources, and numerical solver diagnostics.
#' @export
multisource_cov_tl <- function(X_target, sources, center = TRUE) {
  target <- .ld_tl_target_moments(X_target, center = center)
  target_covariance <- target$covariance
  target_variance <- target$variance
  p <- nrow(target_covariance)
  resolved_sources <- .ld_multi_source_resolve(sources, p)
  source_names <- resolved_sources$source_names
  source_covariances <- resolved_sources$covariances
  n_sources <- length(source_covariances)

  gram <- cpp_multi_source_gram(target_covariance, source_covariances)
  dimnames(gram) <- list(source_names, source_names)
  linear <- rep(target_variance, n_sources)
  names(linear) <- source_names
  solution <- .ld_multi_source_qp(gram, linear)
  source_weights <- solution$coefficients
  names(source_weights) <- source_names
  target_weight <- 1 - sum(source_weights)

  covariance <- cpp_multi_source_combine(
    target_covariance,
    source_covariances,
    source_weights
  )
  selected_ure <- target_variance +
    sum(source_weights * (gram %*% source_weights)) -
    2 * sum(linear * source_weights)

  structure(
    list(
      covariance = covariance,
      source_weights = source_weights,
      target_weight = target_weight,
      active_sources = source_names[source_weights > solution$tolerance],
      gram = gram,
      linear = linear,
      target_variance = target_variance,
      selected_ure = selected_ure,
      objective = solution$objective,
      solver = solution$solver,
      numerical_tolerance = solution$tolerance
    ),
    class = "multisource_cov_tl"
  )
}

#' Multi-source covariance-path eigenspace transfer learning
#'
#' Learn one convex source direction jointly from all supplied source
#' covariances, then select the amount of transfer along that direction using
#' held-out target reconstruction scores. In each target training fold, let
#' `T` be the fold covariance, `D_s = S_s - T`, and
#' `A_sr = <D_s, D_r>_F`. The source composition is
#' \deqn{\widehat\pi=\arg\min_{\pi\geq0,\ 1^T\pi=1}\pi^TA\pi.}
#' For direct covariance weights `w` on the supplied grid, the eigenspace path
#' is fitted to
#' \deqn{C(w)=(1-w)T+w\sum_s\widehat\pi_sS_s.}
#'
#' The source summaries are resolved once and reused in every target fold.
#' Target covariance and source composition are re-estimated independently in
#' every training fold. CovTL is a separate estimator and is not fitted or
#' returned by this function.
#'
#' `method = "one_se"` returns the most transferred candidate whose paired
#' held-out target-score deficit is within one standard error of the empirical
#' winner. `method = "min"` returns the candidate with minimum held-out
#' reconstruction risk, equivalently maximum captured target score.
#'
#' @param X_target Target individual-level numeric matrix with observations in
#'   rows.
#' @param sources Named list of at least two `ldtl_source_moments` objects or
#'   source covariance matrices.
#' @param rank Eigenspace size. A number strictly between zero and one is the
#'   cumulative explained-variance threshold computed once from the full target
#'   data. A positive integer is a fixed rank. The resulting integer rank is
#'   shared by every target fold and path candidate.
#' @param folds Number of target validation folds.
#' @param transfer_weight_grid Strictly increasing direct covariance weights in
#'   `[0, 1]`, including zero and one.
#' @param center If `TRUE`, center each training and validation fold at the
#'   corresponding training-fold target mean.
#' @param fold_id Optional fold labels. Supply a shared vector to reuse the
#'   exact target folds across methods.
#' @param method Path selector. The default paired one-standard-error rule
#'   favors more transfer among competitive candidates; `"min"` selects the
#'   minimum reconstruction-risk candidate.
#' @param standard_error_multiplier Nonnegative multiplier for `method =
#'   "one_se"`. The default is one.
#' @param eigen_solver Eigensolver used after the integer rank is fixed.
#'   `"full"` uses the existing CppMatrix full eigendecomposition. `"rspectra"`
#'   independently activates RSpectra and computes only the leading `rank`
#'   eigenpairs. A fractional `rank` still requires one full-target
#'   eigendecomposition to determine the common integer cutoff.
#'
#' @return An object inheriting from `eigen_tl`. It contains the selected
#'   covariance and eigenspace, full and fold-specific source compositions,
#'   and complete held-out scores.
#' @export
multisource_eigen_tl <- function(
    X_target,
    sources,
    rank,
    folds = 5L,
    transfer_weight_grid = c(
      0, .05, .10, .20, .35, .50, .65, .80, .90, .95, 1
    ),
    center = TRUE,
    fold_id = NULL,
    method = c("one_se", "min"),
    standard_error_multiplier = 1,
    eigen_solver = c("full", "rspectra")
) {
  method <- match.arg(method)
  eigen_solver <- match.arg(eigen_solver)
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
  rank_spec <- .ld_tl_rank_spec(rank, p)
  resolved_sources <- .ld_multi_source_resolve(sources, p)
  source_names <- resolved_sources$source_names
  source_covariances <- resolved_sources$covariances
  n_sources <- length(source_covariances)
  path_weights <- .ld_path_transfer_weight_grid(transfer_weight_grid)
  fold_id <- .ld_path_fold_id(n_target, folds, fold_id)
  fold_levels <- sort(unique(fold_id))
  n_folds <- length(fold_levels)
  n_candidates <- length(path_weights)

  if (isTRUE(center)) {
    full_mean <- colMeans(X_target)
    X_full <- sweep(X_target, 2L, full_mean, "-")
  } else {
    X_full <- X_target
  }
  target_full <- .ld_tl_target_moments(X_full, center = FALSE)
  target_eig <- NULL
  if (rank_spec$type == "fixed") {
    effective_rank <- rank_spec$value
  } else {
    target_eig <- .ld_eigen(target_full$covariance)
    effective_rank <- .ld_tl_effective_rank(
      rank_spec,
      target_eig$values,
      p
    )
  }

  observation_scores <- matrix(NA_real_, n_target, n_candidates)
  fold_source_weights <- matrix(
    NA_real_,
    nrow = n_folds,
    ncol = n_sources,
    dimnames = list(as.character(fold_levels), source_names)
  )
  fold_source_distances <- numeric(n_folds)
  fold_grams <- vector("list", n_folds)
  fold_effective_ranks <- rep.int(effective_rank, n_folds)

  for (fold_index in seq_along(fold_levels)) {
    validation_indices <- which(fold_id == fold_levels[fold_index])
    training_indices <- which(fold_id != fold_levels[fold_index])
    if (!length(training_indices) || !length(validation_indices)) {
      stop(
        "Every target fold must have nonempty training and validation sets.",
        call. = FALSE
      )
    }
    X_fit <- X_target[training_indices, , drop = FALSE]
    X_validation <- X_target[validation_indices, , drop = FALSE]
    if (isTRUE(center)) {
      training_mean <- colMeans(X_fit)
      X_fit <- sweep(X_fit, 2L, training_mean, "-")
      X_validation <- sweep(X_validation, 2L, training_mean, "-")
    }
    target_fit <- .ld_tl_target_moments(X_fit, center = FALSE)
    fold_target_eig <- .ld_eigen_leading(
      target_fit$covariance,
      effective_rank,
      eigen_solver
    )
    gram <- cpp_multi_source_gram(
      target_fit$covariance,
      source_covariances
    )
    dimnames(gram) <- list(source_names, source_names)
    direction <- .ld_multi_source_composition_qp(gram)
    source_covariance <- cpp_multi_source_combine(
      target_fit$covariance,
      source_covariances,
      direction$coefficients
    )
    fold_source_weights[fold_index, ] <- direction$coefficients
    fold_source_distances[fold_index] <- direction$objective
    fold_grams[[fold_index]] <- gram

    for (candidate_index in seq_along(path_weights)) {
      weight <- path_weights[candidate_index]
      candidate_covariance <- .ld_symmetrize(
        (1 - weight) * target_fit$covariance +
          weight * source_covariance
      )
      candidate_eig <- if (candidate_index == 1L) {
        fold_target_eig
      } else {
        .ld_eigen_leading(
          candidate_covariance,
          effective_rank,
          eigen_solver
        )
      }
      candidate_vectors <-
        candidate_eig$vectors
      validation_projection <- .ld_matrix_multiply(
        X_validation,
        candidate_vectors
      )
      observation_scores[validation_indices, candidate_index] <-
        rowSums(validation_projection^2)
    }
  }

  selection <- .ld_eigen_path_selection(
    observation_scores = observation_scores,
    path_weights = path_weights,
    selector = method,
    standard_error_multiplier = standard_error_multiplier
  )

  full_gram <- cpp_multi_source_gram(
    target_full$covariance,
    source_covariances
  )
  dimnames(full_gram) <- list(source_names, source_names)
  full_direction <- .ld_multi_source_composition_qp(full_gram)
  full_source_covariance <- cpp_multi_source_combine(
    target_full$covariance,
    source_covariances,
    full_direction$coefficients
  )
  selected_weight <- path_weights[selection$selected_index]
  covariance <- .ld_symmetrize(
    (1 - selected_weight) * target_full$covariance +
      selected_weight * full_source_covariance
  )
  fit_eig <- if (selection$selected_index == 1L && !is.null(target_eig)) {
    list(
      values = target_eig$values[seq_len(effective_rank)],
      vectors = target_eig$vectors[, seq_len(effective_rank), drop = FALSE]
    )
  } else {
    .ld_eigen_leading(covariance, effective_rank, eigen_solver)
  }
  vectors <- fit_eig$vectors
  projector <- .ld_symmetrize(
    .ld_matrix_multiply(vectors, vectors, transB = TRUE)
  )

  selected_score <- selection$mean_scores[selection$selected_index]
  target_score <- selection$mean_scores[1L]
  reconstruction_gain <- selected_score - target_score

  structure(
    list(
      covariance = covariance,
      vectors = vectors,
      projector = projector,
      eigenvalues = fit_eig$values,
      effective_rank = effective_rank,
      fold_effective_ranks = fold_effective_ranks,
      selected_index = selection$selected_index,
      selected_weight = selected_weight,
      best_index = selection$best_index,
      best_weight = path_weights[selection$best_index],
      selected_score = selected_score,
      target_score = target_score,
      reconstruction_gain = reconstruction_gain,
      positive_reconstruction_gain = max(reconstruction_gain, 0),
      path_weights = path_weights,
      mean_scores = selection$mean_scores,
      score_gaps = selection$score_gaps,
      paired_standard_errors = selection$paired_standard_errors,
      competitive_set = which(selection$competitive),
      observation_scores = observation_scores,
      source_weights = stats::setNames(
        full_direction$coefficients,
        source_names
      ),
      fold_source_weights = fold_source_weights,
      source_distance_squared = full_direction$objective,
      fold_source_distances = fold_source_distances,
      gram = full_gram,
      fold_grams = fold_grams
    ),
    class = c("multisource_eigen_tl", "eigen_tl")
  )
}

.ld_multi_source_resolve <- function(sources, p) {
  if (!is.list(sources) || length(sources) < 2L) {
    stop(
      "sources must be a list containing at least two source summaries.",
      call. = FALSE
    )
  }
  source_names <- names(sources)
  if (is.null(source_names)) {
    source_names <- paste0("source_", seq_along(sources))
  }
  if (anyNA(source_names) || any(!nzchar(source_names)) ||
      anyDuplicated(source_names)) {
    stop("sources must have unique, nonempty source names.", call. = FALSE)
  }
  names(sources) <- source_names
  n_sources <- length(sources)
  if (n_sources > 15L) {
    stop(
      "At most 15 sources are supported by simplex-face enumeration.",
      call. = FALSE
    )
  }

  covariances <- vector("list", n_sources)
  names(covariances) <- source_names

  for (source_index in seq_along(sources)) {
    if (inherits(sources[[source_index]], "cov_tl") ||
        inherits(sources[[source_index]], "eigen_tl")) {
      stop(
        paste0(
          "sources must contain source moment objects or covariance matrices, ",
          "not fitted transfer-learning objects."
        ),
        call. = FALSE
      )
    }
    resolved <- .ld_tl_resolve_source(sources[[source_index]], p)
    covariances[[source_index]] <- resolved$covariance
  }

  list(
    source_names = source_names,
    covariances = covariances
  )
}

.ld_multi_source_composition_qp <- function(A) {
  n_sources <- nrow(A)
  tolerance <- 100 * .Machine$double.eps * max(1, n_sources)
  candidates <- list()

  for (mask in seq_len(2^n_sources - 1L)) {
    active <- which(as.logical(intToBits(mask)[seq_len(n_sources)]))
    A_active <- A[active, active, drop = FALSE]
    kkt <- rbind(
      cbind(A_active, rep(1, length(active))),
      c(rep(1, length(active)), 0)
    )
    boundary <- .ld_symmetric_linear_solve(
      kkt,
      c(rep(0, length(active)), 1),
      tolerance
    )
    if (!is.null(boundary)) {
      coefficients <- numeric(n_sources)
      coefficients[active] <- boundary[seq_along(active)]
      if (all(coefficients >= -tolerance) &&
          abs(sum(coefficients) - 1) <= tolerance) {
        candidates[[length(candidates) + 1L]] <- coefficients
      }
    }
  }

  if (!length(candidates)) {
    stop("No feasible multi-source composition was found.", call. = FALSE)
  }
  candidate_matrix <- do.call(rbind, candidates)
  candidate_matrix[abs(candidate_matrix) <= tolerance] <- 0
  candidate_matrix[candidate_matrix < 0] <- 0
  candidate_matrix <- candidate_matrix / rowSums(candidate_matrix)
  objectives <- apply(candidate_matrix, 1L, function(coefficients) {
    sum(coefficients * (A %*% coefficients))
  })
  best <- which.min(objectives)

  list(
    coefficients = candidate_matrix[best, ],
    objective = objectives[best],
    tolerance = tolerance,
    solver = "exact_simplex_face_enumeration"
  )
}

.ld_path_transfer_weight_grid <- function(transfer_weight_grid) {
  path_weights <- as.numeric(transfer_weight_grid)
  if (length(path_weights) < 2L || any(!is.finite(path_weights)) ||
      any(path_weights < 0) || any(path_weights > 1) ||
      is.unsorted(path_weights, strictly = TRUE) ||
      path_weights[1L] != 0 || path_weights[length(path_weights)] != 1) {
    stop(
      paste0(
        "transfer_weight_grid must be strictly increasing in [0, 1] ",
        "and include endpoints 0 and 1."
      ),
      call. = FALSE
    )
  }
  path_weights
}

.ld_eigen_path_selection <- function(
    observation_scores,
    path_weights,
    selector,
    standard_error_multiplier
) {
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
  ) / sqrt(nrow(observation_scores))
  paired_standard_errors[best_index] <- 0

  if (selector == "min") {
    selected_index <- best_index
    competitive <- seq_along(path_weights) == best_index
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

  list(
    selected_index = selected_index,
    best_index = best_index,
    mean_scores = mean_scores,
    score_gaps = score_gaps,
    paired_standard_errors = paired_standard_errors,
    competitive = competitive
  )
}

.ld_multi_source_qp <- function(A, b) {
  n_sources <- length(b)
  tolerance <- 100 * .Machine$double.eps * max(1, n_sources)
  candidates <- list(rep(0, n_sources))

  for (mask in seq_len(2^n_sources - 1L)) {
    active <- which(as.logical(intToBits(mask)[seq_len(n_sources)]))
    A_active <- A[active, active, drop = FALSE]
    b_active <- b[active]

    interior <- .ld_symmetric_linear_solve(A_active, b_active, tolerance)
    if (!is.null(interior)) {
      alpha <- numeric(n_sources)
      alpha[active] <- interior
      if (all(alpha >= -tolerance) && sum(alpha) <= 1 + tolerance) {
        candidates[[length(candidates) + 1L]] <- alpha
      }
    }

    kkt <- rbind(
      cbind(A_active, rep(1, length(active))),
      c(rep(1, length(active)), 0)
    )
    boundary <- .ld_symmetric_linear_solve(
      kkt,
      c(b_active, 1),
      tolerance
    )
    if (!is.null(boundary)) {
      alpha <- numeric(n_sources)
      alpha[active] <- boundary[seq_along(active)]
      if (all(alpha >= -tolerance) && abs(sum(alpha) - 1) <= tolerance) {
        candidates[[length(candidates) + 1L]] <- alpha
      }
    }
  }

  candidate_matrix <- do.call(rbind, candidates)
  candidate_matrix[abs(candidate_matrix) <= tolerance] <- 0
  candidate_matrix[candidate_matrix < 0] <- 0
  row_sums <- rowSums(candidate_matrix)
  boundary_rows <- row_sums > 1
  candidate_matrix[boundary_rows, ] <-
    candidate_matrix[boundary_rows, , drop = FALSE] / row_sums[boundary_rows]
  objectives <- apply(candidate_matrix, 1L, function(alpha) {
    sum(alpha * (A %*% alpha)) - 2 * sum(b * alpha)
  })
  best <- which.min(objectives)

  list(
    coefficients = candidate_matrix[best, ],
    objective = objectives[best],
    tolerance = tolerance,
    solver = "exact_simplex_face_enumeration"
  )
}

.ld_symmetric_linear_solve <- function(A, b, tolerance) {
  eig <- eigen(A, symmetric = TRUE)
  solve_tolerance <- tolerance * max(1, max(abs(A)), max(abs(b)))
  keep <- abs(eig$values) > solve_tolerance
  if (!any(keep)) {
    if (sqrt(sum(b^2)) <= solve_tolerance) {
      return(rep(0, length(b)))
    }
    return(NULL)
  }
  vectors <- eig$vectors[, keep, drop = FALSE]
  solution <- vectors %*% (crossprod(vectors, b) / eig$values[keep])
  solution <- as.numeric(solution)
  residual <- A %*% solution - b
  if (sqrt(sum(residual^2)) > solve_tolerance * max(1, sqrt(sum(b^2)))) {
    return(NULL)
  }
  solution
}
