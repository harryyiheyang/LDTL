#' Joint multi-source covariance transfer learning
#'
#' Combine source-specific [cov_tl()] fits by solving one joint multi-target
#' shrinkage problem. The individual `lambda` values are not averaged. Let `T`
#' be the common target covariance and `D_s = S_s - T` for source covariance
#' `S_s`. The method computes
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
#' @param fits Named list of at least two source-specific `cov_tl` fits. Every
#'   fit must use exactly the same target covariance and preprocessing.
#'
#' @return An object of class `multi_source_tl` containing the joint covariance,
#'   source coefficients, target coefficient, Gram matrix, the URE linear term,
#'   selected URE value, active sources, numerical solver diagnostics, and the
#'   original fits.
#' @export
multi_source_tl <- function(fits) {
  if (!is.list(fits) || length(fits) < 2L) {
    stop("fits must be a list containing at least two fitted objects.", call. = FALSE)
  }
  source_names <- names(fits)
  if (is.null(source_names)) {
    source_names <- paste0("source_", seq_along(fits))
  }
  if (anyNA(source_names) || any(!nzchar(source_names)) ||
      anyDuplicated(source_names)) {
    stop("fits must have unique, nonempty source names.", call. = FALSE)
  }
  names(fits) <- source_names

  if (!all(vapply(fits, inherits, logical(1), what = "cov_tl"))) {
    stop(
      paste0(
        "fits must contain only cov_tl objects. Joint multi-source EigenTL ",
        "requires cross-source held-out path terms that are not stored in ",
        "separate eigen_tl fits."
      ),
      call. = FALSE
    )
  }
  n_sources <- length(fits)
  if (n_sources > 15L) {
    stop(
      "multi_source_tl supports at most 15 sources because it enumerates simplex faces.",
      call. = FALSE
    )
  }

  target_covariance <- fits[[1L]]$target_covariance
  if (!is.matrix(target_covariance) || nrow(target_covariance) < 1L ||
      nrow(target_covariance) != ncol(target_covariance) ||
      any(!is.finite(target_covariance))) {
    stop("Every fit must contain a finite square target_covariance.", call. = FALSE)
  }
  p <- nrow(target_covariance)
  covariance_tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(abs(target_covariance)))
  source_covariances <- vector("list", n_sources)
  names(source_covariances) <- source_names
  target_variances <- numeric(n_sources)
  source_variances <- rep(NA_real_, n_sources)

  for (source_index in seq_along(fits)) {
    source <- fits[[source_index]]$source_covariance
    target <- fits[[source_index]]$target_covariance
    target_variance <- fits[[source_index]]$target_variance
    source_variance <- fits[[source_index]]$source_variance
    if (!is.matrix(source) || !identical(dim(source), c(p, p)) ||
        any(!is.finite(source))) {
      stop(
        "Every fit must contain a compatible finite source_covariance.",
        call. = FALSE
      )
    }
    if (!is.matrix(target) || !identical(dim(target), c(p, p)) ||
        any(!is.finite(target)) ||
        max(abs(target - target_covariance)) > covariance_tolerance) {
      stop("All fits must use the same target covariance.", call. = FALSE)
    }
    if (!is.numeric(target_variance) || length(target_variance) != 1L ||
        !is.finite(target_variance) || target_variance < 0) {
      stop(
        "Every fit must contain a nonnegative finite target_variance.",
        call. = FALSE
      )
    }
    if (!is.null(source_variance) &&
        (!is.numeric(source_variance) || length(source_variance) != 1L ||
         (!is.na(source_variance) &&
          (!is.finite(source_variance) || source_variance < 0)))) {
      stop(
        "Each source_variance must be nonnegative and finite or NA.",
        call. = FALSE
      )
    }
    source_covariances[[source_index]] <- source
    target_variances[source_index] <- target_variance
    if (!is.null(source_variance)) {
      source_variances[source_index] <- source_variance
    }
  }

  variance_tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(target_variances))
  if (max(abs(target_variances - target_variances[1L])) > variance_tolerance) {
    stop("All fits must use the same target variance estimate.", call. = FALSE)
  }
  target_variance <- target_variances[1L]

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
      weights = source_weights,
      target_weight = target_weight,
      active_sources = source_names[source_weights > solution$tolerance],
      gram = gram,
      linear = linear,
      target_variance = target_variance,
      source_variance = stats::setNames(source_variances, source_names),
      selected_ure = selected_ure,
      objective = solution$objective,
      solver = solution$solver,
      numerical_tolerance = solution$tolerance,
      aggregation = "joint_multi_target_ure",
      fit_family = "cov_tl",
      source_names = source_names,
      target_covariance = target_covariance,
      source_covariances = source_covariances,
      fits = fits
    ),
    class = "multi_source_tl"
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
