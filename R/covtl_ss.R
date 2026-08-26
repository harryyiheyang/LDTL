#' Multi-source covariance transfer from summary statistics
#'
#' Estimate shared convex weights for locus-specific source LD matrices from
#' independent null-locus summary vectors. For locus `l`, the working model is
#' \deqn{z_l \sim N(0, \Sigma_l(\alpha))}
#' \deqn{\Sigma_l(\alpha) = \alpha_0 I + \sum_{k=1}^K \alpha_k R_{lk}}
#' with nonnegative weights summing to one and a small upper bound on the
#' identity weight. The estimator minimizes the summed Gaussian negative
#' log-likelihood across loci using Fisher scoring on the capped simplex.
#'
#' `CovTL_SS()` does not modify or replace [cov_tl()]. The latter requires
#' individual-level target data, whereas `CovTL_SS()` uses null-locus summary
#' statistics only.
#'
#' Each element of `sources[[l]]` may be a covariance/LD matrix or a truncated
#' eigensystem `list(U = U, D = D)`. Eigensystems may be full rank or low rank.
#' Every reconstructed source is symmetrized and converted to a correlation
#' matrix before fitting. Source names and their order are shared across loci.
#'
#' The optimizer uses at most `max_iterations` Fisher updates. Each update
#' evaluates at most `max_line_search` step sizes. At the last permitted step,
#' a finite positive-definite candidate is accepted even when the likelihood
#' does not decrease, matching the intended bounded-search implementation.
#'
#' @param z Numeric null-locus z-score vector, or a list of such vectors. Each
#'   locus contributes one multivariate observation.
#' @param sources For one locus, a named list of source LD matrices or
#'   `list(U, D)` eigensystems. For multiple loci, a list of these named
#'   source lists in the same order as `z`.
#' @param identity_max Fixed upper bound for the identity coefficient. The
#'   default is `0.015`.
#' @param max_iterations Positive integer maximum number of Fisher updates.
#' @param max_line_search Positive integer maximum number of line-search
#'   evaluations per Fisher update.
#' @param tolerance Positive convergence tolerance for the constrained score
#'   gap and coefficient updates.
#' @param profile_tolerance Positive tolerance for the one-dimensional
#'   identity profiles used to initialize the optimizer.
#'
#' @return An object of class `CovTL_SS` containing shared weights, fitted
#'   locus covariance matrices, Gaussian negative log-likelihood, convergence
#'   diagnostics, the free-coordinate Fisher information matrix, and its
#'   minimum eigenvalue and condition number.
#' @export
CovTL_SS <- function(
    z,
    sources,
    identity_max = 0.015,
    max_iterations = 20L,
    max_line_search = 3L,
    tolerance = 1e-6,
    profile_tolerance = 1e-4
) {
  input <- .ld_ss_resolve_input(z, sources)
  z_list <- input$z
  A_list <- input$sources
  source_names <- input$source_names
  K <- length(source_names)

  identity_max <- .ld_ss_scalar(
    identity_max,
    "identity_max",
    lower = 0,
    upper = 1,
    upper_open = TRUE
  )
  max_iterations <- .ld_ss_positive_integer(
    max_iterations,
    "max_iterations"
  )
  max_line_search <- .ld_ss_positive_integer(
    max_line_search,
    "max_line_search"
  )
  tolerance <- .ld_ss_scalar(
    tolerance,
    "tolerance",
    lower = 0,
    lower_open = TRUE
  )
  profile_tolerance <- .ld_ss_scalar(
    profile_tolerance,
    "profile_tolerance",
    lower = 0,
    lower_open = TRUE
  )

  starts <- vector("list", K + 1L)
  values <- rep(Inf, K + 1L)
  for (j in seq_len(K)) {
    source_weights <- numeric(K)
    source_weights[j] <- 1
    starts[[j]] <- .ld_ss_profile_identity(
      source_weights,
      z_list,
      A_list,
      identity_max,
      profile_tolerance
    )
    values[j] <- starts[[j]]$negative_log_likelihood
  }
  starts[[K + 1L]] <- .ld_ss_profile_identity(
    rep(1 / K, K),
    z_list,
    A_list,
    identity_max,
    profile_tolerance
  )
  values[K + 1L] <- starts[[K + 1L]]$negative_log_likelihood
  if (!any(is.finite(values))) {
    stop(
      "No positive-definite source/identity initialization was found.",
      call. = FALSE
    )
  }

  start_index <- which.min(values)
  weights <- starts[[start_index]]$weights
  negative_log_likelihood <- values[start_index]
  d <- K + 1L
  converged <- FALSE
  updates <- 0L
  forced_updates <- 0L
  line_search_evaluations <- 0L
  score_gap <- Inf
  likelihood_trace <- negative_log_likelihood

  for (iteration in seq_len(max_iterations)) {
    fisher_step <- .ld_ss_full_fisher(weights, z_list, A_list)
    gradient <- fisher_step$gradient
    fisher <- fisher_step$fisher

    vertex <- numeric(d)
    source_index <- which.min(gradient[seq_len(K)])
    vertex[source_index] <- 1
    if (gradient[d] < gradient[source_index]) {
      vertex[source_index] <- 1 - identity_max
      vertex[d] <- identity_max
    }
    score_gap <- sum(gradient * (weights - vertex))
    if (score_gap <= tolerance * max(1, abs(negative_log_likelihood))) {
      converged <- TRUE
      break
    }

    proposal <- .ld_ss_capped_simplex_qp(
      fisher,
      as.numeric(fisher %*% weights - gradient),
      identity_max
    )
    if (max(abs(proposal - weights)) <= tolerance) {
      converged <- TRUE
      break
    }

    accepted <- FALSE
    for (line_index in seq_len(max_line_search)) {
      step_size <- 2^(-(line_index - 1L))
      candidate <- (1 - step_size) * weights + step_size * proposal
      candidate_nll <- .ld_ss_nll(candidate, z_list, A_list)
      line_search_evaluations <- line_search_evaluations + 1L
      if (is.finite(candidate_nll) &&
          (candidate_nll < negative_log_likelihood ||
           line_index == max_line_search)) {
        if (line_index == max_line_search &&
            candidate_nll >= negative_log_likelihood) {
          forced_updates <- forced_updates + 1L
        }
        weights <- candidate
        negative_log_likelihood <- candidate_nll
        accepted <- TRUE
        break
      }
    }
    if (!accepted) {
      stop(
        "Fisher proposal did not produce a finite positive-definite update.",
        call. = FALSE
      )
    }
    updates <- updates + 1L
    likelihood_trace <- c(likelihood_trace, negative_log_likelihood)
  }

  names(weights) <- c(source_names, "Identity")
  covariance <- lapply(A_list, .ld_ss_combine, weights = weights)
  names(covariance) <- names(z_list)
  information <- .ld_ss_free_information(weights, A_list)

  structure(
    list(
      weights = weights,
      source_weights = weights[source_names],
      identity_weight = unname(weights["Identity"]),
      covariance = covariance,
      negative_log_likelihood = negative_log_likelihood,
      likelihood_trace = likelihood_trace,
      iterations = updates,
      converged = converged,
      score_gap = score_gap,
      forced_updates = forced_updates,
      line_search_evaluations = line_search_evaluations,
      fisher_information = information$matrix,
      fisher_min_eigenvalue = information$minimum,
      fisher_condition_number = information$condition,
      source_names = source_names,
      locus_dimensions = vapply(z_list, length, integer(1L)),
      identity_max = identity_max,
      max_iterations = max_iterations,
      max_line_search = max_line_search,
      tolerance = tolerance,
      profile_tolerance = profile_tolerance,
      initialization = start_index
    ),
    class = "CovTL_SS"
  )
}

.ld_ss_resolve_input <- function(z, sources) {
  if (is.numeric(z) && is.null(dim(z))) {
    z <- list(locus_1 = z)
    sources <- list(locus_1 = sources)
  } else if (is.list(z) && length(z)) {
    if (is.null(names(z))) {
      names(z) <- paste0("locus_", seq_along(z))
    }
  } else {
    stop("z must be a numeric vector or a nonempty list of vectors.", call. = FALSE)
  }
  if (!is.list(sources) || length(sources) != length(z)) {
    stop("sources must contain one source list for every z locus.", call. = FALSE)
  }

  z_list <- vector("list", length(z))
  names(z_list) <- names(z)
  A_list <- vector("list", length(z))
  names(A_list) <- names(z)
  source_names <- NULL

  for (locus_index in seq_along(z)) {
    z_locus <- as.numeric(z[[locus_index]])
    if (length(z_locus) < 2L || any(!is.finite(z_locus))) {
      stop(
        "Every z locus must contain at least two finite values.",
        call. = FALSE
      )
    }
    z_list[[locus_index]] <- z_locus
    locus_sources <- sources[[locus_index]]
    if (!is.list(locus_sources) || !length(locus_sources)) {
      stop("Every locus must contain at least one source LD input.", call. = FALSE)
    }
    if (length(locus_sources) > 15L) {
      stop("At most 15 sources are supported.", call. = FALSE)
    }

    locus_names <- names(locus_sources)
    if (is.null(locus_names)) {
      locus_names <- paste0("source_", seq_along(locus_sources))
    }
    if (anyNA(locus_names) || any(!nzchar(locus_names)) ||
        anyDuplicated(locus_names)) {
      stop("Source names must be unique and nonempty within each locus.", call. = FALSE)
    }
    names(locus_sources) <- locus_names
    if (is.null(source_names)) {
      source_names <- locus_names
    } else {
      if (!setequal(locus_names, source_names)) {
        stop("All loci must contain the same named sources.", call. = FALSE)
      }
      locus_sources <- locus_sources[source_names]
    }

    resolved <- vector("list", length(source_names) + 1L)
    names(resolved) <- c(source_names, "Identity")
    for (source_index in seq_along(source_names)) {
      resolved[[source_index]] <- .ld_ss_resolve_source(
        locus_sources[[source_index]],
        length(z_locus),
        paste0("sources[[", locus_index, "]][[", source_index, "]]"
        )
      )
    }
    resolved[[length(resolved)]] <- diag(length(z_locus))
    A_list[[locus_index]] <- resolved
  }

  list(z = z_list, sources = A_list, source_names = source_names)
}

.ld_ss_resolve_source <- function(source, p, name) {
  if (is.list(source) && !is.null(source$U) && !is.null(source$D)) {
    U <- as.matrix(source$U)
    storage.mode(U) <- "double"
    D <- as.numeric(source$D)
    if (nrow(U) != p || ncol(U) != length(D) || !length(D) ||
        any(!is.finite(U)) || any(!is.finite(D)) || any(D < 0) ||
        !any(D > 0)) {
      stop(
        name,
        " must have finite p by r U and r nonnegative D with a positive value.",
        call. = FALSE
      )
    }
    source <- .ld_matrix_multiply(
      sweep(U, 2L, D, "*"),
      U,
      transB = TRUE
    )
  } else {
    source <- .ld_as_square_matrix(source, name)
    if (!all(dim(source) == c(p, p)) || any(!is.finite(source))) {
      stop(name, " must be a finite p by p matrix.", call. = FALSE)
    }
  }

  source <- .ld_symmetrize(source)
  diagonal <- diag(source)
  if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    stop(name, " must have a positive finite diagonal.", call. = FALSE)
  }
  scale <- sqrt(diagonal)
  source <- source / tcrossprod(scale)
  source <- .ld_symmetrize(source)
  diag(source) <- 1
  minimum <- .ld_min_eigen(source)
  tolerance <- sqrt(.Machine$double.eps) * max(1, p)
  if (!is.finite(minimum) || minimum < -tolerance) {
    stop(name, " must be positive semidefinite after cov2cor scaling.", call. = FALSE)
  }
  source
}

.ld_ss_scalar <- function(
    x,
    name,
    lower = -Inf,
    upper = Inf,
    lower_open = FALSE,
    upper_open = FALSE
) {
  if (length(x) != 1L || !is.finite(x)) {
    stop(name, " must be a finite scalar.", call. = FALSE)
  }
  invalid_lower <- if (lower_open) x <= lower else x < lower
  invalid_upper <- if (upper_open) x >= upper else x > upper
  if (invalid_lower || invalid_upper) {
    stop(name, " is outside its permitted range.", call. = FALSE)
  }
  as.numeric(x)
}

.ld_ss_positive_integer <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x != round(x) || x < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  as.integer(x)
}

.ld_ss_combine <- function(A, weights) {
  covariance <- weights[1L] * A[[1L]]
  if (length(A) > 1L) {
    for (j in 2:length(A)) {
      covariance <- covariance + weights[j] * A[[j]]
    }
  }
  .ld_symmetrize(covariance)
}

.ld_ss_nll <- function(weights, z_list, A_list) {
  value <- 0
  for (locus_index in seq_along(z_list)) {
    covariance <- .ld_ss_combine(A_list[[locus_index]], weights)
    factor <- tryCatch(
      chol(covariance),
      error = function(error) NULL
    )
    if (is.null(factor)) {
      return(Inf)
    }
    solved <- backsolve(
      factor,
      backsolve(factor, z_list[[locus_index]], transpose = TRUE)
    )
    value <- value + 0.5 * (
      2 * sum(log(diag(factor))) +
        sum(z_list[[locus_index]] * solved)
    )
  }
  value
}

.ld_ss_profile_identity <- function(
    source_weights,
    z_list,
    A_list,
    identity_max,
    profile_tolerance
) {
  if (identity_max == 0) {
    weights <- c(source_weights, 0)
    return(list(
      weights = weights,
      negative_log_likelihood = .ld_ss_nll(weights, z_list, A_list)
    ))
  }
  fit <- stats::optimize(
    function(identity_weight) {
      weights <- c(
        (1 - identity_weight) * source_weights,
        identity_weight
      )
      .ld_ss_nll(weights, z_list, A_list)
    },
    interval = c(0, identity_max),
    tol = profile_tolerance
  )
  list(
    weights = c((1 - fit$minimum) * source_weights, fit$minimum),
    negative_log_likelihood = fit$objective
  )
}

.ld_ss_full_fisher <- function(weights, z_list, A_list) {
  d <- length(weights)
  gradient <- numeric(d)
  fisher <- matrix(0, d, d)

  for (locus_index in seq_along(z_list)) {
    A <- A_list[[locus_index]]
    covariance <- .ld_ss_combine(A, weights)
    factor <- chol(covariance)
    precision <- chol2inv(factor)
    solved <- as.numeric(precision %*% z_list[[locus_index]])
    products <- vector("list", d)
    for (j in seq_len(d)) {
      gradient[j] <- gradient[j] + 0.5 * (
        sum(precision * A[[j]]) -
          sum(solved * (A[[j]] %*% solved))
      )
      products[[j]] <- .ld_matrix_multiply(precision, A[[j]])
      for (r in seq_len(j)) {
        value <- 0.5 * sum(products[[j]] * t(products[[r]]))
        fisher[j, r] <- fisher[j, r] + value
        if (r != j) {
          fisher[r, j] <- fisher[r, j] + value
        }
      }
    }
  }
  list(gradient = gradient, fisher = fisher)
}

.ld_ss_capped_simplex_qp <- function(fisher, linear, identity_max) {
  d <- length(linear)
  K <- d - 1L
  identity_index <- d
  tolerance <- 1e-9
  candidates <- list()

  for (mask in seq_len(2^K - 1L)) {
    active <- which(as.logical(intToBits(mask)[seq_len(K)]))
    fisher_active <- fisher[active, active, drop = FALSE]
    kkt <- rbind(
      cbind(fisher_active, rep(1, length(active))),
      c(rep(1, length(active)), 0)
    )

    fit <- .ld_symmetric_linear_solve(
      kkt,
      c(linear[active], 1),
      1e-12
    )
    if (!is.null(fit)) {
      weights <- numeric(d)
      weights[active] <- fit[seq_along(active)]
      if (all(weights >= -tolerance) &&
          abs(sum(weights) - 1) <= tolerance) {
        candidates[[length(candidates) + 1L]] <- weights
      }
    }

    linear_boundary <- linear[active] -
      fisher[active, identity_index] * identity_max
    fit <- .ld_symmetric_linear_solve(
      kkt,
      c(linear_boundary, 1 - identity_max),
      1e-12
    )
    if (!is.null(fit)) {
      weights <- numeric(d)
      weights[active] <- fit[seq_along(active)]
      weights[identity_index] <- identity_max
      if (all(weights >= -tolerance) &&
          abs(sum(weights) - 1) <= tolerance) {
        candidates[[length(candidates) + 1L]] <- weights
      }
    }

    active_with_identity <- c(active, identity_index)
    fisher_active <- fisher[
      active_with_identity,
      active_with_identity,
      drop = FALSE
    ]
    kkt_identity <- rbind(
      cbind(fisher_active, rep(1, length(active_with_identity))),
      c(rep(1, length(active_with_identity)), 0)
    )
    fit <- .ld_symmetric_linear_solve(
      kkt_identity,
      c(linear[active_with_identity], 1),
      1e-12
    )
    if (!is.null(fit)) {
      weights <- numeric(d)
      weights[active_with_identity] <- fit[seq_along(active_with_identity)]
      if (all(weights >= -tolerance) &&
          weights[identity_index] <= identity_max + tolerance &&
          abs(sum(weights) - 1) <= tolerance) {
        candidates[[length(candidates) + 1L]] <- weights
      }
    }
  }

  if (!length(candidates)) {
    stop("No feasible capped-simplex Fisher proposal was found.", call. = FALSE)
  }
  candidate_matrix <- do.call(rbind, candidates)
  candidate_matrix[abs(candidate_matrix) < tolerance] <- 0
  objectives <- apply(candidate_matrix, 1L, function(weights) {
    0.5 * sum(weights * (fisher %*% weights)) - sum(linear * weights)
  })
  as.numeric(candidate_matrix[which.min(objectives), ])
}

.ld_ss_free_information <- function(weights, A_list) {
  K <- length(weights) - 1L
  information <- matrix(0, K, K)

  for (A in A_list) {
    covariance <- .ld_ss_combine(A, weights)
    precision <- chol2inv(chol(covariance))
    differences <- lapply(
      seq_len(K),
      function(j) A[[j]] - A[[K + 1L]]
    )
    products <- lapply(
      differences,
      function(difference) .ld_matrix_multiply(precision, difference)
    )
    for (j in seq_len(K)) {
      for (r in seq_len(j)) {
        value <- 0.5 * sum(products[[j]] * t(products[[r]]))
        information[j, r] <- information[j, r] + value
        if (r != j) {
          information[r, j] <- information[r, j] + value
        }
      }
    }
  }
  dimnames(information) <- list(names(weights)[seq_len(K)], names(weights)[seq_len(K)])
  eigenvalues <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
  minimum <- min(eigenvalues)
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(eigenvalues)))
  condition <- if (minimum <= tolerance) Inf else max(eigenvalues) / minimum
  list(matrix = information, minimum = minimum, condition = condition)
}
