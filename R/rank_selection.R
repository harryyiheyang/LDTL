#' Select a target-only transfer-PCA rank by projected bootstrap
#'
#' Select a stable compression dimension from `X_target` without assuming a
#' finite-spike covariance model. The requested cumulative fractions are
#' operational variance-retention targets, not estimates of a true algebraic
#' rank. This distinction permits the selected rank to be a substantial
#' fraction of the ambient dimension.
#'
#' The native algorithm first computes one adaptive randomized target
#' eigenspace. Directions below `anchor_fraction - guard_fraction` are treated
#' as a fixed, well-estimated core. The guard band around the anchor and all
#' subsequent directions remain flexible. Each bootstrap sample fits the
#' flexible directions on its in-bag individuals and evaluates cumulative
#' captured variance on its out-of-bag individuals. All calculations use
#' precomputed target scores in the flexible band; neither a bootstrap `p` by
#' `p` covariance nor a repeated full spectral decomposition is formed.
#'
#' For a requested fraction `q`, the selected rank is the smallest candidate
#' whose lower out-of-bag cumulative-variance quantile is at least `q`.
#' Flexible-subspace stability is then reported as an audit status at that
#' rank. It deliberately does not increase the selected rank: normalized
#' projector loss can fall mechanically when a candidate approaches the full
#' flexible space, so searching upward for stability would favor excessively
#' large dimensions. The calculation is conditional on the fixed core and
#' randomized reference band.
#'
#' @param X_target Target individual-level matrix, with observations in rows.
#' @param target_fraction_grid Increasing cumulative target-variance goals to
#'   report. The default returns a path from 50 to 90 percent.
#' @param selection_fraction Fraction whose selected rank is returned in
#'   `rank`. It is automatically added to the reported path.
#' @param anchor_fraction Computational anchor fraction. This is not the final
#'   variance target.
#' @param guard_fraction Fraction below the anchor left flexible to avoid
#'   treating the anchor boundary as known. The fixed core ends at
#'   `anchor_fraction - guard_fraction`.
#' @param n_boot Number of individual bootstrap samples.
#' @param confidence One-sided confidence level for the cumulative-variance
#'   lower quantile.
#' @param stability_quantile Upper bootstrap quantile used for subspace loss.
#' @param stability_tolerance Maximum accepted normalized flexible-projector
#'   loss at `stability_quantile`.
#' @param basis_margin Extra reference variance beyond the largest requested
#'   fraction. The adaptive randomized basis expands until it captures this
#'   buffer or reaches the sample algebraic-rank limit.
#' @param initial_rank Initial randomized reference rank.
#' @param rank_block Minimum reference-rank expansion when more variance is
#'   needed.
#' @param oversample Randomized range-finder oversampling dimension.
#' @param bootstrap_extra Extra flexible directions beyond the reference rank
#'   needed for the largest requested fraction.
#' @param power_iterations Number of randomized subspace power iterations.
#' @param center Whether to subtract the target sample mean and recenter every
#'   bootstrap sample.
#' @param n_threads Number of OpenMP bootstrap threads. Zero uses the runtime
#'   default.
#' @param seed Nonnegative integer-valued seed. Each bootstrap replication has
#'   its own deterministic stream, so results do not depend on thread count.
#'
#' @return An object of class `tl_rank_bootstrap`. `rank` is the selection for
#'   `selection_fraction`; `rank_path` reports all requested fractions and
#'   selection statuses; `diagnostics` contains the full candidate curves.
#' @export
select_tl_rank_bootstrap <- function(
    X_target,
    target_fraction_grid = c(.5, .6, .7, .8, .9),
    selection_fraction = .8,
    anchor_fraction = .5,
    guard_fraction = .05,
    n_boot = 100L,
    confidence = .9,
    stability_quantile = .9,
    stability_tolerance = .1,
    basis_margin = .05,
    initial_rank = 64L,
    rank_block = 64L,
    oversample = 20L,
    bootstrap_extra = 64L,
    power_iterations = 2L,
    center = TRUE,
    n_threads = 0L,
    seed = 20260821
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
  .ld_rank_fraction(anchor_fraction, "anchor_fraction", open = TRUE)
  .ld_rank_fraction(selection_fraction, "selection_fraction", open = TRUE)
  .ld_rank_fraction(confidence, "confidence", open = TRUE)
  .ld_rank_fraction(stability_quantile, "stability_quantile", open = TRUE)
  if (length(guard_fraction) != 1L || !is.finite(guard_fraction) ||
      guard_fraction < 0 || guard_fraction >= anchor_fraction) {
    stop(
      "guard_fraction must be nonnegative and smaller than anchor_fraction.",
      call. = FALSE
    )
  }
  if (length(basis_margin) != 1L || !is.finite(basis_margin) ||
      basis_margin < 0 || basis_margin >= 1) {
    stop("basis_margin must be in [0, 1).", call. = FALSE)
  }
  if (length(stability_tolerance) != 1L ||
      !is.finite(stability_tolerance) || stability_tolerance < 0 ||
      stability_tolerance > 1) {
    stop("stability_tolerance must be in [0, 1].", call. = FALSE)
  }

  fractions <- sort(unique(c(
    as.numeric(target_fraction_grid),
    as.numeric(selection_fraction)
  )))
  if (!length(fractions) || any(!is.finite(fractions)) ||
      any(fractions <= 0) || any(fractions >= 1) ||
      any(fractions < anchor_fraction)) {
    stop(
      "target fractions must lie in [anchor_fraction, 1).",
      call. = FALSE
    )
  }
  integer_parameters <- list(
    n_boot = n_boot,
    initial_rank = initial_rank,
    rank_block = rank_block,
    oversample = oversample,
    bootstrap_extra = bootstrap_extra,
    power_iterations = power_iterations,
    n_threads = n_threads
  )
  minimums <- c(2, 1, 1, 0, 0, 0, 0)
  for (index in seq_along(integer_parameters)) {
    .ld_rank_integer(
      integer_parameters[[index]],
      names(integer_parameters)[index],
      minimums[index]
    )
  }
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed != round(seed) || seed > 2^53) {
    stop("seed must be an integer-valued scalar between 0 and 2^53.", call. = FALSE)
  }

  largest_fraction <- max(fractions)
  basis_fraction <- min(.995, largest_fraction + basis_margin)
  if (basis_fraction <= largest_fraction) {
    stop(
      "basis_margin leaves no reference buffer below one; reduce the largest target fraction.",
      call. = FALSE
    )
  }
  native <- cpp_projected_bootstrap_rank(
    X_target = X_target,
    center = isTRUE(center),
    n_boot = as.integer(n_boot),
    anchor_fraction = anchor_fraction,
    guard_fraction = guard_fraction,
    target_fraction = largest_fraction,
    basis_fraction = basis_fraction,
    initial_rank = as.integer(initial_rank),
    rank_block = as.integer(rank_block),
    oversample = as.integer(oversample),
    bootstrap_extra = as.integer(bootstrap_extra),
    power_iterations = as.integer(power_iterations),
    n_threads = as.integer(n_threads),
    seed = as.double(seed)
  )

  ranks <- native$core_rank + seq_len(native$flexible_keep)
  reference_cumulative <- cumsum(native$eigenvalues) /
    native$total_variance
  reference_candidate <- reference_cumulative[ranks]
  bootstrap_fraction <- native$bootstrap_cumulative_fraction
  bootstrap_loss <- native$bootstrap_flexible_subspace_loss
  lower_fraction <- apply(
    bootstrap_fraction,
    2L,
    stats::quantile,
    probs = 1 - confidence,
    names = FALSE,
    type = 8
  )
  median_fraction <- apply(
    bootstrap_fraction,
    2L,
    stats::median
  )
  mean_loss <- colMeans(bootstrap_loss)
  upper_loss <- apply(
    bootstrap_loss,
    2L,
    stats::quantile,
    probs = stability_quantile,
    names = FALSE,
    type = 8
  )
  diagnostics <- data.frame(
    rank = ranks,
    reference_cumulative_fraction = reference_candidate,
    bootstrap_lower_cumulative_fraction = lower_fraction,
    bootstrap_median_cumulative_fraction = median_fraction,
    bootstrap_mean_flexible_loss = mean_loss,
    bootstrap_upper_flexible_loss = upper_loss,
    stable = upper_loss <= stability_tolerance
  )

  rank_path <- do.call(rbind, lapply(fractions, function(fraction) {
    reaches <- which(lower_fraction >= fraction)
    if (!length(reaches)) {
      selected <- length(ranks)
      status <- "target_not_reached_in_reference_band"
    } else {
      selected <- reaches[1L]
      status <- if (upper_loss[selected] <= stability_tolerance) {
        "stable"
      } else {
        "target_reached_but_unstable"
      }
    }
    data.frame(
      target_fraction = fraction,
      rank = ranks[selected],
      lower_cumulative_fraction = lower_fraction[selected],
      upper_flexible_loss = upper_loss[selected],
      status = status,
      stringsAsFactors = FALSE
    )
  }))
  selection_index <- match(selection_fraction, rank_path$target_fraction)

  structure(
    list(
      rank = rank_path$rank[selection_index],
      selection_fraction = selection_fraction,
      selection_status = rank_path$status[selection_index],
      rank_path = rank_path,
      diagnostics = diagnostics,
      eigenvalues = native$eigenvalues,
      vectors = native$vectors,
      total_variance = native$total_variance,
      core_rank = native$core_rank,
      anchor_rank = native$anchor_rank,
      reference_rank = native$reference_rank,
      flexible_keep = native$flexible_keep,
      basis_fraction_requested = basis_fraction,
      basis_fraction_captured = native$basis_fraction_captured,
      bootstrap_cumulative_fraction = bootstrap_fraction,
      bootstrap_flexible_subspace_loss = bootstrap_loss,
      n = native$n,
      p = native$p,
      n_boot = native$n_boot,
      center = isTRUE(center),
      confidence = confidence,
      stability_quantile = stability_quantile,
      stability_tolerance = stability_tolerance,
      anchor_fraction = anchor_fraction,
      guard_fraction = guard_fraction,
      threads_used = native$threads_used,
      openmp = native$openmp,
      source_data_used = FALSE,
      conditional_on_fixed_core = TRUE,
      selection_definition =
        "smallest rank whose bootstrap out-of-bag lower cumulative-variance quantile reaches the user target; stability is audited but does not inflate rank"
    ),
    class = "tl_rank_bootstrap"
  )
}

#' @export
print.tl_rank_bootstrap <- function(x, ...) {
  cat("Target-only projected-bootstrap rank selection\n")
  cat("Selected rank:", x$rank, "for fraction", x$selection_fraction, "\n")
  cat("Status:", x$selection_status, "\n")
  cat(
    "Fixed core / anchor / reference ranks:",
    x$core_rank, "/", x$anchor_rank, "/", x$reference_rank, "\n"
  )
  print(x$rank_path, row.names = FALSE)
  invisible(x)
}

.ld_rank_fraction <- function(value, name, open) {
  valid <- length(value) == 1L && is.finite(value)
  if (isTRUE(open)) {
    valid <- valid && value > 0 && value < 1
  } else {
    valid <- valid && value >= 0 && value <= 1
  }
  if (!valid) {
    stop(name, " must be a finite fraction strictly between zero and one.", call. = FALSE)
  }
  invisible(value)
}

.ld_rank_integer <- function(value, name, minimum) {
  if (length(value) != 1L || !is.finite(value) || value < minimum ||
      value != round(value)) {
    stop(name, " must be an integer at least ", minimum, ".", call. = FALSE)
  }
  invisible(value)
}
