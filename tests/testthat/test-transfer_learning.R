test_that("cov_tl uses the analytic covariance URE weight", {
  X <- cbind(
    c(2, -2, 1, -1),
    c(1, 1, -1, -1)
  )
  S_target <- crossprod(X) / nrow(X)
  S_source <- matrix(c(2.2, 0.2, 0.2, 1.1), 2, 2)

  row_norm_squared <- rowSums(X^2)
  variance <-
    (sum(row_norm_squared^2) - nrow(X) * sum(S_target^2)) /
    (nrow(X) * (nrow(X) - 1))
  distance_squared <- sum((S_source - S_target)^2)
  expected_lambda <- min(variance / distance_squared, 1)

  fit <- cov_tl(X, S_source, center = FALSE)

  expect_s3_class(fit, "cov_tl")
  expect_equal(fit$target_covariance, S_target)
  expect_equal(fit$target_variance, variance)
  expect_equal(fit$distance_squared, distance_squared)
  expect_equal(fit$lambda, expected_lambda)
  expect_equal(
    fit$covariance,
    (1 - expected_lambda) * S_target + expected_lambda * S_source
  )
  expect_equal(fit$alpha, fit$lambda / (1 - fit$lambda))
  expect_true(fit$exact_known_mean_ure)
})

test_that("cov_tl handles the URE boundary and optional centering", {
  X <- cbind(c(2, -1, 1), c(1, 0, -1))
  S_target <- crossprod(X) / nrow(X)

  boundary <- cov_tl(X, S_target)
  expect_equal(boundary$lambda, 1)
  expect_identical(boundary$alpha, Inf)

  centered <- cov_tl(X, S_target, center = TRUE)
  X_centered <- sweep(X, 2, colMeans(X), "-")
  expect_false(centered$exact_known_mean_ure)
  expect_equal(
    centered$target_covariance,
    crossprod(X_centered) / nrow(X_centered)
  )
})

test_that("transfer-learning inputs are validated", {
  X <- matrix(rnorm(20), 10, 2)

  expect_error(cov_tl(X, diag(3)), "same number of variables")
  expect_error(
    eigen_tl(X, diag(2), rank = 0, n_source = 20),
    "rank must be between"
  )
  expect_error(
    eigen_tl(X, diag(2), rank = 1, n_source = 20, method = "median"),
    "arg"
  )
})

test_that("Frobenius distance is stable for nearly equal large matrices", {
  A <- diag(1e8, 4)
  B <- A
  B[1, 2] <- 1
  B[2, 1] <- 1

  expect_equal(
    LDRegularization:::.ld_tl_frobenius_distance_squared(A, B),
    norm(A - B, type = "F")^2
  )
  expect_equal(
    LDRegularization:::.ld_tl_frobenius_distance_squared(A, B),
    2
  )
})

test_that("Frobenius geometry also applies to nonsymmetric matrices", {
  A <- matrix(1:6, nrow = 2)
  B <- matrix(c(2, 0, 1, 4, -1, 3), nrow = 2)

  inner_product <- sum(diag(t(A) %*% B))
  expect_equal(inner_product, sum(A * B))
  expect_equal(
    LDRegularization:::.ld_tl_frobenius_distance_squared(A, B),
    sum((A - B)^2)
  )
})

test_that("source_moments_method matches direct known-mean calculations", {
  X_source <- cbind(
    c(2, -2, 1, -1, 0),
    c(1, 1, -1, -1, 0),
    c(0, 1, 2, -1, -2)
  )
  n <- nrow(X_source)
  expected_covariance <- crossprod(X_source) / n
  row_norm_squared <- rowSums(X_source^2)
  expected_fourth <- sum(row_norm_squared^2)
  expected_variance <-
    (expected_fourth - n * sum(expected_covariance^2)) /
    (n * (n - 1))

  fit <- source_moments_method(
    X_source,
    center = FALSE,
    block_size = 2,
    n_threads = 2
  )

  expect_s3_class(fit, "ld_source_moments")
  expect_equal(fit$covariance, expected_covariance, tolerance = 1e-14)
  expect_equal(fit$fourth_sum, expected_fourth, tolerance = 1e-14)
  expect_equal(fit$variance, expected_variance, tolerance = 1e-14)
  expect_true(fit$exact_known_mean_moments)
  expect_identical(fit$normalization, "1/n")
})

test_that("source moment blocks support integer data, centering, and cached S0", {
  X_source <- matrix(
    as.integer(c(0, 1, 2, 1, 0, 2, 2, 1, 0, 0, 1, 2)),
    nrow = 4
  )
  centered <- sweep(X_source, 2, colMeans(X_source), "-")
  S_source <- crossprod(centered) / nrow(centered)

  computed <- source_moments_method(
    X_source,
    center = TRUE,
    block_size = 1,
    n_threads = 1
  )
  cached <- source_moments_method(
    X_source,
    S_source = S_source,
    center = TRUE,
    block_size = 3,
    n_threads = 2
  )

  expect_equal(computed$covariance, S_source, tolerance = 1e-14)
  expect_equal(cached$covariance, S_source, tolerance = 1e-14)
  expect_equal(computed$fourth_sum, cached$fourth_sum, tolerance = 1e-14)
  expect_equal(computed$variance, cached$variance, tolerance = 1e-14)
  expect_true(computed$covariance_computed)
  expect_false(cached$covariance_computed)
  expect_false(computed$exact_known_mean_moments)
})

test_that("cov_tl selects its finite-source branch from the source input", {
  X_target <- cbind(
    c(2, -2, 1, -1),
    c(1, 1, -1, -1)
  )
  X_source <- cbind(
    c(2.2, -2.1, 0.8, -0.9, 0.1, -0.1),
    c(0.8, 1.2, -1.1, -0.9, 0.2, -0.2)
  )
  source_fit <- source_moments_method(
    X_source,
    center = FALSE,
    block_size = 2
  )
  target <- LDRegularization:::.ld_tl_target_moments(
    X_target,
    center = FALSE
  )
  distance <- sum((source_fit$covariance - target$covariance)^2)
  expected_denominator <- max(
    target$variance + source_fit$variance,
    distance
  )

  finite_source <- cov_tl(
    X_target,
    source_fit,
    center = FALSE
  )
  summary_only <- cov_tl(
    X_target,
    source_fit$covariance,
    center = FALSE
  )
  legacy <- LDRegularization:::cov_tl1(
    X_target,
    source_fit$covariance,
    center = FALSE
  )

  expect_s3_class(finite_source, "cov_tl")
  expect_identical(finite_source$weight_method, "finite_source")
  expect_equal(finite_source$denominator, expected_denominator)
  expect_equal(finite_source$lambda, target$variance / expected_denominator)
  expect_identical(summary_only$weight_method, "summary_only_ure")
  expect_equal(summary_only$lambda, legacy$lambda)
  expect_s3_class(legacy, "cov_tl1")
})

test_that("source moment inputs reject missing and incompatible data", {
  X <- matrix(1:12, 4, 3)
  X[2, 2] <- NA
  expect_error(source_moments_method(X), "finite values")
  expect_error(
    source_moments_method(matrix(1:12, 4, 3), S_source = diag(2)),
    "same number of variables"
  )
})

test_that("deprecated transfer names are not exported", {
  exports <- getNamespaceExports("LDRegularization")

  expect_false("cov_tl1" %in% exports)
  expect_false("eigspac_tl" %in% exports)
  expect_false("path_tpca_max_score" %in% exports)
  expect_false("path_tpca_one_se" %in% exports)
  expect_true("eigen_tl" %in% exports)
  expect_true("multi_source_tl" %in% exports)
})

test_that("EigenTL min path uses fold-adjusted effective source sizes", {
  X_target <- rbind(
    c(2.0, 0.0),
    c(-2.0, 0.0),
    c(1.0, 1.0),
    c(-1.0, -1.0),
    c(0.5, -0.5),
    c(-0.5, 0.5)
  )
  S_source <- diag(c(3, 0.5))
  fold_id <- c(1, 2, 3, 1, 2, 3)
  zeta <- c(0, 0.5, 1)

  fit <- eigen_tl(
    X_target,
    S_source,
    rank = 1,
    n_source = 12,
    source_fraction_grid = zeta,
    center = FALSE,
    fold_id = fold_id,
    method = "min"
  )

  expect_s3_class(fit, "eigen_tl")
  expect_identical(fit$method, "min")
  expect_equal(fit$selected_index, which.max(fit$mean_scores))
  expect_equal(fit$selected_zeta, zeta[fit$selected_index])
  expect_equal(fit$fold_weights[, 1], rep(0, 3))
  expect_equal(fit$fold_weights[, 3], rep(12 / (4 + 12), 3))
  expect_equal(fit$full_data_weights, 12 * zeta / (6 + 12 * zeta))
  expect_equal(
    fit$selected_weight,
    fit$full_data_weights[fit$selected_index]
  )
  expect_equal(sum(diag(fit$projector)), 1, tolerance = 1e-10)
  expect_equal(fit$projector %*% fit$projector, fit$projector,
               tolerance = 1e-10)
  expect_false(fit$source_fourth_moment_used)
  expect_equal(
    fit$reconstruction_gain,
    fit$selected_score - fit$target_score
  )
  expect_error(
    eigen_tl(
      X_target, S_source, rank = 1, fold_id = fold_id, method = "min"
    ),
    "n_source is required"
  )
})

test_that("Teacher-B paired one-SE rule matches its individual score formula", {
  X_target <- rbind(
    c(2.0, 0.1),
    c(-1.8, -0.2),
    c(1.1, 0.8),
    c(-0.9, -1.0),
    c(0.4, -0.7),
    c(-0.5, 0.6)
  )
  X_source <- rbind(
    c(2.3, 0.4),
    c(-2.1, -0.3),
    c(1.5, 0.2),
    c(-1.4, -0.1),
    c(0.8, 0.5),
    c(-0.7, -0.4),
    c(0.2, 0.7),
    c(-0.3, -0.6)
  )
  source_fit <- source_moments_method(
    X_source,
    center = FALSE,
    block_size = 3
  )
  fold_id <- c(1, 2, 3, 1, 2, 3)
  zeta <- c(0, 0.25, 0.5, 1)

  min_fit <- eigen_tl(
    X_target,
    source_fit,
    rank = 1,
    source_fraction_grid = zeta,
    center = FALSE,
    fold_id = fold_id,
    method = "min"
  )
  one_se <- eigen_tl(
    X_target,
    source_fit,
    rank = 1,
    source_fraction_grid = zeta,
    center = FALSE,
    fold_id = fold_id
  )

  paired <- one_se$observation_scores[, one_se$best_index] -
    one_se$observation_scores
  expected_se <- apply(paired, 2, stats::sd) / sqrt(nrow(X_target))
  expected_se[one_se$best_index] <- 0
  tolerance <- .Machine$double.eps * max(1, abs(one_se$mean_scores))
  expected_competitive <- which(
    one_se$score_gaps <= expected_se + tolerance
  )

  expect_s3_class(one_se, "eigen_tl")
  expect_identical(one_se$method, "one_se")
  expect_equal(one_se$observation_scores, min_fit$observation_scores)
  expect_equal(one_se$paired_standard_errors, expected_se)
  expect_equal(one_se$competitive_set, expected_competitive)
  expect_equal(one_se$selected_index, max(expected_competitive))
  expect_gte(one_se$selected_index, one_se$best_index)
  expect_equal(one_se$n_source, nrow(X_source))
  expect_false(one_se$source_fourth_moment_used)
  expect_error(
    eigen_tl(
      X_target,
      source_fit,
      rank = 1,
      n_source = nrow(X_source) + 1,
      fold_id = fold_id
    ),
    "does not match"
  )
})

test_that("EigenTL inputs enforce the proposed grid and fold rules", {
  X <- matrix(seq_len(18), 6, 3)
  S <- diag(3)

  expect_error(
    eigen_tl(
      X, S, rank = 1, n_source = 10,
      source_fraction_grid = c(0, 0.5), fold_id = rep(1:2, 3), method = "min"
    ),
    "include endpoints"
  )
  expect_error(
    eigen_tl(
      X, S, rank = 1, n_source = 10,
      source_fraction_grid = c(0, 0.5, 0.5, 1), fold_id = rep(1:2, 3),
      method = "min"
    ),
    "strictly increasing"
  )
  expect_error(
    eigen_tl(
      X, S, rank = 1, n_source = 10,
      fold_id = rep(1:2, 3), standard_error_multiplier = -1
    ),
    "nonnegative"
  )
  expect_error(
    eigen_tl(
      X, S, rank = 1, n_source = 10, fold_id = rep(1, 6), method = "min"
    ),
    "at least two"
  )
})

test_that("multi_source_tl uses CovTL shrinkage gains", {
  X <- rbind(
    c(2.0, 0.0),
    c(-2.0, 0.0),
    c(1.0, 1.0),
    c(-1.0, -1.0),
    c(0.5, -0.5),
    c(-0.5, 0.5)
  )
  fits <- list(
    EUR = cov_tl(X, diag(c(3.0, 0.5)), center = FALSE),
    AFR = cov_tl(X, diag(c(2.0, 1.0)), center = FALSE),
    AMR = cov_tl(X, diag(c(1.5, 1.5)), center = FALSE)
  )

  fit <- multi_source_tl(fits)
  expected_weights <- vapply(fits, function(x) x$lambda, numeric(1))
  expected_weights <- expected_weights / sum(expected_weights)
  expected_covariance <- matrix(0, 2, 2)
  for (source_index in seq_along(fits)) {
    expected_covariance <- expected_covariance +
      expected_weights[source_index] * fits[[source_index]]$covariance
  }

  expect_s3_class(fit, "multi_source_tl")
  expect_identical(fit$fit_family, "cov_tl")
  expect_identical(fit$aggregation, "shrinkage_gain")
  expect_equal(fit$weights, expected_weights)
  expect_equal(fit$covariance, expected_covariance)
  expect_equal(sum(fit$weights), 1)
  expect_equal(
    fit$effective_source_weights,
    fit$weights * vapply(fits, function(x) x$lambda, numeric(1))
  )
  expect_equal(
    fit$target_weight + sum(fit$effective_source_weights),
    1
  )
})

test_that("multi_source_tl uses EigenTL reconstruction gains", {
  X <- rbind(
    c(2.0, 0.1),
    c(-1.8, -0.2),
    c(1.1, 0.8),
    c(-0.9, -1.0),
    c(0.4, -0.7),
    c(-0.5, 0.6)
  )
  fold_id <- c(1, 2, 3, 1, 2, 3)
  source <- list(
    EUR = diag(c(3.0, 0.5)),
    AFR = matrix(c(2.0, 0.4, 0.4, 1.0), 2, 2),
    AMR = diag(c(1.2, 1.1))
  )
  fits <- lapply(source, function(S) {
    eigen_tl(
      X,
      S,
      rank = 1,
      n_source = 20,
      source_fraction_grid = c(0, 0.5, 1),
      center = FALSE,
      fold_id = fold_id
    )
  })
  gains <- c(EUR = 0.1, AFR = 0.3, AMR = -0.1)
  for (source_name in names(fits)) {
    fits[[source_name]]$reconstruction_gain <- gains[source_name]
  }

  fit <- multi_source_tl(fits)

  expect_s3_class(fit, "multi_source_tl")
  expect_identical(fit$fit_family, "eigen_tl")
  expect_identical(fit$aggregation, "reuse_cv_reconstruction_gain")
  expect_equal(fit$signed_gains, gains)
  expect_equal(fit$weights, c(EUR = 0.25, AFR = 0.75, AMR = 0))
  expect_equal(sum(diag(fit$projector)), 1, tolerance = 1e-10)
  expect_equal(
    fit$projector %*% fit$projector,
    fit$projector,
    tolerance = 1e-10
  )
  expect_equal(sum(fit$weights), 1)
})

test_that("multi_source_tl validates fit families, folds, and priors", {
  X <- matrix(rnorm(24), 12, 2)
  cov_fit <- cov_tl(X, diag(2))
  fold_id <- rep(1:3, 4)
  eigen_fit <- eigen_tl(
    X,
    diag(2),
    rank = 1,
    n_source = 20,
    fold_id = fold_id
  )

  expect_error(
    multi_source_tl(list(EUR = cov_fit, AFR = eigen_fit)),
    "only cov_tl objects or only eigen_tl objects"
  )
  expect_error(
    multi_source_tl(list(EUR = cov_fit, AFR = cov_fit), prior = c(1, -1)),
    "nonnegative"
  )

  other_fold <- eigen_fit
  other_fold$fold_id <- rev(other_fold$fold_id)
  expect_error(
    multi_source_tl(list(EUR = eigen_fit, AFR = other_fold)),
    "identical target folds"
  )

  other_method <- eigen_fit
  other_method$method <- "min"
  expect_error(
    multi_source_tl(list(EUR = eigen_fit, AFR = other_method)),
    "same method"
  )
})
