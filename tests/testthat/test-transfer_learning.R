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
  expect_equal(fit$target_variance, variance)
  expect_equal(fit$distance_squared, distance_squared)
  expect_equal(fit$lambda, expected_lambda)
  expect_equal(
    fit$covariance,
    (1 - expected_lambda) * S_target + expected_lambda * S_source
  )
  expect_equal(fit$alpha, fit$lambda / (1 - fit$lambda))
  expect_null(fit$target_covariance)
  expect_null(fit$source_covariance)
  expect_null(fit$n_target)
  expect_null(fit$n_source)
  expect_null(fit$center)
})

test_that("cov_tl handles the URE boundary and optional centering", {
  X <- cbind(c(2, -1, 1), c(1, 0, -1))
  S_target <- crossprod(X) / nrow(X)

  boundary <- cov_tl(X, S_target)
  expect_equal(boundary$lambda, 1)
  expect_identical(boundary$alpha, Inf)

  centered <- cov_tl(X, S_target, center = TRUE)
  expect_true(all(is.finite(centered$covariance)))
})

test_that("transfer-learning inputs are validated", {
  X <- matrix(rnorm(20), 10, 2)

  expect_error(cov_tl(X, diag(3)), "same number of variables")
  expect_error(
    eigen_tl(X, diag(2), rank = 0, n_source = 20),
    "positive integer"
  )
  expect_error(
    eigen_tl(X, diag(2), rank = 1.5, n_source = 20),
    "fixed rank must be an integer"
  )
  expect_error(
    eigen_tl(X, diag(2), rank = 1, n_source = 20, method = "median"),
    "arg"
  )
})

test_that("EigenTL rank accepts a fixed dimension or target variance fraction", {
  X <- rbind(
    c(4, 0.2, 0.1), c(-4, -0.2, -0.1),
    c(3, 0.1, 0), c(-3, -0.1, 0),
    c(2, 0, 0.1), c(-2, 0, -0.1)
  )
  fold_id <- rep(1:3, each = 2)
  source <- diag(c(10, 1, 0.5))
  sources <- list(
    EUR = source,
    AFR = diag(c(8, 1.2, 0.7))
  )

  fixed <- eigen_tl(
    X, source, rank = 2, n_source = 20,
    fold_id = fold_id, center = FALSE, method = "min"
  )
  fractional <- eigen_tl(
    X, source, rank = 0.8, n_source = 20,
    fold_id = fold_id, center = FALSE, method = "min"
  )
  multi_fractional <- multisource_eigen_tl(
    X, sources, rank = 0.8,
    fold_id = fold_id, center = FALSE, method = "min"
  )
  E <- LDTL:::.ld_eigen(crossprod(X) / nrow(X))
  expected_rank <- which(
    cumsum(pmax(E$values, 0)) / sum(pmax(E$values, 0)) >= 0.8
  )[1L]

  expect_equal(fixed$effective_rank, 2L)
  expect_equal(fixed$fold_effective_ranks, rep(2L, 3))
  expect_equal(ncol(fixed$vectors), 2L)
  expect_equal(fractional$effective_rank, expected_rank)
  expect_equal(
    fractional$fold_effective_ranks,
    rep(expected_rank, 3)
  )
  expect_equal(ncol(fractional$vectors), expected_rank)
  expect_equal(multi_fractional$effective_rank, expected_rank)
  expect_equal(
    multi_fractional$fold_effective_ranks,
    rep(expected_rank, 3)
  )
  expect_equal(ncol(multi_fractional$vectors), expected_rank)
})

test_that("RSpectra EigenTL agrees with the full eigensolver", {
  skip_if_not_installed("RSpectra")
  set.seed(29)
  X <- matrix(rnorm(640), 80, 8)
  X <- sweep(X, 2, seq(0.5, 1.2, length.out = 8), "*")
  Z1 <- X + matrix(rnorm(640, sd = 0.15), 80, 8)
  Z2 <- X + matrix(rnorm(640, sd = 0.3), 80, 8)
  source1 <- source_moments(Z1, center = TRUE)
  sources <- list(
    EUR = source1,
    AFR = source_moments(Z2, center = TRUE)
  )
  fold_id <- rep(1:5, length.out = nrow(X))
  grid <- c(0, 0.25, 0.5, 0.75, 1)

  single_full <- eigen_tl(
    X, source1, rank = 3, fold_id = fold_id,
    source_fraction_grid = grid, method = "one_se",
    eigen_solver = "full"
  )
  set.seed(29)
  single_partial <- eigen_tl(
    X, source1, rank = 3, fold_id = fold_id,
    source_fraction_grid = grid, method = "one_se",
    eigen_solver = "rspectra"
  )
  multi_full <- multisource_eigen_tl(
    X, sources, rank = 3, fold_id = fold_id,
    transfer_weight_grid = grid, method = "one_se",
    eigen_solver = "full"
  )
  set.seed(29)
  multi_partial <- multisource_eigen_tl(
    X, sources, rank = 3, fold_id = fold_id,
    transfer_weight_grid = grid, method = "one_se",
    eigen_solver = "rspectra"
  )

  expect_equal(
    single_partial$observation_scores,
    single_full$observation_scores,
    tolerance = 1e-7
  )
  expect_equal(single_partial$selected_index, single_full$selected_index)
  expect_equal(
    single_partial$eigenvalues,
    single_full$eigenvalues,
    tolerance = 1e-8
  )
  expect_equal(
    single_partial$projector,
    single_full$projector,
    tolerance = 1e-7
  )
  expect_equal(
    multi_partial$observation_scores,
    multi_full$observation_scores,
    tolerance = 1e-7
  )
  expect_equal(multi_partial$selected_index, multi_full$selected_index)
  expect_equal(
    multi_partial$eigenvalues,
    multi_full$eigenvalues,
    tolerance = 1e-8
  )
  expect_equal(
    multi_partial$projector,
    multi_full$projector,
    tolerance = 1e-7
  )
})

test_that("Frobenius distance is stable for nearly equal large matrices", {
  A <- diag(1e8, 4)
  B <- A
  B[1, 2] <- 1
  B[2, 1] <- 1

  expect_equal(
    LDTL:::.ld_tl_frobenius_distance_squared(A, B),
    norm(A - B, type = "F")^2
  )
  expect_equal(
    LDTL:::.ld_tl_frobenius_distance_squared(A, B),
    2
  )
})

test_that("Frobenius geometry also applies to nonsymmetric matrices", {
  A <- matrix(1:6, nrow = 2)
  B <- matrix(c(2, 0, 1, 4, -1, 3), nrow = 2)

  inner_product <- sum(diag(t(A) %*% B))
  expect_equal(inner_product, sum(A * B))
  expect_equal(
    LDTL:::.ld_tl_frobenius_distance_squared(A, B),
    sum((A - B)^2)
  )
})

test_that("source_moments matches direct known-mean calculations", {
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

  fit <- source_moments(
    X_source,
    center = FALSE,
    block_size = 2,
    n_threads = 2
  )

  expect_s3_class(fit, "ldtl_source_moments")
  expect_equal(fit$R_source, expected_covariance, tolerance = 1e-14)
  expect_equal(fit$fourth_sum, expected_fourth, tolerance = 1e-14)
  expect_equal(fit$variance, expected_variance, tolerance = 1e-14)
  expect_true(fit$exact_known_mean_moments)
  expect_equal(fit$n_source, n)
})

test_that("source moment blocks support integer data, centering, and cached S0", {
  X_source <- matrix(
    as.integer(c(0, 1, 2, 1, 0, 2, 2, 1, 0, 0, 1, 2)),
    nrow = 4
  )
  centered <- sweep(X_source, 2, colMeans(X_source), "-")
  S_source <- crossprod(centered) / nrow(centered)

  computed <- source_moments(
    X_source,
    center = TRUE,
    block_size = 1,
    n_threads = 1
  )
  cached <- source_moments(
    X_source,
    R_source = S_source,
    center = TRUE,
    block_size = 3,
    n_threads = 2
  )

  expect_equal(computed$R_source, S_source, tolerance = 1e-14)
  expect_equal(cached$R_source, S_source, tolerance = 1e-14)
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
  source_fit <- source_moments(
    X_source,
    center = FALSE,
    block_size = 2
  )
  target <- LDTL:::.ld_tl_target_moments(
    X_target,
    center = FALSE
  )
  distance <- sum((source_fit$R_source - target$covariance)^2)
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
    source_fit$R_source,
    center = FALSE
  )
  expect_s3_class(finite_source, "cov_tl")
  expect_identical(finite_source$weight_method, "finite_source")
  expect_equal(finite_source$denominator, expected_denominator)
  expect_equal(finite_source$lambda, target$variance / expected_denominator)
  expect_identical(summary_only$weight_method, "summary_only_ure")
  expect_equal(
    summary_only$lambda,
    min(target$variance / distance, 1)
  )
})

test_that("source moment inputs reject missing and incompatible data", {
  X <- matrix(1:12, 4, 3)
  X[2, 2] <- NA
  expect_error(source_moments(X), "finite values")
  expect_error(
    source_moments(matrix(1:12, 4, 3), R_source = diag(2)),
    "same number of variables"
  )
})

test_that("deprecated transfer names are not exported", {
  exports <- getNamespaceExports("LDTL")

  expect_false("cov_tl1" %in% exports)
  expect_false("eigspac_tl" %in% exports)
  expect_false("path_tpca_max_score" %in% exports)
  expect_false("path_tpca_one_se" %in% exports)
  expect_true("eigen_tl" %in% exports)
  expect_true("cov_tl" %in% exports)
  expect_true("multisource_cov_tl" %in% exports)
  expect_true("multisource_eigen_tl" %in% exports)
  expect_true("source_moments" %in% exports)
  expect_true("source_moments_plink" %in% exports)
  expect_false("multi_source_tl" %in% exports)
  expect_false("multi_source_eigen_tl" %in% exports)
  expect_false("source_moments_method" %in% exports)
  expect_false("source_moments_raw" %in% exports)
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
  expect_null(fit$target_covariance)
  expect_null(fit$source_covariance)
  expect_null(fit$n_target)
  expect_null(fit$n_source)
  expect_null(fit$center)
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
  source_fit <- source_moments(
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
  expect_equal(one_se$observation_scores, min_fit$observation_scores)
  expect_equal(one_se$paired_standard_errors, expected_se)
  expect_equal(one_se$competitive_set, expected_competitive)
  expect_equal(one_se$selected_index, max(expected_competitive))
  expect_gte(one_se$selected_index, one_se$best_index)
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

test_that("multisource_cov_tl solves the joint CovTL regression", {
  X <- matrix(0, 6, 3)
  X[cbind(seq_len(6), rep(seq_len(3), each = 2))] <-
    rep(c(sqrt(6), -sqrt(6)), 3)
  target <- 2 * diag(3)
  directions <- list(
    EUR = diag(c(5, 0, 0)),
    AFR = diag(c(0, 5, 0)),
    AMR = diag(c(0, 0, 5))
  )
  sources <- lapply(directions, function(direction) target + direction)

  fit <- multisource_cov_tl(X, sources, center = FALSE)
  expected_weights <- rep(fit$target_variance / 25, 3)
  names(expected_weights) <- names(sources)
  expected_covariance <-
    (1 - sum(expected_weights)) * target
  for (source_name in names(sources)) {
    expected_covariance <- expected_covariance +
      expected_weights[source_name] * sources[[source_name]]
  }

  expect_s3_class(fit, "multisource_cov_tl")
  expect_equal(unname(fit$gram), 25 * diag(3))
  expect_equal(fit$source_weights, expected_weights)
  expect_equal(fit$target_weight, 1 - sum(expected_weights))
  expect_equal(fit$covariance, expected_covariance)
  expect_equal(
    fit$target_weight + sum(fit$source_weights),
    1
  )
  expect_null(fit$target_covariance)
  expect_null(fit$source_covariance)
  expect_null(fit$source_covariances)
  expect_null(fit$sources)
  expect_null(fit$n_target)
  expect_null(fit$n_source)
  expect_null(fit$center)
})

test_that("multisource_cov_tl handles the simplex boundary and singular Gram", {
  X <- rbind(c(2, 0), c(-2, 0), c(0, 2), c(0, -2))
  target <- 2 * diag(2)
  direction <- diag(c(0.5, 0))
  sources <- lapply(
    c(EUR = 1, AFR = 1, AMR = 1),
    function(unused) target + direction
  )

  fit <- multisource_cov_tl(X, sources, center = FALSE)

  expect_equal(sum(fit$source_weights), 1)
  expect_true(all(fit$source_weights >= 0))
  expect_equal(fit$target_weight, 0)
  expect_equal(fit$covariance, target + direction)
})

test_that("multisource_cov_tl accepts legacy source moments without retaining them", {
  X <- rbind(
    c(2, 0), c(-2, 0), c(0, 2), c(0, -2), c(1, 1), c(-1, -1)
  )
  target <- crossprod(X) / nrow(X)
  values <- list(
    EUR = c(variance = 0.1, fourth_sum = 101, n = 100),
    AFR = c(variance = 0.2, fourth_sum = 202, n = 200)
  )
  sources <- lapply(values, function(value) {
    structure(
      list(
        covariance = target,
        variance = unname(value["variance"]),
        fourth_sum = unname(value["fourth_sum"]),
        n = unname(value["n"]),
        p = 2,
        exact_known_mean_moments = TRUE
      ),
      class = "ld_source_moments"
    )
  })

  fit <- multisource_cov_tl(X, sources, center = FALSE)

  expect_equal(diag(fit$gram), c(EUR = 0, AFR = 0))
  expect_null(fit$source_variance)
  expect_null(fit$source_fourth_sum)
  expect_null(fit$n_source)
  expect_null(fit$sources)
})

test_that("multisource_cov_tl rejects incompatible source summaries and fits", {
  X <- matrix(rnorm(24), 12, 2)
  cov_fit <- cov_tl(X, diag(2))

  expect_error(
    multisource_cov_tl(X, list(EUR = diag(2), AFR = diag(3))),
    "same number of variables"
  )
  expect_error(
    multisource_cov_tl(X, list(EUR = cov_fit, AFR = diag(2))),
    "not fitted transfer-learning objects"
  )
})

test_that("multi-source composition recovers an exact convex source mixture", {
  S1 <- diag(c(4, 1, 2))
  S2 <- diag(c(1, 4, 2))
  S3 <- diag(c(7, 7, 1))
  T <- 0.4 * S1 + 0.6 * S2
  sources <- list(S1 = S1, S2 = S2, S3 = S3)
  A <- LDTL:::cpp_multi_source_gram(T, sources)

  fit <- LDTL:::.ld_multi_source_composition_qp(A)
  R <- LDTL:::cpp_multi_source_combine(
    T,
    sources,
    fit$coefficients
  )

  expect_equal(fit$coefficients, c(0.4, 0.6, 0), tolerance = 1e-10)
  expect_equal(sum(fit$coefficients), 1)
  expect_equal(R, T, tolerance = 1e-10)
  expect_equal(fit$objective, 0, tolerance = 1e-10)
})

test_that("multi-source EigenTL relearns target moments and direction by fold", {
  X <- rbind(
    c(2.0, 0.1, 0.0),
    c(-1.8, -0.2, 0.1),
    c(1.2, 0.9, -0.1),
    c(-1.0, -1.1, 0.2),
    c(0.3, -0.8, 1.2),
    c(-0.4, 0.7, -1.0),
    c(1.5, 0.2, 0.3),
    c(-1.3, -0.1, -0.4),
    c(0.2, 1.4, 0.5)
  )
  sources <- list(
    EUR = diag(c(3.0, 0.8, 0.5)),
    AFR = diag(c(0.8, 3.0, 0.5)),
    AMR = matrix(c(2.0, 0.4, 0.0, 0.4, 2.0, 0.0, 0.0, 0.0, 0.7), 3)
  )
  fold_id <- rep(1:3, 3)
  weights <- c(0, 0.5, 1)

  fit <- multisource_eigen_tl(
    X,
    sources,
    rank = 1,
    transfer_weight_grid = weights,
    center = FALSE,
    fold_id = fold_id,
    method = "min"
  )
  cov_fit <- multisource_cov_tl(X, sources, center = FALSE)

  expected_variances <- numeric(3)
  expected_fold_weights <- matrix(NA_real_, 3, 3)
  expected_covtl_weights <- matrix(NA_real_, 3, 3)
  for (v in 1:3) {
    X_fit <- X[fold_id != v, , drop = FALSE]
    target <- LDTL:::.ld_tl_target_moments(
      X_fit,
      center = FALSE
    )
    A <- LDTL:::cpp_multi_source_gram(
      target$covariance,
      sources
    )
    direction <- LDTL:::.ld_multi_source_composition_qp(A)
    covtl <- LDTL:::.ld_multi_source_qp(
      A,
      rep(target$variance, 3)
    )
    expected_variances[v] <- target$variance
    expected_fold_weights[v, ] <- direction$coefficients
    expected_covtl_weights[v, ] <- covtl$coefficients
  }

  expect_s3_class(fit, "multisource_eigen_tl")
  expect_s3_class(fit, "eigen_tl")
  expect_equal(fit$path_weights, weights)
  expect_equal(fit$fold_target_variances, expected_variances)
  expect_equal(unname(fit$fold_source_weights), expected_fold_weights)
  expect_equal(
    unname(fit$fold_covtl_source_weights),
    expected_covtl_weights
  )
  expect_equal(
    fit$fold_covtl_transfer_weights,
    rowSums(expected_covtl_weights)
  )
  expect_equal(unname(rowSums(fit$fold_source_weights)), rep(1, 3))
  expect_equal(sum(fit$source_weights), 1)
  expect_true(all(fit$source_weights >= 0))
  expect_equal(fit$selected_index, which.max(fit$mean_scores))
  expect_equal(fit$selected_weight, weights[fit$selected_index])
  target_full <- crossprod(X) / nrow(X)
  source_full <- matrix(0, ncol(X), ncol(X))
  for (s in seq_along(sources)) {
    source_full <- source_full + fit$source_weights[s] * sources[[s]]
  }
  expect_equal(
    fit$covariance,
    (1 - fit$selected_weight) * target_full +
      fit$selected_weight * source_full
  )
  expect_equal(fit$covtl_covariance, cov_fit$covariance)
  expect_equal(fit$covtl_source_weights, cov_fit$source_weights)
  expect_equal(
    unname(fit$covtl_source_weights / fit$covtl_transfer_weight),
    unname(fit$source_weights),
    tolerance = 1e-10
  )
  expect_equal(sum(diag(fit$projector)), 1, tolerance = 1e-10)
  expect_equal(
    fit$projector %*% fit$projector,
    fit$projector,
    tolerance = 1e-10
  )
  expect_null(fit$target_covariance)
  expect_null(fit$source_covariance)
  expect_null(fit$source_covariances)
  expect_null(fit$sources)
  expect_null(fit$n_target)
  expect_null(fit$n_source)
  expect_null(fit$center)
})

test_that("multi-source EigenTL one-SE and min share held-out paths", {
  X <- rbind(
    c(2.0, 0.1), c(-1.8, -0.2), c(1.1, 0.8),
    c(-0.9, -1.0), c(0.4, -0.7), c(-0.5, 0.6)
  )
  sources <- list(
    EUR = diag(c(3.0, 0.5)),
    AFR = diag(c(0.5, 3.0)),
    AMR = matrix(c(2, 0.4, 0.4, 1.5), 2)
  )
  fold_id <- c(1, 2, 3, 1, 2, 3)
  weights <- c(0, 0.25, 0.5, 1)

  set.seed(12)
  min_fit <- multisource_eigen_tl(
    X, sources, rank = 1, transfer_weight_grid = weights,
    center = FALSE, fold_id = fold_id, method = "min"
  )
  set.seed(12)
  one_se <- multisource_eigen_tl(
    X, sources, rank = 1, transfer_weight_grid = weights,
    center = FALSE, fold_id = fold_id, method = "one_se"
  )

  expect_equal(one_se$observation_scores, min_fit$observation_scores)
  expect_equal(one_se$fold_source_weights, min_fit$fold_source_weights)
  expect_equal(one_se$best_index, min_fit$selected_index)
  expect_gte(one_se$selected_index, one_se$best_index)
  expect_equal(one_se$selected_weight, weights[one_se$selected_index])
  expect_error(
    multisource_eigen_tl(
      X, sources, rank = 1, transfer_weight_grid = c(0, 0.5),
      fold_id = fold_id
    ),
    "include endpoints"
  )
})
