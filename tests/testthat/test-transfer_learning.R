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

  fit <- cov_tl(X, S_source)

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

test_that("eigspac_tl uses the analytic first-order tangent-risk weight", {
  X <- cbind(
    c(2, -2, 1, -1),
    c(1, 1, -1, -1)
  )
  S_target <- crossprod(X) / nrow(X)
  S_source <- matrix(c(2.5, 0.2, 0.2, 1), 2, 2)
  gap <- S_target[1, 1] - S_target[2, 2]
  cross_scores <- X[, 1] * X[, 2]
  tangent_variance <-
    2 * sum(cross_scores^2) / (nrow(X) * (nrow(X) - 1) * gap^2)
  tangent_distance <- 2 * (S_source[1, 2] / gap)^2
  expected_lambda <- min(tangent_variance / tangent_distance, 1)

  fit <- eigspac_tl(X, S_source, rank = 1)

  expect_s3_class(fit, "eigspac_tl")
  expect_equal(fit$tangent_variance, tangent_variance)
  expect_equal(fit$tangent_distance_squared, tangent_distance)
  expect_equal(fit$lambda, expected_lambda)
  expect_equal(fit$projector, tcrossprod(fit$vectors))
  expect_equal(sum(diag(fit$projector)), 1)
  expect_false(fit$conditional_ure_if_pilot_independent)
})

test_that("eigspac_tl supports an independent pilot derivative", {
  X <- cbind(
    c(2, -2, 1, -1),
    c(1, 1, -1, -1)
  )
  S_source <- matrix(c(2.5, 0.2, 0.2, 1), 2, 2)
  S_pilot <- diag(c(3, 1))

  fit <- eigspac_tl(X, S_source, rank = 1, S_pilot = S_pilot)

  expect_true(fit$pilot_supplied)
  expect_true(fit$conditional_ure_if_pilot_independent)
  expect_equal(fit$pilot_eigengap, 2)
})

test_that("transfer-learning inputs are validated", {
  X <- matrix(rnorm(20), 10, 2)

  expect_error(cov_tl(X, diag(3)), "same number of variables")
  expect_error(eigspac_tl(X, diag(2), rank = 0), "rank must be between")
  expect_error(eigspac_tl(X, diag(2), rank = 2), "rank must be between")
  expect_error(
    eigspac_tl(X, diag(2), rank = 1, S_pilot = diag(2)),
    "eigengap"
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
