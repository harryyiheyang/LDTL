ar_correlation <- function(p, rho) {
  stats::toeplitz(rho^(0:(p - 1L)))
}

eigen_pair <- function(R, rank = nrow(R)) {
  eig <- eigen(R, symmetric = TRUE)
  keep <- seq_len(rank)
  list(U = eig$vectors[, keep, drop = FALSE], D = eig$values[keep])
}

test_that("CovTL_SS fits shared summary-statistic source weights", {
  set.seed(20260826)
  p <- 6L
  L <- 12L
  true_weights <- c(`1000G` = 0.20, UKB = 0.45, AoU = 0.34, Identity = 0.01)
  z <- vector("list", L)
  sources <- vector("list", L)

  for (l in seq_len(L)) {
    shift <- seq(-0.01, 0.01, length.out = L)[l]
    sources[[l]] <- list(
      `1000G` = ar_correlation(p, 0.50 + shift),
      UKB = ar_correlation(p, 0.48 + shift),
      AoU = ar_correlation(p, 0.51 + shift)
    )
    target <- true_weights[1L] * sources[[l]][[1L]] +
      true_weights[2L] * sources[[l]][[2L]] +
      true_weights[3L] * sources[[l]][[3L]] +
      true_weights[4L] * diag(p)
    z[[l]] <- as.numeric(t(chol(target)) %*% rnorm(p))
  }

  fit <- CovTL_SS(z, sources)

  expect_s3_class(fit, "CovTL_SS")
  expect_equal(sum(fit$weights), 1, tolerance = 1e-10)
  expect_true(all(fit$weights >= 0))
  expect_lte(fit$identity_weight, 0.015 + 1e-10)
  expect_equal(names(fit$source_weights), c("1000G", "UKB", "AoU"))
  expect_length(fit$covariance, L)
  expect_true(all(vapply(fit$covariance, is.matrix, logical(1L))))
  expect_true(is.finite(fit$negative_log_likelihood))
  expect_true(all(is.finite(fit$fisher_information)))
  expect_equal(dim(fit$fisher_information), c(3L, 3L))
  expect_lte(fit$iterations, 20L)
  expect_lte(fit$line_search_evaluations, 3L * fit$iterations)

  expected <- fit$identity_weight * diag(p)
  for (j in seq_along(fit$source_weights)) {
    expected <- expected + fit$source_weights[j] * sources[[1L]][[j]]
  }
  expect_equal(fit$covariance[[1L]], expected, tolerance = 1e-10)
})

test_that("CovTL_SS accepts full-rank and truncated U-D source pairs", {
  set.seed(17)
  p <- 8L
  R1 <- ar_correlation(p, 0.50)
  R2 <- ar_correlation(p, 0.48)
  R3 <- ar_correlation(p, 0.51)
  z <- rnorm(p)

  matrix_fit <- CovTL_SS(
    z,
    list(`1000G` = R1, UKB = R2, AoU = R3)
  )
  pair_fit <- CovTL_SS(
    z,
    list(
      `1000G` = eigen_pair(R1),
      UKB = eigen_pair(R2),
      AoU = eigen_pair(R3)
    )
  )
  truncated_fit <- CovTL_SS(
    z,
    list(
      truncated = eigen_pair(R1, rank = 3L),
      full = eigen_pair(R2)
    )
  )

  expect_equal(pair_fit$weights, matrix_fit$weights, tolerance = 1e-8)
  expect_equal(
    pair_fit$negative_log_likelihood,
    matrix_fit$negative_log_likelihood,
    tolerance = 1e-8
  )
  expect_s3_class(truncated_fit, "CovTL_SS")
  expect_true(is.finite(truncated_fit$negative_log_likelihood))
  expect_equal(sum(truncated_fit$weights), 1, tolerance = 1e-10)
})

test_that("CovTL_SS validates locus dimensions, source names, and controls", {
  z <- list(a = c(0.1, -0.2, 0.3), b = c(0.2, 0.1, -0.1))
  R <- ar_correlation(3L, 0.5)

  expect_error(
    CovTL_SS(z, list(list(EUR = R), list(AFR = R))),
    "same named sources"
  )
  expect_error(
    CovTL_SS(c(0.1, -0.2), list(EUR = diag(3))),
    "p by p"
  )
  bad <- R
  bad[1, 1] <- 0
  expect_error(
    CovTL_SS(c(0.1, -0.2, 0.3), list(EUR = bad)),
    "positive finite diagonal"
  )
  expect_error(
    CovTL_SS(c(0.1, -0.2, 0.3), list(EUR = R), identity_max = 1),
    "permitted range"
  )
  expect_error(
    CovTL_SS(c(0.1, -0.2, 0.3), list(EUR = R), max_iterations = 0),
    "positive integer"
  )
})
