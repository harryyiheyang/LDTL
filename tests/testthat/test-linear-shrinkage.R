test_that("linear shrinkage uses the MSE plug-in rule below its boundary", {
  S <- matrix(
    c(
      1, 0.2, 0.1,
      0.2, 1, 0.05,
      0.1, 0.05, 1
    ),
    nrow = 3,
    byrow = TRUE
  )
  target <- diag(3)
  alpha_mse <- LDRegularization:::.ld_mse_shrinkage_intensity(
    S,
    n = 10000,
    target = target
  )

  fit <- linear_shrinkage(S = S, n = 10000)
  expected <- (1 - alpha_mse) * S + alpha_mse * target

  expect_lt(alpha_mse, 0.05)
  expect_equal(as.numeric(fit), as.numeric(expected), tolerance = 1e-12)
  expect_equal(dim(fit), dim(expected))
})

test_that("linear shrinkage respects a smaller supplied boundary", {
  S <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  fit <- linear_shrinkage(S = S, n = 100, alpha = 0.001)
  expected <- 0.999 * S + 0.001 * diag(2)

  expect_equal(as.numeric(fit), as.numeric(expected), tolerance = 1e-12)
  expect_equal(dim(fit), dim(expected))
  expect_error(linear_shrinkage(S = S), "n must")
})
