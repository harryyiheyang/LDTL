test_that("POET preserves small nonzero LD entries before factor decomposition", {
  S <- diag(8)
  S[1, 2] <- 5e-5
  S[2, 1] <- 5e-5

  comp <- LDTL:::.ld_poet_components(
    S,
    n = 200,
    cutoff_method = "D.ratio",
    k_min = 2,
    k_max = 4
  )

  expect_equal(comp$S, S, tolerance = 1e-14)
})

test_that("POET residual diagonal uses the closest positive diagonal value", {
  E <- diag(c(-1, 0, 0.2, 0.5))
  fixed <- LDTL:::.ld_fix_residual_diag(
    E,
    fallback_diag = rep(1, 4)
  )

  expect_equal(diag(fixed), c(0.2, 0.2, 0.2, 0.5))

  no_positive <- diag(c(-2, -1, 0))
  fallback <- LDTL:::.ld_fix_residual_diag(
    no_positive,
    fallback_diag = c(1, 2, 3)
  )
  expect_equal(diag(fallback), rep(1, 3))
})

test_that("POET applies the S-POET positive-part spike correction", {
  d <- c(2, 1.5, 1, 0.7, 0.5, 0.3)
  correction <- LDTL:::.ld_poet_spike_correction(
    d,
    factors = 2,
    n = 100,
    trace_S = 6
  )

  denom <- 6 - 2 - 2 * 6 / 100
  c_hat <- (6 - 2 - 1.5) / denom
  bias <- c_hat * 6 / 100
  expect_equal(correction$c_hat, c_hat, tolerance = 1e-14)
  expect_equal(correction$bias, bias, tolerance = 1e-14)
  expect_equal(correction$values, pmax(d[1:2] - bias, 0))
  expect_error(
    LDTL:::.ld_poet_spike_correction(
      d,
      factors = 5,
      n = 3,
      trace_S = 6
    ),
    "p - K"
  )
})

test_that("POET reuses a supplied full eigendecomposition", {
  set.seed(8)
  X <- matrix(stats::rnorm(4000), nrow = 200, ncol = 20)
  S <- stats::cor(X)
  eig <- LDTL:::.ld_eigen(S)

  computed <- poet_thresholding(S = S, n = 200)
  supplied <- poet_thresholding(S = S, n = 200, eig = eig)

  expect_equal(supplied, computed, tolerance = 1e-12)
  expect_error(
    poet_thresholding(
      S = S,
      n = 200,
      eig = list(values = eig$values[1:5], vectors = eig$vectors[, 1:5])
    ),
    "p eigenvalues"
  )
})

test_that("POET accepts a fixed common factor rank", {
  set.seed(18)
  X <- matrix(stats::rnorm(6000), nrow = 300, ncol = 20)
  S <- stats::cor(X)
  eig <- LDTL:::.ld_eigen(S)
  comp <- LDTL:::.ld_poet_components(
    S,
    n = 300,
    cutoff_method = "ACT",
    k_min = 5,
    k_max = NULL,
    eig = eig,
    factors = 7
  )

  expect_equal(comp$factors, 7L)
  expect_equal(
    poet_thresholding(S = S, n = 300, eig = eig, factors = 7),
    poet_thresholding(S = S, n = 300, eig = eig, factors = 7,
                      factor_method = "D.ratio"),
    tolerance = 1e-12
  )
  expect_error(
    poet_thresholding(S = S, n = 300, factors = 20),
    "between 1 and nrow"
  )
})

test_that("POET linear shrinkage uses the MSE plug-in intensity", {
  E <- matrix(
    c(
      1, 0.2, 0.1,
      0.2, 1, 0.05,
      0.1, 0.05, 1
    ),
    nrow = 3,
    byrow = TRUE
  )
  off_squared <- sum(E^2) - sum(diag(E)^2)
  diagonal_products <- sum(diag(E))^2 - sum(diag(E)^2)
  expected <- (off_squared + diagonal_products) / 9999 / off_squared

  alpha_mse <- LDTL:::.ld_mse_shrinkage_intensity(E, n = 10000)

  expect_equal(alpha_mse, expected, tolerance = 1e-14)
  expect_lt(alpha_mse, 0.05)
  expect_identical(formals(poet_linear_shrinkage)$alpha, 0.05)
})

test_that("ACT uses adjusted correlation eigenvalues", {
  d <- c(10, 8, 6, 5, 4, 3, 2, 1, 0.8, 0.6, 0.4, 0.2)

  adjusted <- LDTL:::.ld_act_adjusted_eigenvalues(
    d,
    n = 100,
    k_max = 8
  )
  factors <- LDTL:::.ld_factor_count_from_values(
    d,
    k_min = 5,
    k_max = 8,
    cutoff_method = "ACT",
    n = 100
  )

  expect_equal(adjusted[1], 7.784903258518539, tolerance = 1e-12)
  expect_equal(factors, 7L)
})

test_that("ACT retains the established POET lower factor bound", {
  d <- rep(1, 12)

  factors <- LDTL:::.ld_factor_count_from_values(
    d,
    k_min = 5,
    k_max = 8,
    cutoff_method = "ACT",
    n = 100
  )

  expect_equal(factors, 5L)
})

test_that("POET factor-search upper bound uses 90 percent eigenvalue mass", {
  flat <- rep(1, 50)
  dominant <- c(1000, rep(1, 49))

  expect_equal(
    LDTL:::.ld_poet_factor_max(flat, k_min = 5),
    45L
  )
  expect_equal(
    LDTL:::.ld_poet_factor_max(dominant, k_min = 5),
    5L
  )
})

test_that("D.ratio remains available without changing its selection rule", {
  d <- c(10, 8, 6, 5, 4, 3, 2, 1, 0.8, 0.6, 0.4, 0.2)

  factors <- LDTL:::.ld_factor_count_from_values(
    d,
    k_min = 5,
    k_max = 8,
    cutoff_method = "D.ratio"
  )

  expect_equal(factors, 8L)
})

test_that("ACT is the primary POET factor-selection method", {
  expect_identical(
    formals(poet_thresholding)$factor_method,
    quote(c("ACT", "D.ratio"))
  )
  expect_identical(
    formals(poet_tapering)$factor_method,
    quote(c("ACT", "D.ratio"))
  )
  expect_identical(
    formals(poet_linear_shrinkage)$factor_method,
    quote(c("ACT", "D.ratio"))
  )
  expect_identical(
    formals(poet_banding)$factor_method,
    quote(c("ACT", "D.ratio"))
  )
  expect_identical(
    formals(poet_nonlinear_shrinkage)$factor_method,
    quote(c("ACT", "D.ratio"))
  )
  expect_null(formals(poet_thresholding)$eig)
  expect_null(formals(poet_tapering)$eig)
  expect_null(formals(poet_linear_shrinkage)$eig)
  expect_null(formals(poet_banding)$eig)
  expect_null(formals(poet_nonlinear_shrinkage)$eig)
})

test_that("POET runs with both ACT and D.ratio", {
  set.seed(21)
  X <- matrix(stats::rnorm(4000), nrow = 200, ncol = 20)
  S <- stats::cor(X)

  R_act <- poet_thresholding(S = S, n = 200)
  R_dratio <- poet_thresholding(
    S = S,
    n = 200,
    factor_method = "D.ratio"
  )

  expect_equal(dim(R_act), c(20L, 20L))
  expect_equal(R_act, t(R_act), tolerance = 1e-12)
  expect_equal(diag(R_act), rep(1, 20), tolerance = 1e-12)
  expect_true(all(is.finite(R_act)))
  expect_equal(dim(R_dratio), c(20L, 20L))
  expect_true(all(is.finite(R_dratio)))
})
