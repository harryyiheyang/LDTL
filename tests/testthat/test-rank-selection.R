test_that("projected bootstrap returns a high-rank cumulative-variance path", {
  set.seed(101)
  n <- 160
  p <- 40
  # A deliberately non-spiked, slowly decaying spectrum. The selector should
  # not collapse this design to a conventional 3--5 factor rank.
  population_values <- seq(2, 0.6, length.out = p)
  X <- matrix(rnorm(n * p), n, p) %*% diag(sqrt(population_values))

  fit <- select_tl_rank_bootstrap(
    X,
    n_boot = 12,
    initial_rank = 12,
    rank_block = 8,
    oversample = 6,
    bootstrap_extra = 6,
    power_iterations = 1,
    n_threads = 2,
    seed = 44
  )

  expect_s3_class(fit, "tl_rank_bootstrap")
  expect_gt(fit$rank, 5)
  expect_equal(fit$rank_path$target_fraction, c(.5, .6, .7, .8, .9))
  expect_true(all(diff(fit$rank_path$rank) >= 0))
  expect_true(all(fit$rank_path$rank <= fit$reference_rank))
  expect_equal(
    dim(fit$bootstrap_cumulative_fraction),
    c(12, fit$flexible_keep)
  )
  expect_equal(
    dim(fit$bootstrap_flexible_subspace_loss),
    c(12, fit$flexible_keep)
  )
  expect_true(all(is.finite(fit$bootstrap_cumulative_fraction)))
  expect_true(all(fit$bootstrap_cumulative_fraction >= 0))
  expect_true(all(fit$bootstrap_cumulative_fraction <= 1 + 1e-10))
  expect_true(all(apply(
    fit$bootstrap_cumulative_fraction,
    1,
    function(value) all(diff(value) >= -1e-12)
  )))
  reference_rank_80 <- which(
    cumsum(fit$eigenvalues) / fit$total_variance >= .8
  )[1]
  expect_gte(fit$rank, reference_rank_80)
  expect_false(fit$source_data_used)
  expect_true(fit$conditional_on_fixed_core)
})

test_that("projected bootstrap is reproducible across OpenMP thread counts", {
  set.seed(202)
  X <- matrix(rnorm(100 * 24), 100, 24)
  arguments <- list(
    X_target = X,
    target_fraction_grid = c(.5, .7, .8),
    n_boot = 10,
    initial_rank = 10,
    rank_block = 6,
    oversample = 5,
    bootstrap_extra = 5,
    power_iterations = 1,
    seed = 812
  )
  one_thread <- do.call(
    select_tl_rank_bootstrap,
    c(arguments, list(n_threads = 1))
  )
  two_threads <- do.call(
    select_tl_rank_bootstrap,
    c(arguments, list(n_threads = 2))
  )

  expect_equal(one_thread$rank_path, two_threads$rank_path, tolerance = 1e-12)
  expect_equal(
    one_thread$bootstrap_cumulative_fraction,
    two_threads$bootstrap_cumulative_fraction,
    tolerance = 1e-12
  )
  expect_equal(
    one_thread$bootstrap_flexible_subspace_loss,
    two_threads$bootstrap_flexible_subspace_loss,
    tolerance = 1e-12
  )
})

test_that("rank selection validates operational fractions and controls", {
  X <- matrix(rnorm(120), 20, 6)
  expect_error(
    select_tl_rank_bootstrap(X, target_fraction_grid = c(.4, .8)),
    "anchor_fraction"
  )
  expect_error(
    select_tl_rank_bootstrap(X, guard_fraction = .5),
    "guard_fraction"
  )
  expect_error(
    select_tl_rank_bootstrap(X, selection_fraction = 1),
    "selection_fraction"
  )
  expect_error(
    select_tl_rank_bootstrap(X, n_boot = 1),
    "n_boot"
  )
  expect_error(
    select_tl_rank_bootstrap(X, seed = -1),
    "seed"
  )
})
