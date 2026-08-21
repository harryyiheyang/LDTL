write_plink_raw_fixture <- function(X, variants) {
  file <- tempfile(fileext = ".raw")
  header <- c(
    "FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE",
    paste0(variants, "_A")
  )
  lines <- character(nrow(X) + 1L)
  lines[1L] <- paste(header, collapse = "\t")
  for (i in seq_len(nrow(X))) {
    dosage <- ifelse(is.na(X[i, ]), "NA", format(X[i, ], trim = TRUE))
    lines[i + 1L] <- paste(
      c(i, i, 0, 0, 0, -9, dosage),
      collapse = "\t"
    )
  }
  writeLines(lines, file)
  file
}

test_that("source_moments_raw matches the in-memory source moments", {
  X <- matrix(
    c(
      0, 1, 2,
      1, 1, NA,
      2, 0, 1,
      0, NA, 2,
      1, 2, 0,
      2, 1, 1
    ),
    nrow = 6L,
    byrow = TRUE
  )
  variants <- c("rs_one", "rs_two", "rs_three")
  file <- write_plink_raw_fixture(X, variants)
  on.exit(unlink(file))

  mu <- colMeans(X, na.rm = TRUE)
  for (j in seq_len(ncol(X))) {
    X[is.na(X[, j]), j] <- mu[j]
  }
  X <- sweep(X, 2L, mu, "-")
  scale <- sqrt(colSums(X^2) / nrow(X))
  X <- sweep(X, 2L, scale, "/")

  streamed <- source_moments_raw(
    file,
    block_size = 2L,
    n_threads = 1L
  )
  direct <- source_moments_method(X, center = FALSE, n_threads = 1L)

  expect_s3_class(streamed, "ld_source_moments")
  expect_equal(streamed$variants, variants)
  expect_equal(streamed$mean, mu, tolerance = 1e-14)
  expect_equal(streamed$scale, scale, tolerance = 1e-14)
  expect_equal(streamed$observed, c(6, 5, 5))
  expect_equal(streamed$missing_calls, 2)
  expect_equal(streamed$covariance, direct$covariance, tolerance = 1e-13)
  expect_equal(streamed$fourth_sum, direct$fourth_sum, tolerance = 1e-12)
  expect_equal(streamed$variance, direct$variance, tolerance = 1e-12)
  expect_equal(diag(streamed$covariance), rep(1, 3), tolerance = 1e-13)
  expect_equal(streamed$passes, 2)
  expect_lt(streamed$working_bytes, 1000)

  cached <- source_moments_raw(
    file,
    S_source = direct$covariance,
    block_size = 2L,
    n_threads = 1L
  )
  expect_false(cached$covariance_computed)
  expect_equal(cached$covariance, direct$covariance, tolerance = 1e-14)
  expect_equal(cached$fourth_sum, direct$fourth_sum, tolerance = 1e-12)
  expect_equal(cached$variance, direct$variance, tolerance = 1e-12)
  expect_lt(cached$working_bytes, streamed$working_bytes)
  expect_error(
    source_moments_raw(file, S_source = diag(2)),
    "same number of variables"
  )

  target <- matrix(rnorm(24), 8, 3)
  fit <- cov_tl(X_target = target, source = streamed, center = TRUE)
  expect_true(all(is.finite(fit$covariance)))
})

test_that("source_moments_raw validates PLINK raw rows", {
  file <- tempfile(fileext = ".raw")
  on.exit(unlink(file))
  writeLines(c(
    "FID IID PAT MAT SEX PHENOTYPE v1_A v2_A",
    "1 1 0 0 0 -9 0"
  ), file)

  expect_error(source_moments_raw(file), "too few dosage fields")
  expect_error(source_moments_raw(file, block_size = 0), "positive integer")
  expect_error(source_moments_raw(file, n_threads = -1), "zero or a positive")
})

test_that("source_moments_raw rejects a monomorphic variant", {
  X <- cbind(c(0, 1, 2), c(1, 1, 1))
  file <- write_plink_raw_fixture(X, c("v1", "v2"))
  on.exit(unlink(file))

  expect_error(source_moments_raw(file), "nonpositive variance")
})
