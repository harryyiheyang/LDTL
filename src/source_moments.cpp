#include <Rcpp.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

int resolve_threads(const int requested) {
#ifdef _OPENMP
  const int available = std::max(1, omp_get_max_threads());
  if (requested <= 0) {
    return available;
  }
  return std::max(1, std::min(requested, available));
#else
  (void)requested;
  return 1;
#endif
}

template <int RTYPE>
inline double matrix_value(
    const typename Rcpp::Matrix<RTYPE>& x,
    const R_xlen_t index
) {
  return static_cast<double>(x[index]);
}

template <>
inline double matrix_value<INTSXP>(
    const Rcpp::Matrix<INTSXP>& x,
    const R_xlen_t index
) {
  const int value = x[index];
  return value == NA_INTEGER ? NA_REAL : static_cast<double>(value);
}

template <int RTYPE>
List source_moments_impl(
    const typename Rcpp::Matrix<RTYPE>& x,
    const bool center,
    const int block_size,
    const int requested_threads,
    const bool compute_covariance
) {
  const R_xlen_t n = x.nrow();
  const R_xlen_t p = x.ncol();
  const int threads = resolve_threads(requested_threads);
  NumericVector means(p, 0.0);
  std::atomic<bool> invalid(false);

  // Validate the input and compute column means without constructing a
  // centered copy of the biobank-scale matrix.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
  for (R_xlen_t j = 0; j < p; ++j) {
    long double total = 0.0L;
    const R_xlen_t offset = n * j;
    for (R_xlen_t i = 0; i < n; ++i) {
      const double value = matrix_value<RTYPE>(x, offset + i);
      if (!R_FINITE(value)) {
        invalid.store(true, std::memory_order_relaxed);
      } else if (center) {
        total += static_cast<long double>(value);
      }
    }
    if (center) {
      means[j] = static_cast<double>(total / static_cast<long double>(n));
    }
  }
  if (invalid.load(std::memory_order_relaxed)) {
    stop("X_source must contain only finite values; impute missing genotypes before computing moments.");
  }

  NumericMatrix covariance;
  if (compute_covariance) {
    covariance = NumericMatrix(p, p);
  }

  long double fourth_sum = 0.0L;
  const R_xlen_t step = std::max<R_xlen_t>(1, block_size);

  // Individual blocking bounds the working set and reduces the number of
  // OpenMP regions for large cohorts. No n-by-p copy or individual outer
  // product is formed.
  for (R_xlen_t start = 0; start < n; start += step) {
    const R_xlen_t end = std::min(n, start + step);
    const R_xlen_t block_n = end - start;
    std::vector<long double> row_norm_squared(block_n, 0.0L);

    // R matrices are column-major. Accumulating contiguous pieces of columns
    // into thread-local row-norm buffers avoids the cache-unfriendly pattern
    // of walking all columns separately for every individual.
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
    {
      std::vector<long double> local_norm(block_n, 0.0L);
#pragma omp for schedule(static)
      for (R_xlen_t j = 0; j < p; ++j) {
        const R_xlen_t offset = n * j;
        for (R_xlen_t i = start; i < end; ++i) {
          const long double value = static_cast<long double>(
              matrix_value<RTYPE>(x, offset + i) - means[j]
          );
          local_norm[i - start] += value * value;
        }
      }
#pragma omp critical
      {
        for (R_xlen_t i = 0; i < block_n; ++i) {
          row_norm_squared[i] += local_norm[i];
        }
      }
    }
#else
    for (R_xlen_t j = 0; j < p; ++j) {
      const R_xlen_t offset = n * j;
      for (R_xlen_t i = start; i < end; ++i) {
        const long double value = static_cast<long double>(
            matrix_value<RTYPE>(x, offset + i) - means[j]
        );
        row_norm_squared[i - start] += value * value;
      }
    }
#endif

    long double block_fourth_sum = 0.0L;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : block_fourth_sum) num_threads(threads)
#endif
    for (R_xlen_t i = 0; i < block_n; ++i) {
      block_fourth_sum += row_norm_squared[i] * row_norm_squared[i];
    }
    fourth_sum += block_fourth_sum;

    if (compute_covariance) {
      // Parallelize over columns of the upper triangle. Each worker owns all
      // updates for its (j, k) entries, so no locks or thread-local p-by-p
      // matrices are required.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
      for (R_xlen_t j = 0; j < p; ++j) {
        for (R_xlen_t k = j; k < p; ++k) {
          long double cross = 0.0L;
          const R_xlen_t offset_j = n * j;
          const R_xlen_t offset_k = n * k;
          for (R_xlen_t i = start; i < end; ++i) {
            const long double left = static_cast<long double>(
                matrix_value<RTYPE>(x, offset_j + i) - means[j]
            );
            const long double right = static_cast<long double>(
                matrix_value<RTYPE>(x, offset_k + i) - means[k]
            );
            cross += left * right;
          }
          covariance[j + p * k] += static_cast<double>(cross);
        }
      }
    }

    checkUserInterrupt();
  }

  double variance_raw = NA_REAL;
  if (compute_covariance) {
    const long double n_long = static_cast<long double>(n);
    long double frobenius_squared = 0.0L;
    for (R_xlen_t j = 0; j < p; ++j) {
      for (R_xlen_t k = j; k < p; ++k) {
        const double value = covariance[j + p * k] / static_cast<double>(n);
        covariance[j + p * k] = value;
        covariance[k + p * j] = value;
        const long double square =
            static_cast<long double>(value) * static_cast<long double>(value);
        frobenius_squared += (j == k) ? square : 2.0L * square;
      }
    }
    const long double denominator = n_long * (n_long - 1.0L);
    variance_raw = static_cast<double>(
        (fourth_sum - n_long * frobenius_squared) / denominator
    );
  }

  return List::create(
      _["covariance"] = compute_covariance ? static_cast<SEXP>(covariance) : R_NilValue,
      _["fourth_sum"] = static_cast<double>(fourth_sum),
      _["variance_raw"] = variance_raw,
      _["mean"] = means,
      _["n"] = static_cast<double>(n),
      _["p"] = static_cast<double>(p),
      _["block_size"] = block_size,
      _["threads_used"] = threads,
      _["openmp"] = static_cast<bool>(
#ifdef _OPENMP
          true
#else
          false
#endif
      )
  );
}

}  // namespace

// [[Rcpp::export]]
List cpp_source_moments(
    SEXP X_source,
    const bool center,
    const int block_size,
    const int n_threads,
    const bool compute_covariance
) {
  if (!Rf_isMatrix(X_source)) {
    stop("X_source must be a numeric or integer matrix with individuals in rows.");
  }
  const IntegerVector dimensions = Rf_getAttrib(X_source, R_DimSymbol);
  if (dimensions[0] < 2 || dimensions[1] < 1) {
    stop("X_source must have at least 2 rows and 1 column.");
  }
  if (block_size < 1) {
    stop("block_size must be a positive integer.");
  }

  switch (TYPEOF(X_source)) {
    case REALSXP:
      return source_moments_impl<REALSXP>(
          NumericMatrix(X_source), center, block_size, n_threads,
          compute_covariance
      );
    case INTSXP:
      return source_moments_impl<INTSXP>(
          IntegerMatrix(X_source), center, block_size, n_threads,
          compute_covariance
      );
    default:
      stop("X_source must be a numeric or integer matrix.");
  }
  return List::create();
}
