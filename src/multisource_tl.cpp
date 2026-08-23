#include <Rcpp.h>

#include <cfloat>
#include <cmath>
#include <cstddef>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cpp_multi_source_gram(
    const NumericMatrix target,
    const List sources
) {
  const R_xlen_t p = target.nrow();
  if (p < 1 || target.ncol() != p) {
    stop("target must be a nonempty square numeric matrix.");
  }

  const R_xlen_t n_sources = sources.size();
  if (n_sources < 2) {
    stop("sources must contain at least two covariance matrices.");
  }

  std::vector<NumericMatrix> source_matrices;
  source_matrices.reserve(static_cast<std::size_t>(n_sources));
  for (R_xlen_t source_index = 0; source_index < n_sources; ++source_index) {
    SEXP source_sexp = sources[source_index];
    if (!Rf_isMatrix(source_sexp) || TYPEOF(source_sexp) != REALSXP) {
      stop("Every source covariance must be a numeric matrix.");
    }
    NumericMatrix source(source_sexp);
    if (source.nrow() != p || source.ncol() != p) {
      stop("Every source covariance must have the same dimensions as target.");
    }
    source_matrices.push_back(source);
  }

  const R_xlen_t entries = p * p;
  std::vector<long double> gram(
    static_cast<std::size_t>(n_sources * n_sources),
    0.0L
  );
  std::vector<long double> differences(
    static_cast<std::size_t>(n_sources),
    0.0L
  );

  for (R_xlen_t entry = 0; entry < entries; ++entry) {
    const double target_value = target[entry];
    if (!R_FINITE(target_value)) {
      stop("target must contain only finite values.");
    }

    for (R_xlen_t source_index = 0;
         source_index < n_sources;
         ++source_index) {
      const double source_value = source_matrices[source_index][entry];
      if (!R_FINITE(source_value)) {
        stop("Every source covariance must contain only finite values.");
      }
      differences[source_index] =
        static_cast<long double>(source_value) -
        static_cast<long double>(target_value);
    }

    for (R_xlen_t left = 0; left < n_sources; ++left) {
      for (R_xlen_t right = 0; right <= left; ++right) {
        gram[left + n_sources * right] +=
          differences[left] * differences[right];
      }
    }

    if ((entry & 1048575) == 0) {
      checkUserInterrupt();
    }
  }

  NumericMatrix output(n_sources, n_sources);
  for (R_xlen_t left = 0; left < n_sources; ++left) {
    for (R_xlen_t right = 0; right <= left; ++right) {
      const double value = static_cast<double>(
        gram[left + n_sources * right]
      );
      output(left, right) = value;
      output(right, left) = value;
    }
  }
  return output;
}

// [[Rcpp::export]]
NumericMatrix cpp_multi_source_combine(
    const NumericMatrix target,
    const List sources,
    const NumericVector source_weights
) {
  const R_xlen_t p = target.nrow();
  const R_xlen_t n_sources = sources.size();
  if (p < 1 || target.ncol() != p) {
    stop("target must be a nonempty square numeric matrix.");
  }
  if (n_sources < 2 || source_weights.size() != n_sources) {
    stop("source_weights must contain one value per source covariance.");
  }

  std::vector<NumericMatrix> source_matrices;
  source_matrices.reserve(static_cast<std::size_t>(n_sources));
  long double source_weight_sum = 0.0L;
  for (R_xlen_t source_index = 0; source_index < n_sources; ++source_index) {
    SEXP source_sexp = sources[source_index];
    if (!Rf_isMatrix(source_sexp) || TYPEOF(source_sexp) != REALSXP) {
      stop("Every source covariance must be a numeric matrix.");
    }
    NumericMatrix source(source_sexp);
    if (source.nrow() != p || source.ncol() != p) {
      stop("Every source covariance must have the same dimensions as target.");
    }
    const double weight = source_weights[source_index];
    if (!R_FINITE(weight) || weight < 0.0) {
      stop("source_weights must be nonnegative and finite.");
    }
    source_weight_sum += static_cast<long double>(weight);
    source_matrices.push_back(source);
  }
  if (source_weight_sum > 1.0L + 100.0L * DBL_EPSILON) {
    stop("source_weights must sum to at most one.");
  }
  const long double target_weight = 1.0L - source_weight_sum;

  NumericMatrix output(p, p);
  for (R_xlen_t column = 0; column < p; ++column) {
    for (R_xlen_t row = 0; row <= column; ++row) {
      const R_xlen_t upper = row + p * column;
      const R_xlen_t lower = column + p * row;
      long double upper_value =
        target_weight * static_cast<long double>(target[upper]);
      long double lower_value =
        target_weight * static_cast<long double>(target[lower]);
      for (R_xlen_t source_index = 0;
           source_index < n_sources;
           ++source_index) {
        const long double weight = static_cast<long double>(
          source_weights[source_index]
        );
        upper_value += weight * source_matrices[source_index][upper];
        lower_value += weight * source_matrices[source_index][lower];
      }
      const double value = static_cast<double>(
        (upper_value + lower_value) / 2.0L
      );
      output[upper] = value;
      output[lower] = value;
    }
    if ((column & 255) == 0) {
      checkUserInterrupt();
    }
  }
  return output;
}
