// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

struct ApproximateEigenspace {
  arma::vec values;
  arma::mat vectors;
  bool success;
};

int rank_threads(const int requested) {
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

std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

arma::mat normal_matrix(
    const arma::uword rows,
    const arma::uword columns,
    std::mt19937_64& generator
) {
  std::normal_distribution<double> normal(0.0, 1.0);
  arma::mat result(rows, columns);
  for (arma::uword index = 0; index < result.n_elem; ++index) {
    result[index] = normal(generator);
  }
  return result;
}

bool orthonormal_basis(const arma::mat& input, arma::mat& basis) {
  arma::mat triangular;
  return arma::qr_econ(basis, triangular, input) &&
    basis.n_cols > 0 && basis.is_finite();
}

ApproximateEigenspace randomized_eigenspace(
    const arma::mat& scaled_data,
    const int keep,
    const int oversample,
    const int power_iterations,
    std::mt19937_64& generator
) {
  const int maximum = static_cast<int>(
    std::min(scaled_data.n_rows, scaled_data.n_cols)
  );
  const int retained = std::max(1, std::min(keep, maximum));
  const int sketch_size = std::min(
    maximum,
    retained + std::max(0, oversample)
  );
  arma::mat omega = normal_matrix(
    scaled_data.n_cols,
    static_cast<arma::uword>(sketch_size),
    generator
  );
  arma::mat left;
  if (!orthonormal_basis(scaled_data * omega, left)) {
    return {arma::vec(), arma::mat(), false};
  }
  arma::mat right;
  for (int iteration = 0; iteration <= power_iterations; ++iteration) {
    if (!orthonormal_basis(scaled_data.t() * left, right)) {
      return {arma::vec(), arma::mat(), false};
    }
    if (iteration < power_iterations) {
      if (!orthonormal_basis(scaled_data * right, left)) {
        return {arma::vec(), arma::mat(), false};
      }
    }
  }

  const arma::mat scores = scaled_data * right;
  arma::vec small_values;
  arma::mat small_vectors;
  if (!arma::eig_sym(
        small_values,
        small_vectors,
        arma::symmatu(scores.t() * scores)
      )) {
    return {arma::vec(), arma::mat(), false};
  }
  const arma::uvec descending = arma::sort_index(small_values, "descend");
  const arma::uvec selected = descending.head(
    static_cast<arma::uword>(retained)
  );
  arma::vec values = small_values.elem(selected);
  values.transform([](double value) { return std::max(value, 0.0); });
  arma::mat vectors = right * small_vectors.cols(selected);
  return {values, vectors, values.is_finite() && vectors.is_finite()};
}

int first_fraction_rank(
    const arma::vec& values,
    const double total,
    const double fraction
) {
  double cumulative = 0.0;
  for (arma::uword index = 0; index < values.n_elem; ++index) {
    cumulative += std::max(values[index], 0.0);
    if (cumulative / total >= fraction) {
      return static_cast<int>(index + 1);
    }
  }
  return static_cast<int>(values.n_elem);
}

}  // namespace

// [[Rcpp::export]]
List cpp_projected_bootstrap_rank(
    const arma::mat& X_target,
    const bool center,
    const int n_boot,
    const double anchor_fraction,
    const double guard_fraction,
    const double target_fraction,
    const double basis_fraction,
    const int initial_rank,
    const int rank_block,
    const int oversample,
    const int bootstrap_extra,
    const int power_iterations,
    const int n_threads,
    const double seed
) {
  const int n = static_cast<int>(X_target.n_rows);
  const int p = static_cast<int>(X_target.n_cols);
  if (n < 2 || p < 2 || !X_target.is_finite()) {
    stop("X_target must have at least 2 rows, 2 columns, and only finite values.");
  }
  if (n_boot < 2) {
    stop("n_boot must be at least 2.");
  }

  arma::mat centered_data = X_target;
  if (center) {
    centered_data.each_row() -= arma::mean(centered_data, 0);
  }
  const double total_variance = arma::accu(arma::square(centered_data)) / n;
  if (!R_FINITE(total_variance) || total_variance <= 0.0) {
    stop("X_target has zero or nonfinite total variance.");
  }
  const int algebraic_rank = std::min(
    p - 1,
    center ? n - 1 : n
  );
  if (algebraic_rank < 1) {
    stop("X_target has no estimable nontrivial eigenspace.");
  }

  const arma::mat scaled_data = centered_data / std::sqrt(
    static_cast<double>(n)
  );
  int requested_rank = std::min(
    algebraic_rank,
    std::max(1, initial_rank)
  );
  ApproximateEigenspace reference;
  const std::uint64_t base_seed = static_cast<std::uint64_t>(seed);
  int reference_attempt = 0;
  while (true) {
    std::mt19937_64 generator(splitmix64(
      base_seed + static_cast<std::uint64_t>(reference_attempt)
    ));
    reference = randomized_eigenspace(
      scaled_data,
      requested_rank,
      oversample,
      power_iterations,
      generator
    );
    if (!reference.success) {
      stop("The randomized target eigenspace calculation failed.");
    }
    const double captured = arma::accu(reference.values) / total_variance;
    if (captured >= basis_fraction || requested_rank >= algebraic_rank) {
      break;
    }
    requested_rank = std::min(
      algebraic_rank,
      std::max(
        requested_rank + std::max(1, rank_block),
        static_cast<int>(std::ceil(1.5 * requested_rank))
      )
    );
    ++reference_attempt;
  }

  const double core_fraction = std::max(
    0.0,
    anchor_fraction - guard_fraction
  );
  const int core_rank = core_fraction > 0.0 ?
    first_fraction_rank(reference.values, total_variance, core_fraction) : 0;
  const int anchor_rank = first_fraction_rank(
    reference.values,
    total_variance,
    anchor_fraction
  );
  const int reference_rank = static_cast<int>(reference.values.n_elem);
  if (core_rank >= reference_rank) {
    stop("The flexible bootstrap band is empty; increase basis_fraction.");
  }
  const int flexible_dimension = reference_rank - core_rank;
  const int target_basis_rank = first_fraction_rank(
    reference.values,
    total_variance,
    std::min(target_fraction, basis_fraction)
  );
  const int flexible_keep = std::min(
    flexible_dimension,
    std::max(1, target_basis_rank - core_rank + bootstrap_extra)
  );

  arma::mat core_vectors;
  arma::mat core_scores;
  if (core_rank > 0) {
    core_vectors = reference.vectors.cols(0, core_rank - 1);
    core_scores = centered_data * core_vectors;
  } else {
    core_scores.set_size(n, 0);
  }
  const arma::mat flexible_vectors = reference.vectors.cols(
    core_rank,
    reference_rank - 1
  );
  const arma::mat flexible_scores = centered_data * flexible_vectors;
  const arma::vec row_norm_squared = arma::sum(
    arma::square(centered_data),
    1
  );

  NumericMatrix bootstrap_fraction(n_boot, flexible_keep);
  NumericMatrix bootstrap_loss(n_boot, flexible_keep);
  std::fill(bootstrap_fraction.begin(), bootstrap_fraction.end(), NA_REAL);
  std::fill(bootstrap_loss.begin(), bootstrap_loss.end(), NA_REAL);
  std::atomic<int> failures(0);
  const int threads = rank_threads(n_threads);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
  for (int bootstrap_index = 0; bootstrap_index < n_boot; ++bootstrap_index) {
    std::mt19937_64 generator(splitmix64(
      base_seed + 1000003ULL + static_cast<std::uint64_t>(bootstrap_index)
    ));
    std::uniform_int_distribution<int> draw(0, n - 1);
    arma::vec counts(n, arma::fill::zeros);
    for (int draw_index = 0; draw_index < n; ++draw_index) {
      counts[draw(generator)] += 1.0;
    }
    const arma::rowvec bootstrap_mean = center ?
      (counts.t() * centered_data) / n :
      arma::rowvec(p, arma::fill::zeros);
    const arma::uvec out_of_bag = arma::find(counts == 0.0);
    const int out_of_bag_n = static_cast<int>(out_of_bag.n_elem);
    if (out_of_bag_n < 2) {
      failures.fetch_add(1, std::memory_order_relaxed);
      continue;
    }
    const arma::rowvec out_of_bag_mean = center ?
      arma::mean(centered_data.rows(out_of_bag), 0) :
      arma::rowvec(p, arma::fill::zeros);
    const double out_of_bag_total =
      arma::mean(row_norm_squared.elem(out_of_bag)) -
      (center ? 2.0 * arma::dot(bootstrap_mean, out_of_bag_mean) -
        arma::dot(bootstrap_mean, bootstrap_mean) : 0.0);
    if (!R_FINITE(out_of_bag_total) || out_of_bag_total <= 0.0) {
      failures.fetch_add(1, std::memory_order_relaxed);
      continue;
    }

    double out_of_bag_core_trace = 0.0;
    if (core_rank > 0) {
      const arma::rowvec core_mean = center ?
        (counts.t() * core_scores) / n :
        arma::rowvec(core_rank, arma::fill::zeros);
      arma::mat out_of_bag_core = core_scores.rows(out_of_bag);
      if (center) {
        out_of_bag_core.each_row() -= core_mean;
      }
      out_of_bag_core_trace = arma::accu(
        arma::square(out_of_bag_core)
      ) / out_of_bag_n;
    }

    arma::mat weighted_flexible = flexible_scores;
    arma::rowvec flexible_mean(
      flexible_dimension,
      arma::fill::zeros
    );
    if (center) {
      flexible_mean = (counts.t() * flexible_scores) / n;
      weighted_flexible.each_row() -= flexible_mean;
    }
    weighted_flexible.each_col() %= arma::sqrt(counts / n);
    ApproximateEigenspace bootstrap_fit = randomized_eigenspace(
      weighted_flexible,
      flexible_keep,
      oversample,
      power_iterations,
      generator
    );
    if (!bootstrap_fit.success) {
      failures.fetch_add(1, std::memory_order_relaxed);
      continue;
    }

    arma::mat out_of_bag_flexible = flexible_scores.rows(out_of_bag);
    if (center) {
      out_of_bag_flexible.each_row() -= flexible_mean;
    }
    const arma::mat out_of_bag_components =
      out_of_bag_flexible * bootstrap_fit.vectors;
    const arma::rowvec component_variance = arma::mean(
      arma::square(out_of_bag_components),
      0
    );
    double cumulative = out_of_bag_core_trace;
    for (int increment = 0; increment < flexible_keep; ++increment) {
      cumulative += component_variance[increment];
      bootstrap_fraction(bootstrap_index, increment) =
        cumulative / out_of_bag_total;
      const int added_rank = increment + 1;
      const arma::mat overlap = bootstrap_fit.vectors.submat(
        0,
        0,
        added_rank - 1,
        added_rank - 1
      );
      const double overlap_trace = arma::accu(arma::square(overlap));
      bootstrap_loss(bootstrap_index, increment) = std::max(
        1.0 - overlap_trace / added_rank,
        0.0
      );
    }
  }
  if (failures.load(std::memory_order_relaxed) > 0) {
    stop("One or more projected bootstrap calculations failed.");
  }

  NumericVector eigenvalues(reference.values.n_elem);
  NumericMatrix reference_vectors(p, reference.values.n_elem);
  std::copy(
    reference.values.begin(), reference.values.end(), eigenvalues.begin()
  );
  std::copy(
    reference.vectors.begin(),
    reference.vectors.end(),
    reference_vectors.begin()
  );
  return List::create(
    _["eigenvalues"] = eigenvalues,
    _["vectors"] = reference_vectors,
    _["total_variance"] = total_variance,
    _["bootstrap_cumulative_fraction"] = bootstrap_fraction,
    _["bootstrap_flexible_subspace_loss"] = bootstrap_loss,
    _["core_rank"] = core_rank,
    _["anchor_rank"] = anchor_rank,
    _["reference_rank"] = reference_rank,
    _["flexible_keep"] = flexible_keep,
    _["basis_fraction_captured"] = arma::accu(reference.values) / total_variance,
    _["n"] = n,
    _["p"] = p,
    _["n_boot"] = n_boot,
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
