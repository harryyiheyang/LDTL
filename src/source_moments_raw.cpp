#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

struct TokenBounds {
  std::size_t begin;
  std::size_t end;
};

int raw_resolve_threads(const int requested) {
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

bool next_token(
    const std::string& line,
    std::size_t& position,
    TokenBounds& token
) {
  const std::size_t size = line.size();
  while (position < size &&
         std::isspace(static_cast<unsigned char>(line[position]))) {
    ++position;
  }
  if (position >= size) {
    return false;
  }
  token.begin = position;
  while (position < size &&
         !std::isspace(static_cast<unsigned char>(line[position]))) {
    ++position;
  }
  token.end = position;
  return true;
}

std::string token_string(
    const std::string& line,
    const TokenBounds& token
) {
  return line.substr(token.begin, token.end - token.begin);
}

bool token_missing(
    const std::string& line,
    const TokenBounds& token
) {
  const std::size_t length = token.end - token.begin;
  if (length == 1 && line[token.begin] == '.') {
    return true;
  }
  return length == 2 &&
    line[token.begin] == 'N' && line[token.begin + 1] == 'A';
}

double token_double(
    const std::string& line,
    const TokenBounds& token,
    const std::size_t line_number
) {
  const char* begin = line.c_str() + token.begin;
  char* parsed_end = NULL;
  const double value = std::strtod(begin, &parsed_end);
  if (parsed_end != line.c_str() + token.end || !R_FINITE(value)) {
    stop(
      "Invalid dosage token in PLINK raw file at line " +
      std::to_string(line_number) + "."
    );
  }
  return value;
}

CharacterVector raw_header_variants(
    const std::string& header,
    const std::string& raw_file
) {
  std::size_t position = 0;
  TokenBounds token;
  std::vector<std::string> fields;
  while (next_token(header, position, token)) {
    fields.push_back(token_string(header, token));
  }
  if (fields.size() < 7) {
    stop("PLINK raw file must contain six metadata columns and at least one variant: " + raw_file);
  }
  const char* expected[6] = {
    "FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"
  };
  for (std::size_t j = 0; j < 6; ++j) {
    if (fields[j] != expected[j]) {
      stop("PLINK raw header has unexpected metadata columns: " + raw_file);
    }
  }

  CharacterVector variants(fields.size() - 6);
  for (std::size_t j = 6; j < fields.size(); ++j) {
    const std::string& label = fields[j];
    const std::size_t delimiter = label.rfind('_');
    if (delimiter == std::string::npos || delimiter == 0) {
      stop("PLINK raw variant column lacks a counted-allele suffix: " + label);
    }
    variants[j - 6] = label.substr(0, delimiter);
  }
  return variants;
}

int read_raw_block(
    std::ifstream& input,
    const int p,
    const int block_size,
    std::vector<double>& block,
    std::size_t& line_number,
    long double& missing_calls
) {
  const double missing = std::numeric_limits<double>::quiet_NaN();
  std::string line;
  int rows = 0;

  while (rows < block_size && std::getline(input, line)) {
    ++line_number;
    if (line.empty()) {
      continue;
    }
    std::size_t position = 0;
    TokenBounds token;
    for (int j = 0; j < 6; ++j) {
      if (!next_token(line, position, token)) {
        stop(
          "PLINK raw row has fewer than six metadata fields at line " +
          std::to_string(line_number) + "."
        );
      }
    }
    for (int j = 0; j < p; ++j) {
      if (!next_token(line, position, token)) {
        stop(
          "PLINK raw row has too few dosage fields at line " +
          std::to_string(line_number) + "."
        );
      }
      if (token_missing(line, token)) {
        block[rows + block_size * j] = missing;
        missing_calls += 1.0L;
      } else {
        block[rows + block_size * j] =
          token_double(line, token, line_number);
      }
    }
    if (next_token(line, position, token)) {
      stop(
        "PLINK raw row has too many dosage fields at line " +
        std::to_string(line_number) + "."
      );
    }
    ++rows;
  }
  return rows;
}

void open_raw(
    const std::string& raw_file,
    std::ifstream& input,
    std::string& header
) {
  input.open(raw_file.c_str(), std::ios::in | std::ios::binary);
  if (!input.is_open()) {
    stop("Cannot open PLINK raw file: " + raw_file);
  }
  if (!std::getline(input, header)) {
    stop("PLINK raw file is empty: " + raw_file);
  }
}

}  // namespace

// [[Rcpp::export]]
List cpp_source_moments_raw(
    const std::string raw_file,
    const int block_size,
    const int n_threads
) {
  if (block_size < 1) {
    stop("block_size must be a positive integer.");
  }
  const int threads = raw_resolve_threads(n_threads);

  std::ifstream first_pass;
  std::string header;
  open_raw(raw_file, first_pass, header);
  CharacterVector variants = raw_header_variants(header, raw_file);
  const int p = variants.size();
  std::vector<double> block(
    static_cast<std::size_t>(block_size) * static_cast<std::size_t>(p)
  );
  std::vector<long double> sums(p, 0.0L);
  std::vector<long double> sums_squared(p, 0.0L);
  std::vector<long double> observed(p, 0.0L);
  long double missing_first = 0.0L;
  std::size_t line_number = 1;
  long double n_long = 0.0L;

  while (true) {
    const int rows = read_raw_block(
      first_pass, p, block_size, block, line_number, missing_first
    );
    if (rows == 0) {
      break;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (int j = 0; j < p; ++j) {
      long double block_sum = 0.0L;
      long double block_sum_squared = 0.0L;
      long double block_observed = 0.0L;
      const std::size_t offset = static_cast<std::size_t>(block_size) * j;
      for (int i = 0; i < rows; ++i) {
        const double value = block[offset + i];
        if (R_FINITE(value)) {
          const long double extended = static_cast<long double>(value);
          block_sum += extended;
          block_sum_squared += extended * extended;
          block_observed += 1.0L;
        }
      }
      sums[j] += block_sum;
      sums_squared[j] += block_sum_squared;
      observed[j] += block_observed;
    }
    n_long += static_cast<long double>(rows);
    checkUserInterrupt();
  }
  first_pass.close();

  if (n_long < 2.0L) {
    stop("PLINK raw file must contain at least two samples.");
  }

  NumericVector means(p);
  NumericVector scales(p);
  NumericVector observed_counts(p);
  for (int j = 0; j < p; ++j) {
    if (observed[j] < 1.0L) {
      stop("A PLINK raw variant has no observed dosage: " + as<std::string>(variants[j]));
    }
    const long double mean = sums[j] / observed[j];
    long double centered_sum_squared =
      sums_squared[j] - sums[j] * sums[j] / observed[j];
    const long double tolerance =
      std::numeric_limits<double>::epsilon() *
      std::max(1.0L, sums_squared[j]);
    if (centered_sum_squared < 0.0L &&
        std::abs(centered_sum_squared) <= tolerance) {
      centered_sum_squared = 0.0L;
    }
    if (centered_sum_squared <= 0.0L) {
      stop("A PLINK raw variant has nonpositive variance: " + as<std::string>(variants[j]));
    }
    means[j] = static_cast<double>(mean);
    scales[j] = static_cast<double>(
      std::sqrt(centered_sum_squared / n_long)
    );
    observed_counts[j] = static_cast<double>(observed[j]);
  }

  std::ifstream second_pass;
  std::string second_header;
  open_raw(raw_file, second_pass, second_header);
  if (second_header != header) {
    stop("PLINK raw header changed between file passes.");
  }

  NumericMatrix covariance(p, p);
  long double fourth_sum = 0.0L;
  long double missing_second = 0.0L;
  long double n_second = 0.0L;
  line_number = 1;

  while (true) {
    const int rows = read_raw_block(
      second_pass, p, block_size, block, line_number, missing_second
    );
    if (rows == 0) {
      break;
    }
    std::vector<long double> row_norm_squared(rows, 0.0L);
    for (int j = 0; j < p; ++j) {
      const std::size_t offset = static_cast<std::size_t>(block_size) * j;
      const double mean = means[j];
      const double scale = scales[j];
      for (int i = 0; i < rows; ++i) {
        const double value = block[offset + i];
        const double standardized =
          R_FINITE(value) ? (value - mean) / scale : 0.0;
        block[offset + i] = standardized;
        const long double extended =
          static_cast<long double>(standardized);
        row_norm_squared[i] += extended * extended;
      }
    }
    for (int i = 0; i < rows; ++i) {
      fourth_sum += row_norm_squared[i] * row_norm_squared[i];
    }

    const char upper = 'U';
    const char transpose = 'T';
    const int p_blas = p;
    const int rows_blas = rows;
    const int leading_block = block_size;
    const int leading_covariance = p;
    const double alpha = 1.0;
    const double beta = 1.0;
    F77_CALL(dsyrk)(
      &upper,
      &transpose,
      &p_blas,
      &rows_blas,
      &alpha,
      block.data(),
      &leading_block,
      &beta,
      covariance.begin(),
      &leading_covariance FCONE FCONE
    );

    n_second += static_cast<long double>(rows);
    checkUserInterrupt();
  }
  second_pass.close();

  if (n_second != n_long || missing_second != missing_first) {
    stop("PLINK raw file changed between file passes.");
  }

  long double frobenius_squared = 0.0L;
  for (int j = 0; j < p; ++j) {
    for (int k = j; k < p; ++k) {
      const double value = covariance[j + p * k] /
        static_cast<double>(n_long);
      covariance[j + p * k] = value;
      covariance[k + p * j] = value;
      const long double square =
        static_cast<long double>(value) * static_cast<long double>(value);
      frobenius_squared += (j == k) ? square : 2.0L * square;
    }
  }
  const long double variance_raw =
    (fourth_sum - n_long * frobenius_squared) /
    (n_long * (n_long - 1.0L));
  const double working_bytes =
    static_cast<double>(block_size) * static_cast<double>(p) * sizeof(double) +
    static_cast<double>(p) * static_cast<double>(p) * sizeof(double) +
    static_cast<double>(block_size) * sizeof(long double);

  return List::create(
    _["covariance"] = covariance,
    _["variance_raw"] = static_cast<double>(variance_raw),
    _["fourth_sum"] = static_cast<double>(fourth_sum),
    _["n"] = static_cast<double>(n_long),
    _["p"] = p,
    _["mean"] = means,
    _["scale"] = scales,
    _["observed"] = observed_counts,
    _["missing_calls"] = static_cast<double>(missing_first),
    _["variants"] = variants,
    _["block_size"] = block_size,
    _["threads_used"] = threads,
    _["passes"] = 2,
    _["working_bytes"] = working_bytes
  );
}
