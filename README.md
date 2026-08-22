# LDRegularization

The `LDRegularization` package provides a collection of covariance and correlation matrix regularization methods, with a focus on high-dimensional settings such as linkage disequilibrium (LD) matrices in statistical genetics. The package implements sparse regularization, shrinkage, nonlinear shrinkage, and POET-style estimators.

These tools help ensure positive semi-definiteness and improve estimation accuracy when working with large, noisy covariance or correlation matrices.

## Installation

You can install the `LDRegularization` package from GitHub using the `devtools` package:

```r
devtools::install_github("harryyiheyang/LDRegularization")
```

## Functions

The package currently includes the following main functions:

- poet_thresholding: POET estimator with entrywise thresholding on the idiosyncratic covariance.

- poet_banding: POET estimator with banding (keep entries within a fixed bandwidth).

- poet_tapering: POET estimator with tapering (apply a tapering kernel depending on |i-j|).

- poet_linear_shrinkage: POET estimator with linear shrinkage of the idiosyncratic covariance.

- poet_nonlinear_shrinkage: POET estimator with mixed nonlinear shrinkage of
  the idiosyncratic covariance.

- thresholding: Standalone entrywise MCP thresholding for covariance/correlation matrices.

- banding: Standalone banding with bandwidth `K`.

- tapering: Standalone tapering with bandwidth `K`.

- linear_shrinkage: Standalone linear shrinkage.

- nonlinear_shrinkage: Standalone mixed nonlinear shrinkage for individual-level
  data.

- cov_tl: Tuning-free covariance transfer learning. A source covariance matrix
  uses the summary-only URE rule; an `ld_source_moments` object automatically
  includes finite-source noise.

- source_moments_method: C++/OpenMP individual-block preprocessing of a source
  cohort. It returns a reusable source covariance and scalar fourth-moment
  noise estimate without retaining the source individuals.

- eigen_tl: Fold-adjusted covariance-path eigenspace transfer learning. Its
  default paired one-standard-error selector favors more transfer among
  competitive candidates; `method = "min"` selects minimum reconstruction
  risk.

- multi_source_tl: Almost tuning-free aggregation of source-specific CovTL or
  EigenTL fits using shrinkage gain or held-out reconstruction gain.

Most matrix estimators accept either a precomputed matrix (`S` or `A`) or
individual-level data `X`. Nonlinear shrinkage requires `X`. The POET functions
use ACT factor selection by default, retain `D.ratio` as an option, and accept an
optional reusable full eigendecomposition through `eig`.
Sparse methods use theory-rate defaults only when the user does not supply the
tuning parameter. Standalone thresholding uses
`lambda = 2 * sqrt(log(p) / n)`, while POET thresholding uses
`lambda = 2 * max(sqrt(log(p) / n), 1 / sqrt(p))`. Banding and tapering use
`K = ceiling(n^(1 / (2 * alpha + 1)))` with `alpha = 1` by default, capped at
`floor(p / 2)`. These defaults require `n`; otherwise pass `lambda` or `K`
directly. Standalone and POET linear shrinkage use a Gaussian/Wishart MSE
plug-in intensity capped at `alpha = 0.05`. Nonlinear shrinkage has
`shrinkage = 0` by default, so scripts should set it explicitly, for example
`shrinkage = 0.5`, when a stronger nonlinear component is desired.

Sparse estimators are made positive semidefinite with fixed-support linear shrinkage, using `eig_min = 0` by default. Set a positive `eig_min` only when downstream computations require a strictly positive-definite matrix. FSPD is a final safeguard and does not choose the sparsity tuning parameter.

## Dependencies

The package makes use of efficient matrix operations implemented in CppMatrix, which relies on Rcpp and RcppArmadillo for performance.

## Example

```r
library(LDRegularization)

# Example: POET with thresholding
set.seed(123)
p <- 50
n <- 200
X <- matrix(rnorm(n * p), n, p)
S <- cov(X)
Sigma_hat <- poet_thresholding(S, n = n)
```

Closed-form transfer learning uses target individual-level data and only a
source covariance summary. For the exact known-mean URE interpretation, the
rows of `X_target` must be independent with a known population mean (zero in
this example):

```r
X_target <- X
S_source <- diag(p)

cov_fit <- cov_tl(X_target, S_source)
cov_fit$lambda       # weight on the source covariance
cov_fit$alpha        # lambda / (1 - lambda)
Sigma_tl <- cov_fit$covariance

eig_fit <- eigen_tl(
  X_target,
  S_source,
  rank = 3,
  n_source = 1000
)
U_tl <- eig_fit$vectors
P_tl <- eig_fit$projector
```

`eigen_tl()` uses `method = "one_se"` by default. Set `method = "min"` to
select the raw held-out reconstruction-risk minimizer. Supply a common
`fold_id` when comparing multiple sources.

When source individual-level data are available, preprocess them once and
reuse the resulting summary across target fits:

```r
source_fit <- source_moments_method(
  X_source,
  center = FALSE,
  block_size = NULL,  # automatic 512 MiB working-set target
  n_threads = 20
)

cov_finite <- cov_tl(
  X_target,
  source_fit,
  center = FALSE
)

eig_finite <- eigen_tl(
  X_target,
  source_fit,
  rank = 3,
  center = FALSE
)
```

`source_moments_method()` processes individuals in native blocks and never
forms individual outer products. If a matching `S_source` is already
available, pass it to the moment function; the additional source scan then
computes only `sum_i ||X_source[i, ]||^4`, reducing the extra work to
`O(n_source * p)`. Missing genotypes should be mean-imputed before the scan.
The covariance output is still `p` by `p`, so variant/LD-region blocking is a
separate requirement when `p` itself is too large.

Teacher B's path-adaptive framework is the package's EigenTL implementation.
Each grid value is an effective fraction of the available source sample. Its
weight is adjusted separately inside every target fold, and the selected
fraction is refitted on all target observations:

```r
fold_id <- sample(rep(1:5, length.out = nrow(X_target)))

eig_min <- eigen_tl(
  X_target,
  source_fit,       # or S_source together with n_source
  rank = 3,
  fold_id = fold_id,
  method = "min"
)

eig_one_se <- eigen_tl(
  X_target,
  source_fit,
  rank = 3,
  fold_id = fold_id
)
```

The path uses the source covariance and sample size, but not the fourth moment.
A common subspace dimension should be chosen externally before comparing
sources.

Multiple source-specific fits can be aggregated without pre-mixing source
covariances:

```r
cov_fits <- list(
  EUR = cov_tl(X_target, source_eur),
  AFR = cov_tl(X_target, source_afr),
  AMR = cov_tl(X_target, source_amr)
)
cov_multi <- multi_source_tl(cov_fits)

eig_fits <- list(
  EUR = eigen_tl(X_target, source_eur, rank = 3, fold_id = fold_id),
  AFR = eigen_tl(X_target, source_afr, rank = 3, fold_id = fold_id),
  AMR = eigen_tl(X_target, source_amr, rank = 3, fold_id = fold_id)
)
eig_multi <- multi_source_tl(eig_fits)
```

CovTL candidates are weighted by their normalized shrinkage coefficients,
which act as plug-in single-source risk-gain scores. EigenTL candidates are
weighted by the positive held-out reconstruction-score gain of their selected
path points over zero transfer. Equal source priors are used by default.

The transfer-learning functions prioritize `CppMatrix` for covariance
construction, matrix products, centering, projector construction, and spectral
decomposition. Base R is used only as a guarded fallback or for scalar
reductions not exposed by `CppMatrix`.

The source-target Frobenius distance is evaluated directly as
`norm(S_source - S_target, "F")^2`. This avoids the cancellation that can
occur in the algebraically equivalent expression
`norm(S_source, "F")^2 + norm(S_target, "F")^2 -
2 * tr(t(S_source) %*% S_target)` when the two large matrices are close.

## License

This package is licensed under the MIT License.

## Contact

Yihe Yang
Email: yxy1234@case.edu
