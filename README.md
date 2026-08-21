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

- cov_tl: Tuning-free covariance transfer learning. It estimates the optimal
  source pooling rate from the closed-form covariance URE criterion, without
  cross-validation or a tuning grid.

- source_moments_method: C++/OpenMP individual-block preprocessing of a source
  cohort. It returns a reusable source covariance and scalar fourth-moment
  noise estimate without retaining the source individuals.

- cov_tl_stabilized: Finite-source stabilized covariance transfer learning.
  This is a separate estimator from `cov_tl`; it uses source fourth-moment
  information when an `ld_source_moments` object is supplied.

- eigspac_tl: Tuning-free eigenspace transfer learning. It estimates the optimal
  first-order pooling rate in spectral-projector tangent geometry and returns
  the pooled leading eigenspace.

- eigspac_tl_stabilized: A separate finite-source eigenspace rule. Because a
  scalar source fourth moment does not identify arbitrary target-direction
  tangent noise, it either uses a supplied target-specific tangent variance or
  an explicitly labelled Frobenius/eigengap upper proxy.

Most matrix estimators accept either a precomputed matrix (`S` or `A`) or
individual-level data `X`. Nonlinear shrinkage requires `X`. The POET functions
automatically select the number of latent factors via the ratio-type criterion.
Sparse methods use theory-rate defaults only when the user does not supply the
tuning parameter. Standalone thresholding uses
`lambda = 2 * sqrt(log(p) / n)`, while POET thresholding uses
`lambda = 2 * max(sqrt(log(p) / n), 1 / sqrt(p))`. Banding and tapering use
`K = ceiling(n^(1 / (2 * alpha + 1)))` with `alpha = 1` by default, capped at
`floor(p / 2)`. These defaults require `n`; otherwise pass `lambda` or `K`
directly. Standalone linear shrinkage uses `alpha = 0.05` by default; POET
linear shrinkage uses `alpha = 0.5` by default. Nonlinear shrinkage has
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

eig_fit <- eigspac_tl(X_target, S_source, rank = 3)
U_tl <- eig_fit$vectors
P_tl <- eig_fit$projector
```

For `eigspac_tl()`, an independently estimated `S_pilot` can be supplied to make
the tangent-risk criterion conditionally unbiased for a fixed spectral
derivative. Without it, the function uses the full-sample first-order plug-in
rule. If the mean must be estimated from `X_target`, use `center = TRUE`; the
resulting weight is then a practical plug-in rule. Neither function performs
CV or grid search.

When source individual-level data are available, preprocess them once and
reuse the resulting summary across target fits:

```r
source_fit <- source_moments_method(
  X_source,
  center = FALSE,
  block_size = NULL,  # automatic 512 MiB working-set target
  n_threads = 20
)

cov_stable <- cov_tl_stabilized(
  X_target,
  source_fit,
  center = FALSE
)

eig_stable <- eigspac_tl_stabilized(
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

The legacy and finite-source estimators deliberately have different entry
points. Calling `cov_tl()` or `eigspac_tl()` never silently activates a
stabilized weight.

Teacher B's path-adaptive framework is exposed through two more independent
entry points. Both interpret each grid value as an effective fraction of the
available source sample, adjust its weight separately inside every target
fold, and refit the selected fraction on all target observations:

```r
fold_id <- sample(rep(1:5, length.out = nrow(X_target)))

path_max <- path_tpca_max_score(
  X_target,
  source_fit,       # or S_source together with n_source
  rank = 3,
  fold_id = fold_id
)

path_one_se <- path_tpca_one_se(
  X_target,
  source_fit,
  rank = 3,
  fold_id = fold_id
)
```

`path_tpca_max_score()` selects the raw held-out target-score maximizer.
`path_tpca_one_se()` instead selects the largest source fraction whose paired
per-individual score deficit is within one standard error of that maximizer.
These two path estimators use the source covariance and sample size, but not
the fourth moment. Thus the package now keeps six public estimators separate:
legacy CovTL and EigenTL, their finite-source fourth-moment variants, and the
two Teacher-B covariance-path refits.

Choose a common target subspace dimension before comparing the eigenspace and
path estimators. `select_tl_rank_bootstrap()` treats cumulative variance as an
operational compression goal rather than a finite-spike rank. It fixes a core
below a 50-percent anchor, leaves a guard band flexible, and uses individual
bootstrap samples only in a randomized projected score space:

```r
rank_fit <- select_tl_rank_bootstrap(
  X_target,
  target_fraction_grid = c(.5, .6, .7, .8, .9),
  selection_fraction = .8,
  anchor_fraction = .5,
  guard_fraction = .05,
  n_boot = 100,
  n_threads = 20
)

K <- rank_fit$rank
rank_fit$rank_path

eig_fit <- eigspac_tl(X_target, S_source, rank = K)
path_fit <- path_tpca_one_se(X_target, source_fit, rank = K)
```

The 0.8 target is user-defined; the out-of-bag bootstrap lower quantile is the
data-driven quantity. Flexible directions are fitted on in-bag individuals and
their captured variance is evaluated out of bag. Stability is audited at the
smallest rank reaching each target but does not push the rank upward, because
normalized projector loss can decrease mechanically near the full flexible
space. No source data, bootstrap `p` by `p` covariance, or repeated full
eigendecomposition is used by this selector.

Both transfer-learning functions prioritize `CppMatrix` for covariance
construction, matrix products, centering, projector construction, and spectral
decomposition. Base R is used only as a guarded fallback or for scalar
reductions not exposed by `CppMatrix`. In `eigspac_tl()`, the tangent variance is
computed with a CppMatrix cross-product identity rather than an
`rank * (p - rank)` R loop.

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
