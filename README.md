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

- eigspac_tl: Tuning-free eigenspace transfer learning. It estimates the optimal
  first-order pooling rate in spectral-projector tangent geometry and returns
  the pooled leading eigenspace.

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

Sparse estimators are made positive definite with fixed-support positive-definite (FSPD) linear shrinkage, using `eigenmin = 0.001` by default. FSPD is a final safeguard and does not choose the sparsity tuning parameter.

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
