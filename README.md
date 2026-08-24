# LDTL

`LDTL` provides single-source and multi-source transfer learning for linkage
disequilibrium covariance matrices and eigenspaces. It also retains the
package's POET and conventional covariance-regularization methods.

## Installation

```r
devtools::install_github("harryyiheyang/LDTL")
library(LDTL)
```

## Core transfer-learning methods

The four primary estimators are:

- `cov_tl()`: tuning-free single-source covariance transfer;
- `eigen_tl()`: single-source reconstruction-selected eigenspace path;
- `multisource_cov_tl()`: tuning-free joint multi-source covariance transfer;
- `multisource_eigen_tl()`: a fold-wise joint source composition followed by
  one reconstruction-selected eigenspace path.

Source covariance and fourth-moment summaries can be computed once and reused
for every target population and fold:

```r
source_eur <- source_moments(X_eur)
source_afr <- source_moments(X_afr)
source_amr <- source_moments(X_amr)

sources <- list(EUR = source_eur, AFR = source_afr, AMR = source_amr)
```

For a large PLINK 2 additive raw file, use the streaming interface:

```r
source_eur <- source_moments_plink("eur.raw")
```

Both source-summary functions accept an optional precomputed `R_source`. When
it is supplied, the source covariance is reused and only the fourth-moment
scan is performed.

### Single-source covariance transfer

```r
fit_cov <- cov_tl(X_target, source_eur)
R_cov <- fit_cov$covariance
fit_cov$lambda
```

### Single-source eigenspace transfer

```r
fit_eigen <- eigen_tl(
  X_target,
  source_eur,
  rank = 0.99,
  fold_id = fold_id,
  method = "one_se",
  eigen_solver = "rspectra"
)
P_eigen <- fit_eigen$projector
```

Set `rank` to a positive integer for a fixed eigenspace dimension, or to a
number strictly between zero and one for a target cumulative explained-variance
threshold. For example, `rank = 0.99` learns K99 once from the full target data
and uses that same integer rank in every target fold and path candidate. When
K99 is already available, pass the integer directly. The opt-in
`eigen_solver = "rspectra"` backend computes only those leading eigenpairs;
the default `"full"` backend preserves the original CppMatrix behavior.

`method = "min"` selects the maximum held-out score. The default
`method = "one_se"` selects the most transferred statistically competitive
candidate.

### Multi-source covariance transfer

```r
fit_multi_cov <- multisource_cov_tl(X_target, sources)
R_multi_cov <- fit_multi_cov$covariance
fit_multi_cov$source_weights
```

The estimator solves one nonnegative joint regression with source coefficients
summing to at most one. Source-source cross terms are retained in the small
Gram matrix; source covariances are not pre-mixed by sample size.

### Multi-source eigenspace transfer

```r
fit_multi_eigen <- multisource_eigen_tl(
  X_target,
  sources,
  rank = 0.99,
  fold_id = fold_id,
  method = "one_se",
  eigen_solver = "rspectra"
)
P_multi_eigen <- fit_multi_eigen$projector
```

Every target training fold relearns its target covariance, target variance,
and convex source composition. The resolved source summaries are reused. A
single direct covariance-weight path is then evaluated by held-out target
reconstruction score, so three sources require one eigenspace path rather than
three independent paths.

Fit objects retain the final learned covariance/eigenspace, weights, risks,
scores, and small optimization diagnostics. They do not retain copies of the
input target or source covariance matrices.

## Additional covariance regularizers

The existing methods remain available:

- `linear_shrinkage()` and `nonlinear_shrinkage()`;
- `banding()`, `tapering()`, and `thresholding()`;
- `poet_linear_shrinkage()` and `poet_nonlinear_shrinkage()`;
- `poet_banding()`, `poet_tapering()`, and `poet_thresholding()`.

Matrix multiplication, eigendecomposition, and covariance construction prefer
`CppMatrix`. Source summaries and multi-source Gram matrices use native
streaming kernels to avoid unnecessary large intermediate matrices.
