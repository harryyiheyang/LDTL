#' POET with MCP thresholding and FSPD positive-definite modification
#'
#' Estimate a factor component from the leading eigenpairs of `S`, threshold the
#' idiosyncratic covariance on correlation scale, and apply FSPD linear
#' shrinkage to the sparse idiosyncratic estimator. The selected spike
#' eigenvalues receive the S-POET positive-part bias correction.
#'
#' @param S Optional symmetric covariance or correlation matrix.
#' @param X Optional individual-level data matrix. Used to compute `S` when
#'   `S` is not supplied and to infer `n`.
#' @param scale If `TRUE`, center and standardize `X` before computing `S`.
#'   Default is `FALSE`, assuming `X` has already been scaled.
#' @param n Sample size used to form `S`. Required when `S` is supplied and
#'   `lambda` is not supplied.
#' @param lambda Optional scalar threshold. Overrides the theoretical default.
#'   If missing, `2 * max(sqrt(log(p) / n), 1 / sqrt(p))` is used.
#' @param eig_min Minimum eigenvalue target for FSPD. Default is 0.
#' @param factor_method Method used to select the number of POET factors.
#'   `"ACT"` uses adjusted correlation thresholding and is the default;
#'   `"D.ratio"` uses the eigenvalue difference-ratio rule. The lower bound is
#'   5 and the upper bound is the rank reaching 90 percent cumulative
#'   eigenvalue mass.
#' @param factors Optional fixed number of POET factors. When supplied, it
#'   overrides `factor_method`.
#' @param eig Optional full eigendecomposition of `S`, supplied as a list with
#'   `values` and `vectors`. If `NULL`, it is computed internally.
#'
#' @return A positive-semidefinite POET regularized correlation matrix. Set
#'   `eig_min` to a positive value to require strict positive definiteness.
#' @export
poet_thresholding <- function(
    S = NULL,
    n = NULL,
    lambda = NULL,
    eig_min = 0,
    X = NULL,
    scale = FALSE,
    factor_method = c("ACT", "D.ratio"),
    eig = NULL,
    factors = NULL
) {
  factor_method <- match.arg(factor_method)
  input <- .ld_resolve_input(S = S, X = X, n = n, name = "S", scale = scale)
  S <- input$S
  n <- input$n
  n <- .ld_validate_n(n)
  comp <- .ld_poet_components(
    S,
    n,
    cutoff_method = factor_method,
    k_min = 5,
    k_max = NULL,
    eig = eig,
    factors = factors
  )
  resolved_lambda <- .ld_resolve_poet_threshold(
    nrow(comp$E),
    n,
    lambda
  )

  E_reg <- thresholding(
    comp$E,
    lambda = resolved_lambda,
    eig_min = eig_min
  )

  Shat <- .ld_symmetrize(comp$P + E_reg)
  stats::cov2cor(Shat)
}
