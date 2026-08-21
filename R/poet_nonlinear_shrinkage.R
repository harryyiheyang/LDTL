#' POET with nonlinear shrinkage of the idiosyncratic component
#'
#' Estimate a factor component from the leading eigenpairs of the target
#' correlation matrix, then regularize the idiosyncratic component in the
#' orthogonal complement with nonlinear shrinkage. The nonlinear estimator is
#' blended with the sample residual covariance, controlled by `shrinkage`. The
#' selected spike eigenvalues receive the S-POET positive-part bias correction.
#'
#' @param X Individual-level data matrix with observations in rows.
#' @param S Optional symmetric covariance or correlation matrix. If omitted,
#'   it is computed from `X`.
#' @param n Optional sample size. Inferred from `X` when omitted.
#' @param shrinkage Mixing weight for the nonlinear residual estimator.
#'   Default is 0, which keeps the sample residual and applies the
#'   positive-semidefinite safeguard.
#' @param scale If `TRUE`, center and standardize `X` before estimation.
#'   Default is `FALSE`, assuming `X` has already been scaled.
#' @param factor_method Method used to select the number of POET factors.
#'   `"ACT"` uses adjusted correlation thresholding and is the default;
#'   `"D.ratio"` uses the eigenvalue difference-ratio rule. The lower bound is
#'   5 and the upper bound is the rank reaching 90 percent cumulative
#'   eigenvalue mass.
#' @param eig Optional full eigendecomposition of `S`, supplied as a list with
#'   `values` and `vectors`. If `NULL`, it is computed internally.
#'
#' @return A positive-semidefinite POET regularized correlation matrix.
#' @export
poet_nonlinear_shrinkage <- function(
    X,
    S = NULL,
    n = NULL,
    shrinkage = 0,
    scale = FALSE,
    factor_method = c("ACT", "D.ratio"),
    eig = NULL
) {
  factor_method <- match.arg(factor_method)
  if (missing(X) || is.null(X)) {
    stop("poet_nonlinear_shrinkage requires individual-level X.", call. = FALSE)
  }
  input <- .ld_resolve_input(S = S, X = X, n = n, name = "S", scale = scale)
  n <- .ld_validate_n(input$n)
  comp <- .ld_poet_components(
    input$S,
    n,
    cutoff_method = factor_method,
    k_min = 5,
    k_max = NULL,
    eig = eig
  )

  E_reg <- .ld_nonlinear_residual(
    input$X,
    comp$U,
    shrinkage = shrinkage,
    scale = FALSE
  )
  out <- .ld_fspd(.ld_symmetrize(comp$P + E_reg), eig_min = 0)
  out <- stats::cov2cor(out)
  out[is.na(out)] <- 0
  out <- .ld_symmetrize(out)
  diag(out) <- 1
  out
}
