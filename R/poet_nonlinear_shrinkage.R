#' POET with nonlinear shrinkage of the idiosyncratic component
#'
#' Estimate a factor component from the leading eigenpairs of the target
#' correlation matrix, then regularize the idiosyncratic component in the
#' orthogonal complement with nonlinear shrinkage. The nonlinear estimator is
#' blended with the sample residual covariance, controlled by `shrinkage`.
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
#'
#' @return A positive-semidefinite POET regularized correlation matrix.
#' @export
poet_nonlinear_shrinkage <- function(
    X,
    S = NULL,
    n = NULL,
    shrinkage = 0,
    scale = FALSE
) {
  if (missing(X) || is.null(X)) {
    stop("poet_nonlinear_shrinkage requires individual-level X.", call. = FALSE)
  }
  input <- .ld_resolve_input(S = S, X = X, n = n, name = "S", scale = scale)
  n <- .ld_validate_n(input$n)
  comp <- .ld_poet_components(
    input$S,
    n,
    cutoff_method = "D.ratio",
    k_min = 5,
    k_max = min(15, floor(nrow(input$S) / 2))
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
