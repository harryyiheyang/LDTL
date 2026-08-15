#' POET with tapering and FSPD positive-definite modification
#'
#' Estimate a factor component from the leading eigenpairs of `S`, taper the
#' idiosyncratic covariance, and apply FSPD linear shrinkage to the sparse
#' idiosyncratic estimator.
#'
#' @param S Optional symmetric covariance or correlation matrix.
#' @param X Optional individual-level data matrix. Used to compute `S` when
#'   `S` is not supplied and to infer `n`.
#' @param scale If `TRUE`, center and standardize `X` before computing `S`.
#'   Default is `FALSE`, assuming `X` has already been scaled.
#' @param n Sample size used to form `S`. Required when `S` is supplied and
#'   `K` is not supplied.
#' @param K Optional bandwidth. Overrides the theoretical default.
#' @param alpha Bandable covariance smoothness parameter used in the theoretical
#'   bandwidth. Default is 1.
#' @param eig_min Minimum eigenvalue target for FSPD. Default is 0.
#'
#' @return A positive-semidefinite POET regularized correlation matrix. Set
#'   `eig_min` to a positive value to require strict positive definiteness.
#' @export
poet_tapering <- function(
    S = NULL,
    n = NULL,
    K = NULL,
    alpha = 1,
    eig_min = 0,
    X = NULL,
    scale = FALSE
) {
  input <- .ld_resolve_input(S = S, X = X, n = n, name = "S", scale = scale)
  S <- input$S
  n <- input$n
  n <- .ld_validate_n(n)
  comp <- .ld_poet_components(
    S,
    n,
    cutoff_method = "D.ratio",
    k_min = 5,
    k_max = min(15, floor(nrow(S) / 2))
  )

  E_reg <- tapering(
    comp$E,
    n = n,
    K = K,
    alpha = alpha,
    eig_min = eig_min
  )

  Shat <- .ld_symmetrize(comp$P + E_reg)
  stats::cov2cor(Shat)
}
