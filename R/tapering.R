#' Tapering regularization with FSPD positive-definite modification
#'
#' Keep entries with `|i - j| <= K`, linearly taper entries until
#' `|i - j| <= 2 * K`, and set farther entries to zero. If `K` is not
#' supplied, the default bandwidth is capped at `floor(p / 2)`. Positive
#' definiteness is enforced afterward by fixed-support linear shrinkage.
#'
#' @param A Symmetric covariance or correlation matrix.
#' @param X Optional individual-level data matrix. Used to compute `A` when
#'   `A` is not supplied and to infer `n` for the default bandwidth.
#' @param scale If `TRUE`, center and standardize `X` before computing `A`.
#'   Default is `FALSE`, assuming `X` has already been scaled.
#' @param n Sample size used to form `A`. Required when `K` is not supplied.
#' @param K Optional bandwidth. If supplied, it overrides the theoretical
#'   default.
#' @param alpha Bandable covariance smoothness parameter used in the theoretical
#'   bandwidth. Default is 1.
#' @param eig_min Minimum eigenvalue target for FSPD. Default is 0.
#'
#' @return A positive-semidefinite tapered covariance matrix. Set `eig_min` to a
#'   positive value to require strict positive definiteness.
#' @export
tapering <- function(
    A = NULL,
    n = NULL,
    K = NULL,
    alpha = 1,
    eig_min = 0,
    X = NULL,
    scale = FALSE
) {
  input <- .ld_resolve_input(S = A, X = X, n = n, name = "A", scale = scale)
  A <- input$S
  n <- input$n
  p <- nrow(A)
  A <- .ld_symmetrize(A)

  K <- .ld_resolve_bandwidth(p, n, K, alpha)
  W <- .ld_taper_weight(p, K, 2)
  Areg <- .ld_symmetrize(A * W)

  .ld_fspd(Areg, eig_min = eig_min)
}
