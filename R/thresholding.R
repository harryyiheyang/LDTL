#' Entrywise MCP thresholding with FSPD positive-definite modification
#'
#' Threshold the off-diagonal entries of the correlation version of a covariance
#' matrix and map the result back to covariance scale. By default the threshold
#' is the high-dimensional rate `2 * sqrt(log(p) / n)`. Positive definiteness
#' is enforced afterward by fixed-support linear shrinkage.
#'
#' @param S Symmetric covariance or correlation matrix.
#' @param X Optional individual-level data matrix. Used to compute `S` when
#'   `S` is not supplied and to infer `n` for the default threshold.
#' @param scale If `TRUE`, center and standardize `X` before computing `S`.
#'   Default is `FALSE`, assuming `X` has already been scaled.
#' @param n Sample size used to form `S`. Required when `lambda` is not supplied.
#' @param lambda Optional scalar or matrix threshold. If supplied, it overrides
#'   the theoretical default.
#' @param eig_min Minimum eigenvalue target for FSPD. Default is 0.
#'
#' @return A positive-semidefinite thresholded covariance matrix. Set `eig_min`
#'   to a positive value to require strict positive definiteness.
#' @export
thresholding <- function(
    S = NULL,
    n = NULL,
    lambda = NULL,
    eig_min = 0,
    X = NULL,
    scale = FALSE
) {
  input <- .ld_resolve_input(S = S, X = X, n = n, name = "S", scale = scale)
  S <- input$S
  n <- input$n
  p <- nrow(S)
  S <- .ld_symmetrize(S)

  diag_vals <- sqrt(pmax(diag(S), 0))
  S_cor <- stats::cov2cor(S)
  S_cor[is.na(S_cor)] <- 0
  diag(S_cor) <- 0

  lambda <- .ld_resolve_threshold(p, n, lambda)
  if (length(lambda) == 1L) {
    cut <- lambda
  } else {
    cut <- .ld_as_square_matrix(lambda, "lambda")
    if (!all(dim(cut) == c(p, p))) {
      stop("lambda matrix must have the same dimensions as S.", call. = FALSE)
    }
    diag(cut) <- 0
  }

  Mcor <- .ld_mcp(S_cor, lam = cut, a = 3)
  Mcor <- .ld_symmetrize(Mcor)
  diag(Mcor) <- 1

  Sigma <- t(t(Mcor) * diag_vals) * diag_vals
  Sigma <- .ld_symmetrize(Sigma)
  diag(Sigma) <- diag_vals^2

  .ld_fspd(Sigma, eig_min = eig_min)
}
