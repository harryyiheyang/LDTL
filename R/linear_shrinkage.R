#' Linear shrinkage for a correlation matrix
#'
#' Shrink a sample correlation matrix toward a target matrix. The default
#' shrinkage intensity is intentionally small and tuning-free.
#'
#' @param S Optional symmetric correlation matrix.
#' @param X Optional individual-level data matrix. Used to compute `S` when
#'   `S` is not supplied.
#' @param scale If `TRUE`, center and standardize `X` before computing `S`.
#'   Default is `FALSE`, assuming `X` has already been scaled.
#' @param alpha Shrinkage intensity in `[0, 1]`. Default is 0.05.
#' @param target Optional shrinkage target. Defaults to the identity matrix.
#' @param lambda Deprecated alias for `alpha`, kept for backward
#'   compatibility.
#'
#' @return A positive-semidefinite regularized correlation matrix after
#'   shrinkage.
#' @export
linear_shrinkage <- function(
    S = NULL,
    X = NULL,
    alpha = 0.05,
    target = NULL,
    lambda = NULL,
    scale = FALSE
) {
  input <- .ld_resolve_input(S = S, X = X, name = "S", scale = scale)
  S <- input$S
  p <- ncol(S)
  if (is.null(target)) {
    target <- diag(1, p)
  } else {
    target <- .ld_as_square_matrix(target, "target")
    if (!all(dim(target) == c(p, p))) {
      stop("target must have the same dimensions as S.", call. = FALSE)
    }
  }

  S <- .ld_symmetrize(S)
  target <- .ld_symmetrize(target)

  if (!is.null(lambda)) {
    alpha <- lambda
  }
  if (length(alpha) != 1L || !is.finite(alpha)) {
    stop("alpha must be a finite scalar.", call. = FALSE)
  }
  alpha <- min(max(alpha[1L], 0), 1)

  out <- .ld_symmetrize(alpha * target + (1 - alpha) * S)
  out <- .ld_fspd(out, eig_min = 0)
  out <- stats::cov2cor(out)
  out[is.na(out)] <- 0
  out <- .ld_symmetrize(out)
  diag(out) <- 1
  out
}
