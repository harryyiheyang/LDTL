#' Linear shrinkage for a correlation matrix
#'
#' Shrink a sample correlation matrix toward a target matrix. The shrinkage
#' intensity follows a Gaussian/Wishart MSE plug-in rule capped by `alpha`.
#'
#' @param S Optional symmetric correlation matrix.
#' @param X Optional individual-level data matrix. Used to compute `S` when
#'   `S` is not supplied.
#' @param scale If `TRUE`, center and standardize `X` before computing `S`.
#'   Default is `FALSE`, assuming `X` has already been scaled.
#' @param alpha Boundary in `[0, 1]` for the MSE plug-in shrinkage intensity.
#'   The applied intensity is the smaller of the MSE rule and `alpha`. Default
#'   is 0.05.
#' @param target Optional shrinkage target. Defaults to the identity matrix.
#' @param lambda Deprecated alias for `alpha`, kept for backward
#'   compatibility.
#' @param n Sample size used to form `S`. Inferred from `X`; required when only
#'   `S` is supplied.
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
    scale = FALSE,
    n = NULL
) {
  input <- .ld_resolve_input(
    S = S,
    X = X,
    n = n,
    name = "S",
    scale = scale
  )
  S <- input$S
  n <- .ld_validate_n(input$n)
  p <- ncol(S)
  default_target <- is.null(target)
  if (!default_target) {
    target <- .ld_as_square_matrix(target, "target")
    if (!all(dim(target) == c(p, p))) {
      stop("target must have the same dimensions as S.", call. = FALSE)
    }
  }

  S <- .ld_symmetrize(S)
  if (!default_target) {
    target <- .ld_symmetrize(target)
  }

  if (!is.null(lambda)) {
    alpha <- lambda
  }
  if (length(alpha) != 1L || !is.finite(alpha)) {
    stop("alpha must be a finite scalar.", call. = FALSE)
  }
  alpha <- min(max(alpha[1L], 0), 1)
  alpha_mse <- .ld_mse_shrinkage_intensity(S, n, target = target)
  alpha_used <- min(alpha_mse, alpha)

  if (default_target) {
    out <- (1 - alpha_used) * S
    diag(out) <- diag(out) + alpha_used
    out <- .ld_symmetrize(out)
  } else {
    out <- .ld_symmetrize(alpha_used * target + (1 - alpha_used) * S)
  }
  out <- .ld_fspd(out, eig_min = 0)
  out <- stats::cov2cor(out)
  out[is.na(out)] <- 0
  out <- .ld_symmetrize(out)
  diag(out) <- 1
  out
}
