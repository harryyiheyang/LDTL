#' POET with linear shrinkage of the idiosyncratic covariance
#'
#' Estimate a factor component from the leading eigenpairs of `S` and shrink the
#' idiosyncratic covariance toward its diagonal.
#'
#' @param S Optional symmetric covariance or correlation matrix.
#' @param X Optional individual-level data matrix. Used to compute `S` when
#'   `S` is not supplied and to infer `n`.
#' @param scale If `TRUE`, center and standardize `X` before computing `S`.
#'   Default is `FALSE`, assuming `X` has already been scaled.
#' @param n Sample size used to form `S`.
#' @param alpha Shrinkage intensity in `[0, 1]` toward the diagonal of the
#'   idiosyncratic covariance. Default is 0.5.
#'
#' @return A positive-semidefinite POET regularized correlation matrix.
#' @export
poet_linear_shrinkage <- function(
    S = NULL,
    n = NULL,
    alpha = 0.5,
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
  E <- comp$E
  D <- diag(diag(E), nrow(E), ncol(E))

  alpha <- min(max(alpha[1L], 0), 1)
  E_reg <- (1 - alpha) * E + alpha * D
  E_reg <- .ld_symmetrize(E_reg)
  E_reg <- .ld_fspd(E_reg, eig_min = 0)

  Shat <- .ld_symmetrize(comp$P + E_reg)
  stats::cov2cor(Shat)
}
