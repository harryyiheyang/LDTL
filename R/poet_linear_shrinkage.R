#' POET with linear shrinkage of the idiosyncratic covariance
#'
#' Estimate a factor component from the leading eigenpairs of `S` and shrink the
#' idiosyncratic covariance toward its diagonal. The selected spike eigenvalues
#' receive the S-POET positive-part bias correction.
#'
#' @param S Optional symmetric covariance or correlation matrix.
#' @param X Optional individual-level data matrix. Used to compute `S` when
#'   `S` is not supplied and to infer `n`.
#' @param scale If `TRUE`, center and standardize `X` before computing `S`.
#'   Default is `FALSE`, assuming `X` has already been scaled.
#' @param n Sample size used to form `S`.
#' @param alpha Boundary in `[0, 1]` for the MSE plug-in shrinkage intensity.
#'   The applied intensity is the smaller of the MSE rule and `alpha`. Default
#'   is 0.05.
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
#' @return A positive-semidefinite POET regularized correlation matrix.
#' @export
poet_linear_shrinkage <- function(
    S = NULL,
    n = NULL,
    alpha = 0.05,
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
  E <- comp$E
  D <- diag(diag(E), nrow(E), ncol(E))

  if (length(alpha) != 1L || !is.finite(alpha)) {
    stop("alpha must be a finite scalar.", call. = FALSE)
  }
  alpha <- min(max(alpha[1L], 0), 1)
  alpha_mse <- .ld_mse_shrinkage_intensity(E, n)
  alpha_used <- min(alpha_mse, alpha)
  E_reg <- (1 - alpha_used) * E + alpha_used * D
  E_reg <- .ld_symmetrize(E_reg)
  E_reg <- .ld_fspd(E_reg, eig_min = 0)

  Shat <- .ld_symmetrize(comp$P + E_reg)
  stats::cov2cor(Shat)
}
