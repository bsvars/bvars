

#' @title Samples random numbers from the matrix generalised inverse Gaussian distribution
#' @param Sigma a positive definite symmetric matrix
#' @param lambda a positive scalar
#' @param Psi a positive definite symmetric matrix
#' @param Gamma a positive definite symmetric matrix
#' @return  a matrix a draw from the matrix generalised inverse Gaussian distribution
#' @examples
#' rmgig1(diag(2), 1, diag(2), diag(2))
#' @export
rmgig1 <- function(Sigma, lambda, Psi, Gamma) {
  
  # Check if Sigma, Psi, and Gamma are positive definite
  if (!is.matrix(Sigma) || !is.matrix(Psi) || !is.matrix(Gamma)) {
    stop("Sigma, Psi, and Gamma must be matrices.")
  }
  if (nrow(Sigma) != ncol(Sigma) || nrow(Psi) != ncol(Psi) || nrow(Gamma) != ncol(Gamma)) {
    stop("Sigma, Psi, and Gamma must be square matrices.")
  }
  if (nrow(Sigma) != nrow(Psi) || nrow(Sigma) != nrow(Gamma)) {
    stop("Sigma, Psi, and Gamma must have the same dimensions.")
  }
  if (!isSymmetric(Sigma) || !isSymmetric(Psi) || !isSymmetric(Gamma)) {
    stop("Sigma, Psi, and Gamma must be symmetric matrices.")
  }
  if (any(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values <= 0) ||
      any(eigen(Psi, symmetric = TRUE, only.values = TRUE)$values <= 0) ||
      any(eigen(Gamma, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("Sigma, Psi, and Gamma must be positive definite matrices.")
  }
  out = .Call(`_bvarNWish_do_rmgig1`, Sigma, lambda, Psi, Gamma)
  return(out)
}
