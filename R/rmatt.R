#' Sample from an matrix variate student's t distribution using Armadillo
#'
#' @title Random Generation of Matrix Variate Student's T Distributed Matrices
#' @description Generate a \eqn{p \times n} random matrix from the matrix-variate t-distribution with parameters:
#' a \eqn{p \times n} mean matrix \code{M}, a \eqn{p \times p} row covariance matrix \code{V}, a \eqn{n \times n} column covariance matrix \code{S}, and degrees of freedom \code{df}.
#' The matrix-variate t-distribution is a generalization of the multivariate t-distribution to random matrices, and is often denoted as:
#' \deqn{X \sim \mathcal{T}_{p,n}(\nu, M, V, S)}
#' meaning \eqn{\text{vec}(X)} follows a multivariate t-distribution with mean \eqn{\text{vec}(M)} and covariance matrix \eqn{S \otimes V}.
#' @usage rmatt(M, V, S, df)
#' @param M Mean matrix of size \eqn{p \times n}
#' @param V Row covariance matrix (positive-definite, \eqn{p \times p})
#' @param S Column covariance matrix (positive-definite, \eqn{n \times n})
#' @param df Degrees of freedom (positive integer)
#' @details
#' The matrix-variate t-distribution is a generalization of the multivariate t-distribution to matrix-valued random variables.
#' A random matrix \eqn{X} has a matrix-variate t-distribution with parameters \eqn{M}, \eqn{V}, \eqn{S}, and \eqn{\nu} if its density is proportional to
#' \deqn{ |V|^{-n/2} |S|^{-p/2} \left[1 + \text{tr}\left(V^{-1}(X - M) S^{-1} (X - M)^\top \right)/\nu \right]^{-(\nu + pn)/2} }
#' where \eqn{X} is a \eqn{p \times n} matrix, \eqn{M} is the mean matrix, \eqn{V} is the row covariance matrix,
#' \eqn{S} is the column covariance matrix, and \eqn{\nu} is the degrees of freedom.
#'
#' A standard method to generate a draw from the matrix-variate t-distribution \eqn{\mathcal{T}_{p,n}(M, V, S, \nu)} is as follows:
#' \enumerate{
#'   \item Generate a scalar \eqn{w \sim \chi^2(\nu)}.
#'   \item Generate a matrix \eqn{Z \sim \mathcal{MN}(0, V, S)}, i.e., \eqn{Z = L_V^\top E L_S}, where \eqn{E} is a matrix of standard normal random variables,
#'         and \eqn{L_V}, \eqn{L_S} are the Cholesky factors of \eqn{V} and \eqn{S}, respectively.
#'   \item Return \eqn{X = M + Z / \sqrt{w / \nu}}.
#' }
#' This yields a matrix \eqn{X} with a matrix-variate t-distribution with \eqn{\nu} degrees of freedom.
#'
#' Both \code{V} and \code{S} must be symmetric positive-definite matrices, and \code{df} (degrees of freedom) must be a scalar greater than 0.

#'
#' @return A \eqn{p \times n} matrix drawn from the matrix-variate t-distribution
#'
#' @author Rui Liu \email{rl3023@columbia.edu}, 
#' Andres Ramirez Hassan \email{aramir21@gmail.com} & 
#' Tomasz Woźniak \email{wozniak.tom@pm.me}
#'
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical Analysis* (3rd ed.). Wiley.
#'
#' @examples
#' # Define dimensions
#' p <- 3  # number of rows
#' n <- 4  # number of columns
#' df <- 5 # degrees of freedom
#'
#' # Mean matrix
#' M <- matrix(0, nrow = p, ncol = n)
#'
#' # Row and column covariance matrices
#' V <- diag(p)
#' S <- diag(n)
#'
#' # Generate a sample from the matrix-variate t-distribution
#' rmatt(M, V, S, df)
#'
#' @export
rmatt <- function(M, V, S, df) {
  if (!is.matrix(M)) stop("M must be a matrix.")
  if (!is.matrix(V) || !is.matrix(S)) stop("V and S must be matrices.")
  if (nrow(M) != nrow(V)) stop("Number of rows in M must match dimension of V.")
  if (ncol(M) != nrow(S)) stop("Number of columns in M must match dimension of S.")
  if (nrow(V) != ncol(V)) stop("V must be square.")
  if (nrow(S) != ncol(S)) stop("S must be square.")
  if (max(abs(V - t(V))) > 1e-8) stop("V must be symmetric.")
  if (max(abs(S - t(S))) > 1e-8) stop("S must be symmetric.")
  if (any(eigen(V, symmetric = TRUE)$values <= 0)) stop("V must be positive-definite.")
  if (any(eigen(S, symmetric = TRUE)$values <= 0)) stop("S must be positive-definite.")
  if (!is.numeric(df) || length(df) != 1 || df <= 0) stop("df must be a positive scalar.")
  
  .Call(`_bvarNWish_do_rmatt`, M, V, S, as.numeric(df))
}