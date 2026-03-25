#' Sample from a Wishart distribution using Armadillo
#'
#' @title Random Generation of Wishart Distributed Matrices
#' @description Generate a \eqn{p \times p} matrix from the Wishart distribution with parameters:
#' integer \code{p} specifying the dimension, a \eqn{p \times p} positive definite matrix \code{Sigma}
#' as the scale matrix, and \code{df} as the degrees of freedom.
#' @usage rwishart(Sigma, df) 
#' @param Sigma Positive-definite scale matrix
#' @param df Degrees of freedom
#' @details
#' The Wishart distribution is a multivariate generalization of the chi-squared distribution
#' and is defined over symmetric positive-definite matrices. If \code{X_1, ..., X_{n}} are
#' independent and identically distributed \eqn{p}-dimensional random vectors from a
#' multivariate normal distribution with zero mean and covariance matrix \code{Sigma},
#' then the matrix
#' \deqn{W = \sum_{i=1}^{\text{n}} X_i X_i^\top}
#' follows a Wishart distribution with parameters \code{n} and \code{Sigma},
#' denoted by \eqn{W \sim \mathcal{W}_p(\text{n}, \Sigma)}.
#'
#' The scale matrix \code{Sigma} must be symmetric
#' and positive-definite, and the degrees of freedom \code{df} must be an integer
#' greater than or equal to the dimension of \code{Sigma}.
#' The output is a \eqn{p \times p} symmetric positive-definite matrix.
#' @author Rui Liu \email{rl3023@columbia.edu}, 
#' Andres Ramirez Hassan \email{aramir21@gmail.com} & 
#' Tomasz Woźniak \email{wozniak.tom@pm.me}
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical Analysis* (3rd ed.). Wiley.
#' @examples
#' # Define a 3x3 identity scale matrix
#' Sigma <- diag(3)
#' df <- 5
#' # Generate a sample from the Wishart distribution
#' rwishart(Sigma, df)

#' @return A Wishart-distributed matrix
#' @export
rwishart <- function(Sigma, df) {
  stopifnot(
    is.matrix(Sigma),
    nrow(Sigma) == ncol(Sigma),
    is.numeric(df), length(df) == 1, df >= nrow(Sigma)
  )
  .Call(`_bvars_do_rwishart`, Sigma, as.integer(df))
}
