#' Sample from an Inverse Wishart distribution using Armadillo
#'
#' @title Random Generation of Inverse Wishart Distributed Matrices
#' @description Generate a \eqn{p \times p} matrix from the Inverse Wishart distribution with parameters:
#' integer \code{p} specifying the dimension, a \eqn{p \times p} positive definite matrix \code{Psi}
#' as the scale matrix, and \code{df} as the degrees of freedom.
#' @usage rinvwishart(Psi, df)
#' @param Psi Positive-definite scale matrix
#' @param df Degrees of freedom
#' @details
#' The Inverse Wishart distribution is the distribution of the inverse of a Wishart-distributed matrix.
#' If \eqn{W \sim \mathcal{W}_p(\text{df}, \Psi^{-1})}, then \eqn{W^{-1} \sim \mathcal{W}_p^{-1}(\text{df}, \Psi)}.
#'
#' The scale matrix \code{Psi} must be symmetric and positive-definite, and the degrees of freedom \code{df}
#' must be an integer greater than or equal to the dimension of \code{Psi}.
#' The output is a \eqn{p \times p} symmetric positive-definite matrix drawn from the Inverse Wishart distribution.
#'
#' @return An Inverse Wishart-distributed matrix
#'
#' @author Rui Liu \email{rl3023@columbia.edu}, 
#' Andres Ramirez Hassan \email{aramir21@gmail.com} & 
#' Tomasz Woźniak \email{wozniak.tom@pm.me}
#'
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical Analysis* (3rd ed.). Wiley.
#'
#' @examples
#' # Define a 3x3 identity scale matrix
#' Psi <- diag(3)
#' df <- 5
#' # Generate a sample from the Inverse Wishart distribution
#' rinvwishart(Psi, df)
#'
#' @param Psi A positive-definite scale matrix.
#' @param df Degrees of freedom (must be >= dimension of Psi).
#' @return A matrix sampled from the inverse Wishart distribution.
#' @export
rinvwishart <- function(Psi, df) {
  stopifnot(
    is.matrix(Psi),
    nrow(Psi) == ncol(Psi),
    is.numeric(df), length(df) == 1, df >= nrow(Psi),
    all(abs(Psi - t(Psi)) < 1e-8),          # symmetry check
    all(eigen(Psi, symmetric = TRUE)$values > 0)  # positive-definiteness
  )
  .Call(`_bvarNWish_do_rinvwishart`, Psi, as.integer(df))
}