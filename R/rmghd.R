#' Sample from the Matrix Generalized Hyperbolic distribution
#'
#' @title Random Generation from Matrix Generalized Hyperbolic Distribution (MGHD)
#' @description Generate a matrix from the Matrix Generalized Hyperbolic distribution sequentially:
#' \eqn{\boldsymbol{\Sigma} \sim \mathcal{MGIG}(\kappa, \Psi, \Gamma)}, followed by
#' \eqn{\mathbf{A} \mid \boldsymbol{\Sigma} \sim \mathcal{MN}(M, U, \boldsymbol{\Sigma})}.
#'
#' @usage rmgh(M, U, aux_Sigma, lambda, Psi, Gamma)
#' @param M The \eqn{n \times p} matrix of means.
#' @param U The \eqn{n \times n} row covariance matrix (must be positive definite).
#' @param aux_Sigma A \eqn{p \times p} positive definite matrix used to initialize the MGIG sampler.
#' @param lambda The scalar shape parameter \eqn{\kappa} of the MGIG distribution.
#' @param Psi A \eqn{p \times p} positive definite matrix parameter for MGIG.
#' @param Gamma A \eqn{p \times p} positive definite matrix parameter for MGIG.
#'
#' @details
#' This function draws a sample from the Matrix Generalized Hyperbolic distribution by first sampling
#' \eqn{\boldsymbol{\Sigma} \sim \mathcal{MGIG}_p(\kappa, \Psi, \Gamma)}, then sampling
#' \eqn{\mathbf{A} \mid \boldsymbol{\Sigma} \sim \mathcal{MN}_{n \times p}(M, U, \boldsymbol{\Sigma})}.
#'
#' All matrix arguments must be symmetric and positive definite where required. The matrix \code{aux_Sigma} must be a valid
#' initial value for Cholesky decomposition in the MGIG sampler.
#'
#' @return A matrix sampled from the Matrix Generalized Hyperbolic distribution.
#'
#' @author Rui Liu \email{rl3023@columbia.edu}, 
#' Andres Ramirez Hassan \email{aramir21@gmail.com}, 
#' Tomasz Woźniak \email{wozniak.tom@pm.me}
#'
#' @examples
#' n <- 5; p <- 3
#' M <- matrix(0, n, p)
#' U <- diag(n)
#' aux_Sigma <- diag(p)
#' lambda <- 1
#' Psi <- diag(p)
#' Gamma <- diag(p)
#' rmgh(M, U, aux_Sigma, lambda, Psi, Gamma)
#'
#' @export
rmgh <- function(M, U, aux_Sigma, lambda, Psi, Gamma) {
  stopifnot(
    is.matrix(M),
    is.matrix(U),
    is.matrix(aux_Sigma),
    is.matrix(Psi),
    is.matrix(Gamma),
    nrow(M) == nrow(U),
    ncol(M) == nrow(aux_Sigma),
    all(abs(U - t(U)) < 1e-8),         # symmetry
    all(abs(Psi - t(Psi)) < 1e-8),
    all(abs(Gamma - t(Gamma)) < 1e-8),
    all(eigen(U, symmetric = TRUE)$values > 0), #positive-definiteness
    all(eigen(Psi, symmetric = TRUE)$values > 0),
    all(eigen(Gamma, symmetric = TRUE)$values > 0)
  )
  .Call(`_bvarNWish_sample_mgh`, M, U, aux_Sigma, lambda, Psi, Gamma)
}