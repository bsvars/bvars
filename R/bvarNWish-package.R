#' 
#' @title Bayesian Vector Autoregressions with Normal-Wishart Priors
#' 
#' @description Bayesian Vector Autoregressions with Normal-Wishart Priors and 
#' the Matrix Generalized Inverse Gaussian Gibbs sampler.
#' 
#' 
#' @docType package
#' @name bvarNWish-package
#' @aliases bvarNWish-package bvarNWish
#' @useDynLib bvarNWish, .registration = TRUE
#' 
#' @importFrom Rcpp sourceCpp
#' @importFrom bsvars estimate forecast
#' @importFrom R6 R6Class
#' @import bsvars
#' @import RcppArmadillo
#' @import RcppProgress
#' 
#' @note This package is currently in active development. Your comments,
#' suggestions and requests are warmly welcome!
#' 
#' @author Rui Liu \email{rl3023@columbia.edu}, 
#' Andres Ramirez Hassan \email{aramir21@gmail.com} & 
#' Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#' 
#' Hamura, Irie, Sugasawa (2024) Gibbs Sampler for Matrix Generalized Inverse 
#' Gaussian Distributions, Journal of Computational and Graphical Statistics,
#' 33(2), 331--340, <doi:10.1080/10618600.2023.2258186>.
#' 
#' Chan (2020) Large Bayesian VARs: A Flexible Kronecker Error Covariance Structure,
#' Journal of Business and Economic Statistics, 38(1), 68--79,
#' <doi:10.1080/07350015.2018.1451336>.
#' 
'_PACKAGE'
