#' 
#' @title Bayesian Vector Autoregressions with Flexible Priors
#' 
#' @description Bayesian Vector Autoregressions with Normal-Matrix Generalised 
#' Inverse Gaussian Priors.
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
#' Chan (2020) Large Bayesian VARs: A Flexible Kronecker Error Covariance Structure,
#' Journal of Business and Economic Statistics, 38(1), 68--79,
#' <doi:10.1080/07350015.2018.1451336>.
#' 
#' Hamura, Irie, Sugasawa (2024) Gibbs Sampler for Matrix Generalized Inverse 
#' Gaussian Distributions, Journal of Computational and Graphical Statistics,
#' 33(2), 331--340, <doi:10.1080/10618600.2023.2258186>.
#' 
#' Thabane, Safiul Haq (2004) On the Matrix-Variate Generalized Hyperbolic 
#' Distribution and Its Bayesian Applications, Statistics: A Journal of Theoretical 
#' and Applied Statistics, 38(6), 511--526, <doi:10.1080/02331880412331319279>.
#' 
#' @examples
#' # simple workflow
#' ############################################################
#' # upload data
#' data(us_fiscal_lsuw)
#' 
#' # specify the model and set seed
#' specification  = specify_bvarGIG$new(us_fiscal_lsuw, p = 1)
#' set.seed(123)
#' 
#' # run the burn-in
#' burn_in        = estimate(specification, 5)
#' 
#' # estimate the model
#' posterior      = estimate(burn_in, 10)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' us_fiscal_lsuw |>
#'   specify_bvarGIG$new(p = 1) |>
#'   estimate(S = 5) |> 
#'   estimate(S = 10) -> posterior
#'   
'_PACKAGE'
