#  #####################################################################################
#  R package bvars by Rui Liu, Andres Ramirez Hassan and Tomasz Woźniak Copyright (C) 2026
#
#  This file is part of the R package bvars: Forecasting with Large Bayesian Vector Autoregressions
#
#  The R package bvars is free software: you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 3 or
#  any later version of the License.
#
#  The R package bvars is distributed in the hope that it will be
#  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with the R package bsvars. If that is not the case, please
#  refer to <http://www.gnu.org/licenses/>.
#  #####################################################################################
#
#' @title Bayesian Forecasting with Large Vector Autoregressions
#'
#' @description 
#' Provides fast and efficient procedures for Bayesian estimation and forecasting 
#' using state-of-the-art Vector Autoregressions. This package includes the model 
#' proposed by Chan (2020) <doi:10.1080/07350015.2018.1451336>, that is, a 
#' Bayesian Vector Autoregression with Minnesota priors and a flexible structure 
#' of the error term specification. The latter includes: conditional multivariate 
#' normal or Student’s t distributions, as well as homoskedastic or heteroskedastic 
#' specifications with a common volatility modelled by centred or non-centred 
#' Stochastic Volatility. Additionally, the package facilitates predictive 
#' analyses using density forecasting and forecast-error variance decompositions. 
#' All this is complemented by simple workflows, useful plots and summary 
#' functions, and comprehensive documentation. The 'bvars' package aligns with 
#' R packages 'bsvars' by Woźniak (2024) <doi:10.32614/CRAN.package.bsvars>, 
#' 'bsvarSIGNs' by Wang & Woźniak (2025) <doi:10.32614/CRAN.package.bsvarSIGNs>, 
#' and 'bpvars' by Woźniak (2025) <doi:10.32614/CRAN.package.bpvars> regarding 
#' objects, workflows, and code structure, and they constitute an integrated toolset.
#' 
#' @details 
#' \strong{Models.} 
#' All the BVAR models in this package are specified by two equations, including 
#' the reduced form equation:
#' \deqn{Y = AX + E}
#' where \eqn{Y} is an \code{NxT} matrix of dependent variables, 
#' \eqn{X} is a \code{KxT} matrix of explanatory variables, 
#' \eqn{E} is an \code{NxT} matrix of reduced form error terms, 
#' and \eqn{A} is an \code{NxK} matrix of autoregressive slope coefficients and 
#' parameters on deterministic terms in \eqn{X}.
#' 
#' This package assumes that the error matrix follows a matrix normal distribution:
#' \deqn{E \mid X \sim \mathcal{MN}_{N \times T}(\mathbf{0}, \mathbf{\Sigma}, \mathbf{\Omega})}
#' where \eqn{\Sigma} is the \code{NxN} covariance matrix of the error term at 
#' time \eqn{t}, and \eqn{\Omega} is a \code{TxT} diagonal matrix.
#' 
#' The diagonal elements of \eqn{\Omega} determine the specification of the error
#' term covariance structure. Specifically, the error term at time \eqn{t} follows 
#' the multivariate normal distribution 
#' \deqn{e_t \sim \mathcal{N}_N(\mathbf{0}, \sigma_t^2\lambda_t \mathbf{\Sigma})} 
#' where the scalar processes \eqn{\sigma_t^2} and \eqn{\lambda_t} determine the
#' diagonal elements of \eqn{\Omega}. The process \eqn{\sigma_t^2} specifies 
#' conditional variance and includes three options:
#' \describe{
#'  \item{\eqn{\sigma_t^2 = 1}}{homoskedastic error term}
#'  \item{\eqn{\sigma_t^2}}{estimated and following non-centred stochastic volatility}
#'  \item{\eqn{\sigma_t^2}}{estimated and following centred stochastic volatility}
#' }
#' The process \eqn{\lambda_t} specifies the conditional distribution of the error 
#' term and includes two options:
#' \describe{
#'  \item{\eqn{\lambda_t = 1}}{Gaussian error term specification}
#'  \item{\eqn{\lambda_t}}{estimated and following a priori an inverse gamma 2 
#'        distribution \eqn{\mathcal{IG}2(\nu - 2, \nu)}, where \eqn{\nu > 2} is 
#'        a degrees of freedom parameter}
#' }
#' 
#' \strong{Prior distributions.}
#' The autoregressive matrix \eqn{A} is assigned matrix-variate normal distribution:
#' \deqn{
#' \mathbf{A} \mid \underline{\mathbf{A}}, \mathbf{V}, \boldsymbol{\Sigma}
#' \sim \mathcal{MN}_{N \times K}(\underline{\mathbf{A}}, \boldsymbol{\Sigma}, \mathbf{V})
#' }
#' with the mean matrix \eqn{\underline{\mathbf{A}}}, and covariance matrices 
#' \eqn{\boldsymbol{\Sigma}_{N\times N}} and \eqn{\mathbf{V}_{K\times K}} 
#' defining the row- and column-covariance structures. 
#' 
#' This is complemented by the inverse Wishart prior for the error term covariance \eqn{\boldsymbol{\Sigma}}:
#' \deqn{
#' \boldsymbol{\Sigma} \mid \underline{\mathbf{S}}, \underline{\nu} \sim \mathcal{IW}(\underline{\mathbf{S}}, \underline{\nu})
#' }
#' with the scale matrix \eqn{\underline{\mathbf{S}}} and degrees of freedom \eqn{\underline{\nu}}.
#' 
#' @name bvars-package
#' @aliases bvars-package bvars
#' @docType package
#' @useDynLib bvars, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom bsvars estimate specify_data_matrices compute_variance_decompositions compute_fitted_values
#' @importFrom generics forecast
#' @importFrom R6 R6Class
#' @importFrom RcppTN rtn
#' @importFrom stats quantile sd
#' @import RcppArmadillo
#' @import RcppProgress
#' 
#' @note This package is currently in active development. Your comments,
#' suggestions and requests are warmly welcome!
#' 
#' @author Rui Liu \email{rl3023@columbia.edu}, 
#' Andres Ramirez Hassan \email{aramir21@gmail.com}, 
#' Tomasz Woźniak \email{wozniak.tom@pm.me}

#' @references
#' Chan (2020) Large Bayesian VARs: A Flexible Kronecker Error Covariance Structure,
#' Journal of Business and Economic Statistics, 38(1), 68--79,
#' <doi:10.1080/07350015.2018.1451336>.
#' 
#' @examples
#' # simple workflow
#' ############################################################
#' spec = specify_bvar$new(us_macro_chan)        # specify the model
#' burn = estimate(spec, 5)                      # run the burn-in
#' post = estimate(burn, 10)                     # estimate the model
#' 
#' # workflow with the pipe |>
#' ############################################################
#' us_macro_chan |>
#'   specify_bvar$new() |>
#'   estimate(S = 5) |> 
#'   estimate(S = 10) -> post
#'   
'_PACKAGE'
