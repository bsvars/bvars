#  #####################################################################################
#  R package bvars by Rui Liu, Andres Ramirez Hassan and Tomasz Woźniak Copyright (C) 2026
#
#  This file is part of the R package bvars: Bayesian Vector Autoregressions with 
#  Flexible Priors
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
#' @title Bayesian Vector Autoregressions with Flexible Priors
#'
#' @description Provides fast and efficient procedures for Bayesian analysis of 
#' vector Autoregressions with flexible priors. This package supports a hierarchical structure for the 
#' location parameters in BVAR models. In particular, this package allows the row covariance structure of 
#' the location parameters to follow a matrix generalized inverse Gaussian distribution. 
#' This specification is flexible enough to encompass many previous approaches aimed at increasing
#' flexibility in BVAR models, while also allowing for a higher degree of shrinkage. 
#' Moreover, this package supports a Gibbs sampling algorithm that facilitates posterior inference.
#' 
#' @details 
#' \strong{Models.} All the BVAR models in this package are specified by two equations, including 
#' the reduced form equation:
#' \deqn{Y = AX + E}
#' where \eqn{Y} is an \code{TxN} matrix of dependent variables, 
#' \eqn{X} is a \code{TxK} matrix of explanatory variables, 
#' \eqn{E} is an \code{TxN} matrix of reduced form error terms, 
#' and \eqn{A} is an \code{KxN} matrix of autoregressive slope coefficients and 
#' parameters on deterministic terms in \eqn{X}.
#' 
#' This package assumes that the error matrix follows a matrix normal distribution:
#' \eqn{E \mid X \sim \mathcal{MN}_{T \times N}(\mathbf{0}, \mathbf{\Omega}, \mathbf{\Sigma})}
#' where \eqn{\Omega} and \eqn{\Sigma} are the row and column covariance matrices
#' 
#' The likelihood function is proportional to
#' \deqn{
#' L\left(\mathbf{Y} \mid \mathbf{X} \mathbf{A}, \boldsymbol{\Omega}, \boldsymbol{\Sigma} \right)
#' \propto \det\left( \boldsymbol{\Sigma} \right)^{-T/2}
#' \det\left( \boldsymbol{\Omega} \right)^{-N/2}
#' \exp\left\{
#'   -\frac{1}{2} \mathrm{tr}\left[
#'     \boldsymbol{\Sigma}^{-1} \left( \mathbf{Y} - \mathbf{X} \mathbf{A} \right)^{\top}
#'     \boldsymbol{\Omega}^{-1} \left( \mathbf{Y} - \mathbf{X} \mathbf{A} \right)
#'   \right]
#' \right\}
#' }
#' 
#' This package adopts the following Bayesian specification to perform inference in this model. 
#' In particular, this package assumes that the serial covariance structure has an independent prior 
#' that follows a stochastic volatility process (Carriero et al., 2016; Chan, 2020a): 
#' 
#' \deqn{
#' \boldsymbol{\Omega} = \operatorname{diag}(\omega_1, \omega_2, \dots, \omega_T), \\
#' \omega_t = \lambda_t \sigma_t^2, \\
#' \lambda_t \sim \mathcal{IG}(\alpha_0, \alpha_0 + 1), \\
#' \sigma_t^2 = \exp\{\eta h_t\}, \\
#' h_t = g h_{t-1} + v_t, \\
#' v_t \sim \mathcal{N}(0,1),
#' }
#'
#' where \eqn{\mathbb{E}(\lambda_t) = 1} by construction.
#' 
#' This package further assumes that 
#' \deqn{
#' \mathbf{A} \mid \underline{\mathbf{A}}, \mathbf{V}, \boldsymbol{\Sigma}
#' \sim \mathcal{MN}_{K \times N}(\underline{\mathbf{A}}, \mathbf{V}, \boldsymbol{\Sigma})
#' }
#' that is, the location parameters follow a matrix normal distribution with mean matrix 
#' \eqn{\underline{\mathbf{A}}}, and covariance matrices \eqn{\mathbf{V}_{K\times K}} and 
#' \eqn{\boldsymbol{\Sigma}_{N\times N}} defining the row and column covariance structures, respectively. (Koop and Korobilis, 2010; Chan, 2020b). 
#' 
#' \strong{Prior distributions.}
#' This package supports many specifications of the priors, namely: 
#'  \itemize{
#'  \item The column covariance structure \eqn{\boldsymbol{\Sigma} \sim \mathcal{IW}_N(\underline{\mathbf{S}}, \underline{\nu})}
#'  \item The row covariance structure \eqn{\mathbf{V} \sim \mathcal{MGIG}_K(\underline{\boldsymbol{\Lambda}}, \underline{\boldsymbol{\Psi}}, \underline{\kappa})} 
#'  }
#' 
#' Note that, the MGIG distribution encompasses several well-known distributions. 
#' In particular, it generalizes the univariate generalized inverse Gaussian distribution; 
#' specifically: 
#' \itemize{ 
#' \item \eqn{ K = 1 } recovers the scalar GIG case 
#' \item If \eqn{\underline{\boldsymbol{\Psi}} = \mathbf{0}} and \eqn{ \underline{\kappa} > (K - 1)/2 }
#' then \eqn{\mathbf{V} \sim \mathcal{W}(\underline{\boldsymbol{\Lambda}}, \underline{\kappa})}
#' \item If \eqn{ \underline{\boldsymbol{\Lambda}} = \mathbf{0}}
#' and \eqn{\underline{\kappa} < -(K - 1)/2 } 
#' then \eqn{ \mathbf{V} \sim \mathcal{IW}(\underline{\boldsymbol{\Psi}}, \underline{\kappa})}
#' }
#' 
#' @name bvars-package
#' @aliases bvars-package bvars
#' @docType package
#' @useDynLib bvars, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom bsvars estimate forecast
#' @importFrom R6 R6Class
#' @import bsvars
#' @import RcppArmadillo
#' @import RcppProgress
#' @note This package is currently in active development. Your comments,
#' suggestions and requests are warmly welcome!
#' @author Rui Liu \email{rl3023@columbia.edu}, 
#' Andres Ramirez Hassan \email{aramir21@gmail.com}, 
#' Tomasz Woźniak \email{wozniak.tom@pm.me}

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
#' !!!!!!to be updated!!!!!!!!!!!
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
