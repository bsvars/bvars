
#' @title Bayesian estimation of a Bayesian Vector Autoregression with a Flexible 
#' Prior Distribution via Gibbs sampler
#'
#' @description Estimates a Bayesian Vector Autoregression with a flexible prior
#' specification as proposed by Liu, Ramirez Hassan, and Woźniak. The prior 
#' distribution for the autoregressive parameter matrix \eqn{A} is matrix-variate 
#' normal, error term  covariance \eqn{\Sigma} is inverse Wishart, both as 
#' defined by Woźniak (2016). Additionally, the prior equation-specific 
#' covariance of autoregressive parameters, denoted by \eqn{V} follows a matrix
#' generalised inverse Gaussian distribution as proposed by Barndorff-Nielsen, 
#' Blaesild, Jensen, Jorgensen (1982). The resulting marginal prior distribution
#' for the autoregressive parameters \eqn{A} is the matrix-variate generalized 
#' hyperbolic by Thabane, Safiul Haq (2004). The estimation proceeds by Gibbs 
#' sampler where at each iteration parameters \eqn{A} and \eqn{\Sigma} are sampled
#' from the matrix-variate normal inverse Wishart distribution as in Chan (2020), 
#' while the hyper-parameter matrix \eqn{V} is sampled from the matrix generalised 
#' inverse Gaussian distribution as in Hamura, Irie, Sugasawa (2024). 
#' See section \bold{Details} for the model equations.
#' 
#' @details 
#' The Bayesian Vector Autoregression is given by equation:
#' \deqn{Y = AX + E}
#' where \eqn{Y} is an \code{NxT} matrix of dependent variables, \eqn{X} is a 
#' \code{KxT} matrix of explanatory variables, \eqn{E} is an \code{NxT} matrix 
#' of error terms, and \eqn{A} is an \code{NxK} matrix of autoregressive 
#' coefficients and parameters on deterministic terms in \eqn{X}.The error terms, 
#' \code{E}, are temporally and contemporaneously independent and jointly 
#' normally distributed with zero mean and period-specific covariance \eqn{\Sigma}.
#' 
#' @param specification an object of class \code{BVAR} generated using the 
#' \code{specify_bvar$new()} function.
#' @param S a positive integer, the number of posterior draws to be generated
#' @param thin a positive integer, specifying the frequency of MCMC output thinning
#' @param show_progress a logical value, if \code{TRUE} the estimation progress 
#' bar is visible
#' 
#' @return An object of class \code{PosteriorBVAR} containing the Bayesian 
#' estimation output and containing two elements:
#' 
#'  \code{posterior} a list with a collection of \code{S} draws from the 
#'  posterior distribution generated via Gibbs sampler containing:
#'  \describe{
#'  \item{A}{an \code{NxKxS} array with the posterior draws for matrix \eqn{A}}
#'  \item{Sigma}{an \code{NxNxS} array with the posterior draws for matrix \eqn{\Sigma}}
#'  \item{V}{a \code{KxKxS} array with the posterior draws for the hyper-parameter
#'  matrix \eqn{\Sigma}}
#' }
#' 
#' \code{last_draw} an object of class \code{BVAR} with the last draw of the 
#' current MCMC run as the starting value to be passed to the continuation of 
#' the MCMC estimation using \code{estimate()}. 
#'
#' @seealso \code{\link{specify_bvar}}, \code{\link{specify_posterior_bvar}}
#'
#' @author Rui Liu \email{rl3023@columbia.edu}, 
#' Andres Ramirez Hassan \email{aramir21@gmail.com} & 
#' Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Barndorff-Nielsen, Blaesild, Jensen, Jorgensen (1982) Exponential 
#' Transformation Models, Proceedings of the Royal Society of London. A. 
#' Mathematical and Physical Sciences, 379, 41–-65, <doi:10.1098/rspa.1982.0004>.
#' 
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
#' Woźniak (2016) Bayesian Vector Autoregressions, Australian Economic Review,
#' 49(3), 365--380, <doi:10.1111/1467-8462.12179>.
#' 
#' @method estimate BVAR
#' 
#' @examples
#' # simple workflow
#' ############################################################
#' spec = specify_bvar$new(us_macro_chan)   # specify the model
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
#' @export
estimate.BVAR <- function(specification, S, thin = 1, show_progress = TRUE) {
  
  # get the inputs to estimation
  prior               = specification$get_prior()
  starting_values     = specification$get_starting_values()
  data_matrices       = specification$get_data_matrices()
  centred_sv          = specification$get_centred_sv()
  normal              = specification$get_normal()
  homoskedastic       = specification$get_homoskedastic()
  N                   = dim(data_matrices$Y)[1]
  aux_mix             = matrix(NA, 3, 10)
  if (N <= 100) {aux_mix = sv_aux_mix[,,N]}
  
  # estimation
  qqq   = .Call(
    `_bvars_bvar_mgig_cpp`, 
    S, 
    data_matrices$Y, 
    data_matrices$X, 
    prior, 
    starting_values, 
    aux_mix,
    homoskedastic, 
    centred_sv,
    normal,
    thin, 
    show_progress
  )
  
  # prepare the output
  specification$starting_values$set_starting_values(qqq$last_draw)
  output              = specify_posterior_bvar$new(specification, qqq$posterior)
  
  return(output)
}


#' @inherit estimate.BVAR
#' 
#' @method estimate PosteriorBVAR
#' 
#' @param specification an object of class \code{PosteriorBVAR} generated 
#' using the \code{estimate.BVAR()} function. This setup facilitates the 
#' continuation of the MCMC sampling starting from the last draw of the previous 
#' run.
#' 
#' @export
estimate.PosteriorBVAR <- function(specification, S, thin = 1, show_progress = TRUE) {
  
  # get the inputs to estimation
  prior               = specification$last_draw$prior$get_prior()
  starting_values     = specification$last_draw$starting_values$get_starting_values()
  data_matrices       = specification$last_draw$data_matrices$get_data_matrices()
  centred_sv          = specification$last_draw$get_centred_sv()
  normal              = specification$last_draw$get_normal()
  homoskedastic       = specification$last_draw$get_homoskedastic()
  N                   = dim(data_matrices$Y)[1]
  aux_mix             = matrix(NA, 3, 10)
  if (N <= 100) {aux_mix = sv_aux_mix[,,N]}
  
  # estimation
  qqq   = .Call(
    `_bvars_bvar_mgig_cpp`, 
    S, 
    data_matrices$Y, 
    data_matrices$X, 
    prior, 
    starting_values, 
    aux_mix,
    homoskedastic, 
    centred_sv,
    normal,
    thin, 
    show_progress
  )
  
  # prepare the output
  specification$last_draw$starting_values$set_starting_values(qqq$last_draw)
  output              = specify_posterior_bvar$new(specification$last_draw, qqq$posterior)
  
  return(output)
}