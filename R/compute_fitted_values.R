
#' @method compute_fitted_values PosteriorBVAR
#' 
#' @title Computes posterior draws from data predictive density
#'
#' @description Each of the draws from the posterior estimation of the BVAR
#' model is transformed into a draw from the data predictive density. 
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' \code{PosteriorBVAR} obtained by running the \code{estimate} function.
#' 
#' @return An object of class \code{PosteriorFitted}, that is, an \code{NxTxS} 
#' array with attribute \code{PosteriorFitted} containing \code{S} draws from 
#' the data predictive density.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' spec = specify_bvar$new(us_macro_chan)        # specify the model
#' burn = estimate(spec, 5)                      # run the burn-in
#' post = estimate(burn, 10)                     # estimate the model
#' fitt = compute_fitted_values(post)            # fitted values
#' 
#' # workflow with the pipe |>
#' ############################################################
#' us_macro_chan |>
#'   specify_bvar$new() |>
#'   estimate(S = 5) |> 
#'   estimate(S = 5) |> 
#'   compute_fitted_values() -> fitt
#' 
#' @export
compute_fitted_values.PosteriorBVAR <- function(posterior) {
  
  Y               = posterior$last_draw$data_matrices$Y
  posterior_A     = posterior$posterior$A
  posterior_Sigma = posterior$posterior$Sigma
  
  N               = dim(posterior_A)[1]
  T               = dim(posterior$last_draw$data_matrices$X)[2]
  S               = dim(posterior_A)[3]
  posterior_sigma2 = sqrt(posterior$posterior$sigma2 * posterior$posterior$lambda)
  X               = posterior$last_draw$data_matrices$X
  
  fv              = .Call(`_bvars_compute_fitted_values`, posterior_A, posterior_Sigma, posterior_sigma2, X)

  class(fv)       = "PosteriorFitted"
  S               = dim(posterior_A)[3]      
  dimnames(fv)    = list(rownames(Y), colnames(Y), 1:S)
  
  return(fv)
}
