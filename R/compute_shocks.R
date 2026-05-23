
#' @title Computes posterior draws of shocks
#'
#' @description Each of the draws from the posterior estimation of models from
#' the package \pkg{bvars} is transformed into
#' a draw from the posterior distribution of the structural shocks. 
#' 
#' @param posterior posterior estimation outcome obtained by running the 
#' \code{estimate} function. 
#' 
#' @return An object of class \code{PosteriorShocks}, that is, an \code{NxTxS} array 
#' with attribute \code{PosteriorShocks} containing \code{S} draws of the shocks.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' # simple workflow
#' ############################################################
#' spec = specify_bvar$new(us_macro_chan)        # specify the model
#' burn = estimate(spec, 5)                      # run the burn-in
#' post = estimate(burn, 10)                     # estimate the model
#' shoc = compute_shocks(post)                   # compute shocks
#' plot(shoc)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' us_macro_chan |>
#'   specify_bvar$new() |>
#'   estimate(S = 5) |> 
#'   estimate(S = 10) |> 
#'   compute_shocks() |> plot()
#' 
#' @export
compute_shocks <- function(posterior) {
  UseMethod("compute_shocks", posterior)
}





#' @method compute_shocks PosteriorBVAR
#' 
#' @title Computes posterior draws of  shocks
#'
#' @description Each of the draws from the posterior estimation of the BAVR 
#' model is transformed into a draw from the posterior distribution of the shocks. 
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' \code{PosteriorBVAR} obtained by running the \code{estimate} function.
#' 
#' @return An object of class \code{PosteriorShocks}, that is, an \code{NxTxS} 
#' array with attribute \code{PosteriorShocks} containing \code{S} draws of the shocks.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' # simple workflow
#' ############################################################
#' spec = specify_bvar$new(us_macro_chan)        # specify the model
#' burn = estimate(spec, 5)                      # run the burn-in
#' post = estimate(burn, 10)                     # estimate the model
#' shoc = compute_shocks(post)                   # compute shocks
#' 
#' # workflow with the pipe |>
#' ############################################################
#' us_macro_chan |>
#'   specify_bvar$new() |>
#'   estimate(S = 5) |> 
#'   estimate(S = 10) |> 
#'   compute_shocks() -> shoc
#'   
#' @export
compute_shocks.PosteriorBVAR <- function(posterior) {
  
  posterior_A     = posterior$posterior$A
  Y               = posterior$last_draw$data_matrices$Y
  X               = posterior$last_draw$data_matrices$X
  
  ss              = .Call(`_bvars_compute_shocks`, posterior_A, Y, X)
  
  class(ss)       = "PosteriorShocks"
  S               = dim(ss)[3]     
  dimnames(ss)    = list(rownames(Y), colnames(Y), 1:S)
  
  return(ss)
}
