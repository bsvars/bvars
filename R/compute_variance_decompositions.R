
#' @title Computes posterior draws of the forecast error variance decomposition
#'
#' @method compute_variance_decompositions PosteriorBVAR
#' 
#' @description Each of the draws from the posterior estimation of the Vector 
#' Autoregression is transformed into a draw from the posterior distribution of 
#' the forecast error variance decomposition. 
#' 
#' @param posterior posterior estimation outcome obtained by running the \code{estimate} function. 
#' @param horizon a positive integer number denoting the forecast horizon for 
#' the forecast error variance decomposition computations.
#' 
#' @return An object of class \code{PosteriorFEVD}, that is, an \code{NxNx(horizon+1)xS} 
#' array with attribute \code{PosteriorFEVD} containing \code{S} draws of the 
#' forecast error variance decomposition.
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Kilian, L., & Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In: Structural vector autoregressive analysis. Cambridge University Press.
#' 
#' @examples
#' spec = specify_bvar$new(us_macro_chan)        # specify the model
#' burn = estimate(spec, 5)                      # run the burn-in
#' post = estimate(burn, 10)                     # estimate the model
#' fevd = compute_variance_decompositions(post, horizon = 4)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' us_macro_chan |>
#'   specify_bvar$new() |>
#'   estimate(S = 5) |> 
#'   estimate(S = 10) |> 
#'   compute_variance_decompositions(horizon = 4) -> fevd
#'   
#' @export
compute_variance_decompositions.PosteriorBVAR <- function(posterior, horizon) {
  
  posterior_Sigma = posterior$posterior$Sigma
  posterior_A     = posterior$posterior$A
  N               = dim(posterior_A)[1]
  p               = posterior$last_draw$p
  S               = dim(posterior_A)[3]
  Y               = posterior$last_draw$data_matrices$Y
  
  qqq             = .Call(`_bvars_compute_variance_decompositions`, posterior_Sigma, posterior_A, horizon, p)
  
  fevd            = array(NA, c(N, N, horizon + 1, S), dimnames = list(rownames(Y), rownames(Y), 0:horizon, 1:S))
  for (s in 1:S) fevd[,,,s] = qqq[s][[1]]
  class(fevd)     = "PosteriorFEVD"
  
  return(fevd)
}
