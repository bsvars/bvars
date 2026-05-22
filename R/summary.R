
#' @title Provides posterior summary of  VAR estimation
#'
#' @description Provides posterior mean, standard deviations, as well as 5 and 95 
#' percentiles of the parameters: autoregressive 
#' parameters \eqn{\mathbf{A}}, and the covariance matrix \eqn{\Sigma}.
#' 
#' @param object an object of class PosteriorBVAR obtained using the
#' \code{estimate()} function applied to a Bayesian VAR
#' model specification set by function \code{specify_bvar$new()} containing 
#' draws from the  posterior distribution of the parameters. 
#' @param ... additional arguments affecting the summary produced.
#' 
#' @return A list reporting the posterior mean, standard deviations, as well as 5 and 95 
#' percentiles of the parameters: autoregressive 
#' parameters \eqn{\mathbf{A}}, and the covariance matrix \eqn{\Sigma}..
#' 
#' @method summary PosteriorBVAR
#' 
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' # simple workflow
#' ############################################################
#' spec = specify_bvar$new(us_macro_chan)        # specify the model
#' burn = estimate(spec, 5)                      # run the burn-in
#' post = estimate(burn, 10)                     # estimate the model
#' summary(post)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' us_macro_chan |>
#'   specify_bvar$new() |>
#'   estimate(S = 5)  |> 
#'   estimate(S = 10) |> 
#'   summary()
#' 
#' @export
summary.PosteriorBVAR = function(
    object,
    ...
) {
  
  N         = dim(object$posterior$A)[1]
  p         = object$last_draw$p
  K         = dim(object$last_draw$data_matrices$X)[1]
  d         = K - N * p
  
  out       = list()
  out$A     = list()
  out$Sigma = list()
  
  param     = c("A", "Sigma")
  
  for (n in 1:N) {
    
    Sigma       = matrix(object$posterior$Sigma[n,1:n,], nrow = n)
    out$Sigma[[n]] = matrix(
      cbind(
        apply(Sigma, 1, mean),
        apply(Sigma, 1, sd),
        apply(Sigma, 1, quantile, probs = 0.05),
        apply(Sigma, 1, quantile, probs = 0.95)
      ),
      ncol = 4
    )
    colnames(out$Sigma[[n]]) = c("mean", "sd", "5% quantile", "95% quantile")
    rownames(out$Sigma[[n]]) = paste0("Sigma[", n, ",", 1:n, "]")  
    
    out$A[[n]] = cbind(
      apply(object$posterior$A[n,,], 1, mean),
      apply(object$posterior$A[n,,], 1, sd),
      apply(object$posterior$A[n,,], 1, quantile, probs = 0.05),
      apply(object$posterior$A[n,,], 1, quantile, probs = 0.95)
    )
    colnames(out$A[[n]]) = c("mean", "sd", "5% quantile", "95% quantile")
    
    Anames  = c(
      paste0(
        rep("lag", p * N),
        kronecker((1:p), rep(1, N)),
        rep("_var", p * N),
        kronecker((1:N), rep(1, p))
      ),
      "const"
    )
    if (d > 1) {
      Anames = c(Anames, paste0("exo", 1:(d - 1)))
    }
    rownames(out$A[[n]]) = Anames
  } # END n loop
  
  names(out$Sigma)  = paste0("equation", 1:N)
  names(out$A)      = paste0("equation", 1:N)
  
  return(out)
} # END summary.PosteriorBVAR
