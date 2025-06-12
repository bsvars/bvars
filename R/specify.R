
#' R6 Class Representing PriorBSVAR
#'
#' @description
#' The class PriorBSVAR presents a prior specification for the homoskedastic bsvar model.
#' 
#' @examples 
#' prior = specify_prior_bsvar$new(N = 3, p = 1)  # a prior for 3-variable example with one lag
#' prior$A                                        # show autoregressive prior mean
#' 
#' @export
specify_prior_bvarGIG = R6::R6Class(
  "PriorBVARGIG",
  
  public = list(
    
    #' @field A a real-valued \code{KxN} matrix, the mean matrix \eqn{A_0} of 
    #' the matrix-variate normal prior distribution for the parameter 
    #' matrix \eqn{A}. 
    A          = matrix(),
    
    #' @field S a \code{NxN} positive definite scale matrix \eqn{S_0} of the 
    #' Inverse Wishart prior distribution for the error term covariance 
    #' matrix \eqn{\Sigma}.
    S    = matrix(),
    
    #' @field nu a positive scalar, shape parameter \eqn{\nu_0} of the Inverse 
    #' Wishart prior distribution for the error term covariance 
    #' matrix \eqn{\Sigma}.
    nu   = numeric(),
    
    #' @field Psi a \code{KxK} scale matrix \eqn{\Psi_0} of the matrix 
    #' generalized inverse Gaussian distribution for the column-specific prior 
    #' covariance \eqn{V}
    Psi   = matrix(),
    
    #' @field Gamma a \code{KxK} scale matrix \eqn{\Gamma_0} of the matrix 
    #' generalized inverse Gaussian distribution for the column-specific prior 
    #' covariance \eqn{V}
    Gamma   = matrix(),
    
    #' @field lambda a positive scalar shape parameter \eqn{\lambda_0} of the 
    #' matrix generalized inverse Gaussian distribution for the column-specific 
    #' prior covariance \eqn{V}
    lambda   = numeric(),
    
    #' @description
    #' Create a new prior specification \code{PriorBVARGIG}.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the VAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @param stationary an \code{N} logical vector - its element set to 
    #' \code{FALSE} sets the prior mean for the autoregressive parameters of the 
    #' \code{N}th equation to the white noise process, otherwise to random walk.
    #' @return A new prior specification \code{PriorBVARGIG}.
    #' @examples 
    #' # a prior for 3-variable example with one lag and stationary data
    #' prior = specify_prior_bvarGIG$new(N = 3, p = 1, stationary = rep(TRUE, 3))
    #' prior$A # show autoregressive prior mean
    #' 
    initialize = function(N, p, d = 0, stationary = rep(FALSE, N)){
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      stopifnot("Argument stationary must be a logical vector of length N." = length(stationary) == N & is.logical(stationary))
      
      K             = N * p + 1 + d
      self$A        = cbind(diag(as.numeric(!stationary)), matrix(0, N, K - N))
      self$S        = diag(N)
      self$nu       = N + 1
      self$Psi      = diag(N)
      self$Gamma    = diag(N)
      self$lambda   = N + 1
    }, # END initialize
    
    #' @description
    #' Returns the elements of the prior specification \code{PriorBVARGIG} as 
    #' a \code{list}.
    #' 
    #' @examples 
    #' # a prior for 3-variable example with four lags
    #' prior = specify_prior_bvarGIG$new(N = 3, p = 4)
    #' prior$get_prior() # show the prior as list
    #' 
    get_prior = function(){
      list(
        A        = self$A,
        S        = self$S,
        nu       = self$nu,
        Psi      = self$Psi,
        Gamma    = self$Gamma,
        lambda   = self$lambda
      )
    } # END get_prior
    
  ) # END public
) # END specify_prior_bvarGIG
