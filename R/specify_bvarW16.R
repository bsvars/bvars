
bvar_posterior = function(
    Y,
    X,
    A,
    V,
    V_inv,
    S,
    nu
) {
  
  post_v_inv  = solve(V)
  post_V_inv  = post_v_inv + tcrossprod(X)
  post_V      = solve(post_V_inv)
  post_A      = (A %*% post_V_inv + Y %*% t(X)) %*% post_V
  post_S      = S + tcrossprod(Y) + A %*% V_inv %*% t(A) - post_A %*% post_V_inv %*% t(post_A)
  post_nu     = nu + dim(Y)[2]
  
  return(list(
    A     = post_A,
    V     = post_V,
    V_inv = post_V_inv,
    S     = post_S,
    nu    = post_nu
  ))
}


#' R6 Class Representing \code{PriorBVARW16}
#'
#' @description
#' The class \code{PriorBVARW16} presents a prior specification for the BVAR model.
#'
#' @examples
#' prior = specify_prior_bvarW16$new(us_macro_chan, p = 4)
#' prior$A
#'
#' @export
specify_prior_bvarW16 = R6::R6Class(
  "PriorBVARW16",

  public = list(

    #' @field A a real-valued \code{NxK} matrix, the mean matrix \eqn{A_0} of
    #' the matrix-variate normal prior distribution for the parameter
    #' matrix \eqn{A}.
    A          = matrix(),

    #' @field S a \code{NxN} positive definite scale matrix \eqn{S_0} of the
    #' Inverse Wishart prior distribution for the error term covariance
    #' matrix \eqn{\Sigma}.
    S    = matrix(),

    #' @field V a \code{KxK} positive definite row-specific covariance matrix \eqn{V_0} of the
    #' Matrix-variate Normal prior distribution for the autoregressive matrix
    #' \eqn{A}.
    V    = matrix(),
    
    #' @field V_inv a \code{KxK} inverse row-specific covariance matrix of the
    #' Matrix-variate Normal prior distribution for the autoregressive matrix
    #' \eqn{A}.
    V_inv    = matrix(),
    
    #' @field nu a positive scalar, shape parameter \eqn{\nu_0} of the Inverse
    #' Wishart prior distribution for the error term covariance
    #' matrix \eqn{\Sigma}.
    nu   = numeric(),

    #' @field lambda a \code{4}-vector of Minnesota Prior hyper-parameters.
    lambda   = numeric(),

    #' @field Y_plus an \code{NxT_p} matrix matrix of dummy observation prior.
    Y_plus    = matrix(),

    #' @field X_plus a \code{KxT_p} matrix matrix of dummy observation prior.
    X_plus    = matrix(),

    #' @description
    #' Create a new prior specification \code{PriorBVARW16}.
    #' 
    #' @param data a \code{(T+p)xN} matrix with time series data.
    #' @param p a positive integer - the autoregressive lag order of the VAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @param stationary an \code{N} logical vector - its element set to
    #' \code{FALSE} sets the prior mean for the autoregressive parameters of the
    #' \code{N}th equation to the white noise process, otherwise to random walk.
    #' @param ar_sigma2 a positive \code{N}-vector with the autoregressive variance
    #' estimates for each variable to be used in the Minnesota prior for the autoregressive
    #' parameters.
    #' @param lambda a positive \code{4}-vector with the hyper-parameters of
    #' the Minnesota prior for the autoregressive parameters - the first element
    #' is the overall tightness hyper-parameter, while the second element is the
    #' tightness of the prior on the constant and exogenous variable coefficients.
    #' 
    #' @return A new prior specification \code{PriorBVARW16}.
    #' 
    #' @examples
    #' prior = specify_prior_bvarW16$new(us_macro_chan, p = 4)
    #' prior$A
    #'
    initialize = function(
      data,
      p,
      d = 0,
      stationary = rep(FALSE, ncol(data)),
      lambda = c(0.2^2, 10^2, 1, 1)
    ){
      N             = ncol(data)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      stopifnot("Argument stationary must be a logical vector of length N." = length(stationary) == N & is.logical(stationary))

      ar_sigma2     = apply(data, 2, function(x){sum(ar(x, aic = FALSE, order.max = 4)$resid^2, na.rm = TRUE) / (dim(data)[1] - 5)})

      K             = N * p + 1 + d
      self$lambda   = lambda
      self$A        = cbind(diag(as.numeric(!stationary)), matrix(0, N, K - N))
      self$V        = diag(c(self$lambda[1] / (kronecker(rep(1, p), ar_sigma2) * kronecker((1:p)^2, rep(1, N))), rep(self$lambda[2], 1 + d)))
      self$V_inv    = diag(1/diag(self$V))
      self$S        = diag(N)
      self$nu       = N + 3       # as in Chan (2020, JBES)

      ybar    = colMeans(matrix(data[1:p,], ncol = N))
      Ysoc    = diag(ybar)
      Ysur    = t(ybar)
      Xsoc    = cbind(kronecker(t(rep(1, p)), Ysoc), matrix(0, N, d + 1))
      Xsur    = cbind(kronecker(t(rep(1, p)), Ysur), 1, matrix(0, 1, d))

      self$Y_plus = t(rbind(Ysoc, Ysur))
      self$X_plus = t(rbind(Xsoc, Xsur))

    }, # END initialize

    #' @description
    #' Returns the elements of the prior specification \code{PriorBVARW16} as
    #' a \code{list}.
    #'
    #' @examples
    #' # a prior for 3-variable example with four lags
    #' prior = specify_prior_bvar$new(N = 3, p = 4)
    #' prior$get_prior() # show the prior as list
    #'
    get_prior = function(){
      list(
        A        = self$A,
        V        = self$V,
        V_inv    = self$V_inv,
        S        = self$S,
        nu       = self$nu,
        lambda   = self$lambda,
        Y_plus   = self$Y_plus,
        X_plus   = self$X_plus
      )
    } # END get_prior

  ) # END public
) # END specify_prior_bvar






#' R6 Class Representing \code{JointPriorBVARW16}
#'
#' @description
#' The class \code{JointPriorBVARW16} presents the joint prior for the BVAR model.
#'
#' @export
specify_joint_prior_bvarW16 = R6::R6Class(
  "JointPriorBVARW16",
  
  public = list(
    
    #' @field A a real-valued \code{NxK} matrix, the mean of
    #' the matrix-variate normal posterior distribution for the parameter
    #' matrix \eqn{A}.
    A          = matrix(),
    
    #' @field S a \code{NxN} positive definite scale matrix of the
    #' Inverse Wishart posterior distribution for the error term covariance
    #' matrix \eqn{\Sigma}.
    S    = matrix(),
    
    #' @field V a \code{KxK} positive definite row-specific covariance matrix of the
    #' Matrix-variate Normal posterior distribution for the autoregressive matrix
    #' \eqn{A}.
    V    = matrix(),
    
    #' @field V_inv a \code{KxK} positive definite row-specific precision matrix of the
    #' Matrix-variate Normal posterior distribution for the autoregressive matrix
    #' \eqn{A}.
    V_inv    = matrix(),
    
    #' @field nu a positive scalar, shape parameter of the Inverse
    #' Wishart posterior distribution for the error term covariance
    #' matrix \eqn{\Sigma}.
    nu   = numeric(),
    
    #' @description
    #' Compute the parameters of the analytical posterior distribution and
    #' return in an object of class  \code{JointPriorBVARW16}.
    #' @param data_matrices an object \code{DataMatricesBSVAR} with the data matrices.
    #' @param prior an object \code{PriorBVARW16} with the prior specification.
    #'
    #' @return Parameters of the analytical posterior distribution \code{JointPriorBVARW16}.
    #'
    initialize = function(
      data_matrices,
      prior
    ){
     
      post = bvar_posterior(prior$Y_plus, prior$X_plus, 
                            prior$A, prior$V, prior$V_inv, prior$S, prior$nu)
      
      self$V_inv  = post$V_inv
      self$V      = post$V
      self$A      = post$A
      self$S      = post$S
      self$nu     = post$nu
    }, # END initialize
    
    #' @description
    #' Returns the elements of the analytical posterior \code{JointPriorBVARW16} 
    #' as a \code{list}.
    #'
    get_joint_prior   = function(){
      list(
        A     = self$A,
        V     = self$V,
        V_inv = self$V_inv,
        S     = self$S,
        nu    = self$nu
      )
    } # END get_analytical_posterior
  ) # END public
) # END specify_joint_prior_bvarW16







#' R6 Class Representing \code{AnalyticalPosteriorBVARW16}
#'
#' @description
#' The class \code{AnalyticalPosteriorBVARW16} presents starting values for the BVAR model.
#'
#' @export
specify_analytical_posterior_bvarW16 = R6::R6Class(
  "AnalyticalPosteriorBVARW16",

  public = list(

    #' @field A a real-valued \code{NxK} matrix, the mean of
    #' the matrix-variate normal posterior distribution for the parameter
    #' matrix \eqn{A}.
    A          = matrix(),

    #' @field S a \code{NxN} positive definite scale matrix of the
    #' Inverse Wishart posterior distribution for the error term covariance
    #' matrix \eqn{\Sigma}.
    S    = matrix(),

    #' @field V a \code{KxK} positive definite row-specific covariance matrix of the
    #' Matrix-variate Normal posterior distribution for the autoregressive matrix
    #' \eqn{A}.
    V    = matrix(),

    #' @field V_inv a \code{KxK} positive definite row-specific precision matrix of the
    #' Matrix-variate Normal posterior distribution for the autoregressive matrix
    #' \eqn{A}.
    V_inv    = matrix(),

    #' @field nu a positive scalar, shape parameter of the Inverse
    #' Wishart posterior distribution for the error term covariance
    #' matrix \eqn{\Sigma}.
    nu   = numeric(),

    #' @description
    #' Compute the parameters of the analytical posterior distribution and
    #' return in an object of class  \code{AnalyticalPosteriorBVARW16}.
    #' @param data_matrices an object \code{DataMatricesBSVAR} with the data matrices.
    #' @param prior an object \code{PriorBVARW16} with the prior specification.
    #'
    #' @return Parameters of the analytical posterior distribution \code{AnalyticalPosteriorBVARW16}.
    #'
    initialize = function(
      data_matrices,
      prior
    ){

      post = bvar_posterior(data_matrices$Y, data_matrices$X, 
                     prior$A, prior$V, prior$V_inv, prior$S, prior$nu)
      
      self$V_inv  = post$V_inv
      self$V      = post$V
      self$A      = post$A
      self$S      = post$S
      self$nu     = post$nu
    }, # END initialize

    #' @description
    #' Returns the elements of the analytical posterior \code{AnalyticalPosteriorBVARW16} as a \code{list}.
    #'
    get_analytical_posterior   = function(){
      list(
        A     = self$A,
        V     = self$V,
        V_inv = self$V_inv,
        S     = self$S,
        nu    = self$nu
      )
    } # END get_analytical_posterior
  ) # END public
) # END specify_analytical_posterior_bvarW16



#' R6 Class representing the specification of the \code{BVARW16} model
#'
#' @description
#' The class \code{BVARW16} presents complete specification for the BVAR model.
#'
#' @examples
#' spec = specify_bvarW16$new(us_macro_chan)
#'
#' @export
specify_bvarW16 = R6::R6Class(
  "BVARW16",

  public = list(

    #' @field p a non-negative integer specifying the autoregressive lag order of the model.
    p                      = numeric(),

    #' @field prior an object \code{PriorBVARW16} with the prior specification.
    prior                  = list(),

    #' @field joint_prior an object \code{JointPriorBVARW16} with the prior specification.
    joint_prior           = list(),
    
    #' @field data_matrices an object \code{DataMatricesBSVAR} with the data matrices.
    data_matrices          = list(),

    #' @field analytical_posterior an object \code{AnalyticalPosteriorBVARW16} with the starting values.
    analytical_posterior        = list(),

    #' @description
    #' Create a new specification of the \code{BVARW16} model.
    #' 
    #' @param data a \code{(T+p)xN} matrix with time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @param exogenous a \code{(T+p)xd} matrix of exogenous variables.
    #' @param stationary an \code{N} logical vector - its element set to
    #' \code{FALSE} sets the prior mean for the autoregressive parameters of the
    #' \code{N}th equation to the white noise process, otherwise to random walk.
    #' @param lambda a positive \code{4}-vector with the hyper-parameters of
    #' the Minnesota prior for the autoregressive parameters - the first element
    #' is the overall tightness hyper-parameter, while the second element is the
    #' tightness of the prior on the constant and exogenous variable coefficients.
    #' 
    #' @return A new complete specification for the \code{BVARW16} model.
    initialize = function(
      data,
      p = 1L,
      exogenous = NULL,
      stationary = rep(FALSE, ncol(data)),
      lambda = c(0.2^2, 10^2, 1, 1)
    ) {

      stopifnot("Argument p has to be a positive integer." = ((p %% 1) == 0 & p > 0))
      self$p     = p

      TT            = nrow(data)
      T             = TT - self$p
      N             = ncol(data)
      d             = 0
      if (!is.null(exogenous)) {
        d           = ncol(exogenous)
      }
      K             = N * p + 1 + d

      self$data_matrices  = bsvars::specify_data_matrices$new(data, p, exogenous)
      self$prior          = specify_prior_bvarW16$new(data, p, d, stationary, lambda)
      dm = list(
        Y = self$prior$Y_plus,
        X = self$prior$X_plus
      )
      self$joint_prior          = specify_joint_prior_bvarW16$new(dm, self$prior)
      self$analytical_posterior = specify_analytical_posterior_bvarW16$new(self$data_matrices, self$joint_prior)
      
    }, # END initialize

    #' @description
    #' Returns the data matrices as the \code{DataMatricesBSVAR} object.
    #'
    #' @examples
    #' spec = specify_bvarW16$new(us_macro_chan)
    #' spec$get_data_matrices()
    #'
    get_data_matrices = function() {
      self$data_matrices$get_data_matrices()
    }, # END get_data_matrices

    #' @description
    #' Returns the prior specification as the \code{PriorBVARW16} object.
    #'
    #' @examples
    #' spec = specify_bvarW16$new(us_macro_chan)
    #' spec$get_prior()
    #'
    get_prior = function() {
      self$prior$clone()
    }, # END get_prior

    #' @description
    #' Returns the starting values as the \code{JointPriorBVARW16} object.
    #'
    #' @examples
    #' spec = specify_bvarW16$new(us_macro_chan)
    #' spec$get_joint_prior()
    #'
    get_joint_prior = function() {
      self$joint_prior$clone()
    }, # END get_joint_prior
    
    #' @description
    #' Returns the starting values as the \code{AnalyticalPosteriorBVARW16} object.
    #'
    #' @examples
    #' spec = specify_bvarW16$new(us_macro_chan)
    #' spec$get_analytical_posterior()
    #'
    get_analytical_posterior = function() {
      self$analytical_posterior$clone()
    } # END get_analytical_posterior
  ) # END public
) # END specify_bvar








#' R6 Class Representing \code{PosteriorBVARW16}
#'
#' @description
#' The class \code{PosteriorBVARW16} contains posterior output and the
#' specification including the specification for the BVAR model.
#'
#' @examples
#' spec = specify_bvarW16$new(us_macro_chan)
#' #post = estimate(spec, 5)
#'
#' @export
specify_posterior_bvarW16 = R6::R6Class(
  "PosteriorBVARW16",

  public = list(

    #' @field last_draw an object of class \code{BVARW16} with the specification
    #' of the BVAR model.
    last_draw = list(),

    #' @field posterior a list containing Bayesian estimation output collected
    #' in elements \code{A} and \code{Sigma}.
    posterior = list(),

    #' @description
    #' Create a new posterior output \code{PosteriorBVARW16}.
    #' @param specification_bvar an object of class \code{BVARW16} with the
    #' specification o fthe BVAR model.
    #' @param posterior_bvar a list containing Bayesian estimation output
    #' collected in elements \code{A} and \code{Sigma}.
    #' 
    #' @return A posterior output \code{PosteriorBVARW16}.
    initialize = function(specification_bvar, posterior_bvar) {

      stopifnot("Argument specification_bsvar must be of class BVARW16." = any(class(specification_bvar) == "BVARW16"))
      stopifnot("Argument posterior_bsvar must must contain MCMC output." = is.list(posterior_bvar) & is.array(posterior_bvar$A) & is.array(posterior_bvar$Sigma))

      self$last_draw    = specification_bvar
      self$posterior    = posterior_bvar
    }, # END initialize

    #' @description
    #' Returns a list containing Bayesian estimation output collected in elements
    #' \code{A}, \code{Sigma}.
    #'
    #' @examples
    #' spec = specify_bvarW16$new(us_macro_chan)
    #' #post = estimate(spec, 5)
    #' #post$get_posterior()
    #'
    get_posterior       = function(){
      self$posterior
    }, # END get_posterior

    #' @description
    #' Returns an object of class \code{BVARW16} with the specification of the 
    #' BVAR model.
    #'
    #' @examples
    #' spec = specify_bvarW16$new(us_macro_chan)
    #' #post = estimate(burn, 5)
    #' #post$get_last_draw()
    get_last_draw      = function(){
      self$last_draw$clone()
    } # END get_last_draw

  ) # END public
) # END specify_posterior_bsvarW16
