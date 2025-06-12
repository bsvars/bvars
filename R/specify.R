
#' R6 Class Representing \code{PriorBVARGIG}
#'
#' @description
#' The class \code{PriorBVARGIG} presents a prior specification for the BVAR model.
#' 
#' @examples 
#' prior = specify_prior_bsvar$new(N = 3, p = 1)  # a prior for 3-variable example with one lag
#' prior$A                                        # show autoregressive prior mean
#' 
#' @export
specify_prior_bvarGIG = R6::R6Class(
  "PriorBVARGIG",
  
  public = list(
    
    #' @field A a real-valued \code{NxK} matrix, the mean matrix \eqn{A_0} of 
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
    #' generalized inverse Gaussian distribution for the equation-specific prior 
    #' covariance \eqn{V}
    Psi   = matrix(),
    
    #' @field Gamma a \code{KxK} scale matrix \eqn{\Gamma_0} of the matrix 
    #' generalized inverse Gaussian distribution for the equation-specific prior 
    #' covariance \eqn{V}
    Gamma   = matrix(),
    
    #' @field lambda a positive scalar shape parameter \eqn{\lambda_0} of the 
    #' matrix generalized inverse Gaussian distribution for the equation-specific 
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



#' R6 Class Representing \code{StartingValuesBVARGIG}
#'
#' @description
#' The class \code{StartingValuesBVARGIG} presents starting values for the BVAR model.
#' 
#' @examples 
#' # starting values for a 3-variable BVAR model.
#' sv = specify_starting_values_bvarGIG$new(N = 3, p = 1)
#' 
#' @export
specify_starting_values_bvarGIG = R6::R6Class(
  "StartingValuesBVARGIG",
  
  public = list(
    
    #' @field A an \code{NxK} matrix of starting values for the autoregressive 
    #' matrix \eqn{A}. 
    A             = matrix(),
    
    #' @field Sigma an \code{NxN} matrix of starting values for the error term 
    #' covariance \eqn{\Sigma}. 
    Sigma             = matrix(),
    
    #' @field V a \code{KxK} matrix of starting values for the prior 
    #' equation-specific covariance \eqn{V} of the hierarchical prior distribution
    #' for matrix \eqn{A}. 
    V         = matrix(),
    
    #' @description
    #' Create new starting values \code{StartingValuesBVARGIG}.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the BVAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @return Starting values \code{StartingValuesBVARGIG}.
    #' @examples 
    #' # starting values for a 3-variable BVAR model
    #' sv = specify_starting_values_bvarGIG$new(N = 3, p = 4)
    #' 
    initialize = function(N, p, d = 0){
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      
      K                   = N * p + 1 + d
      self$A              = cbind(diag(runif(N)), matrix(0, N, K - N))
      self$Sigma          = diag(rgamma(N, 1))
      self$V              = diag(rgamma(K, 1))
    }, # END initialize
    
    #' @description
    #' Returns the elements of the starting values \code{StartingValuesBVARGIG} as a \code{list}.
    #' 
    #' @examples 
    #' # starting values for a 3-variable BVAR model
    #' sv = specify_starting_values_bvarGIG$new(N = 3, p = 1)
    #' sv$get_starting_values()   # show starting values as list
    #' 
    get_starting_values   = function(){
      list(
        A                 = self$A,
        Sigma             = self$Sigma,
        V                 = self$V
      )
    }, # END get_starting_values
    
    #' @description
    #' Sets the elements of the starting values \code{StartingValuesBVARGIG} to 
    #' provided values.
    #' @param last_draw a list containing the last draw of elements \code{A} - 
    #' a \code{KxN} matrix, \code{Sigma} - an \code{NxN} matrix, and \code{V} - 
    #' a \code{KxK} matrix.
    #' @return An object of class \code{StartingValuesBVARGIG} including the 
    #' last draw of the current MCMC as the starting value to be passed to the 
    #' continuation of the MCMC estimation using \code{estimate()}.
    #' 
    #' @examples 
    #' # starting values for a 3-variable BVAR model
    #' sv = specify_starting_values_bvarGIG$new(N = 3, p = 1)
    #' 
    #' # Modify the starting values by:
    #' sv_list = sv$get_starting_values()   # getting them as list
    #' sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
    #' sv$set_starting_values(sv_list)      # providing to the class object
    #' 
    set_starting_values   = function(last_draw) {
      self$A            = last_draw$A
      self$Sigma        = last_draw$Sigma
      self$V            = last_draw$V
    } # END set_starting_values
  ) # END public
) # END specify_starting_values_bvarGIG



#' R6 Class representing the specification of the \code{BVARGIG} model
#'
#' @description
#' The class \code{BVARGIG} presents complete specification for the BVAR model.
#' 
#' @examples 
#' data(us_fiscal_lsuw)
#' spec = specify_bvarGIG$new(
#'    data = us_fiscal_lsuw,
#'    p = 4
#' )
#' 
#' @export
specify_bvarGIG = R6::R6Class(
  "BVARGIG",
  
  public = list(
    
    #' @field p a non-negative integer specifying the autoregressive lag order of the model. 
    p                      = numeric(),
    
    #' @field prior an object \code{PriorBVARGIG} with the prior specification. 
    prior                  = list(),
    
    #' @field data_matrices an object \code{DataMatricesBSVAR} with the data matrices.
    data_matrices          = list(),
    
    #' @field starting_values an object \code{StartingValuesBVARGIG} with the starting values.
    starting_values        = list(),
    
    #' @description
    #' Create a new specification of the \code{BVARGIG} model.
    #' @param data a \code{(T+p)xN} matrix with time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @param exogenous a \code{(T+p)xd} matrix of exogenous variables. 
    #' @param stationary an \code{N} logical vector - its element set to
    #' \code{FALSE} sets the prior mean for the autoregressive parameters of the 
    #' \code{N}th equation to the white noise process, otherwise to random walk.
    #' @return A new complete specification for the \code{BVARGIG} model.
    initialize = function(
    data,
    p = 1L,
    exogenous = NULL,
    stationary = rep(FALSE, ncol(data))
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
      
      self$data_matrices   = specify_data_matrices$new(data, p, exogenous)
      self$prior           = specify_prior_bvarGIG$new(N, p, d, stationary)
      self$starting_values = specify_starting_values_bvarGIG$new(N, self$p, d)
    }, # END initialize
    
    #' @description
    #' Returns the data matrices as the \code{DataMatricesBSVAR} object.
    #' 
    #' @examples 
    #' data(us_fiscal_lsuw)
    #' spec = specify_bvarGIG$new(
    #'    data = us_fiscal_lsuw,
    #'    p = 4
    #' )
    #' spec$get_data_matrices()
    #' 
    get_data_matrices = function() {
      self$data_matrices$clone()
    }, # END get_data_matrices
    
    #' @description
    #' Returns the prior specification as the \code{PriorBVARGIG} object.
    #' 
    #' @examples 
    #' data(us_fiscal_lsuw)
    #' spec = specify_bvarGIG$new(
    #'    data = us_fiscal_lsuw,
    #'    p = 4
    #' )
    #' spec$get_prior()
    #' 
    get_prior = function() {
      self$prior$clone()
    }, # END get_prior
    
    #' @description
    #' Returns the starting values as the \code{StartingValuesBVARGIG} object.
    #' 
    #' @examples 
    #' data(us_fiscal_lsuw)
    #' spec = specify_bvarGIG$new(
    #'    data = us_fiscal_lsuw,
    #'    p = 4
    #' )
    #' spec$get_starting_values()
    #' 
    get_starting_values = function() {
      self$starting_values$clone()
    } # END get_starting_values
  ) # END public
) # END specify_bvarGIG

