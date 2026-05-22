
#' R6 Class Representing \code{PriorBVAR}
#'
#' @description
#' The class \code{PriorBVAR} presents a prior specification for the BVAR model.
#' 
#' #' \strong{The Model.} 
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
#' @examples 
#' prior = specify_prior_bvar$new(N = 3, p = 1)  # a prior for 3-variable example with one lag
#' prior$A                                        # show autoregressive prior mean
#' 
#' @export
specify_prior_bvar = R6::R6Class(
  "PriorBVAR",
  
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
    
    #' @field sv_a a positive scalar, the shape parameter of the gamma prior in 
    #' the hierarchical prior for the common stochastic volatility. 
    sv_a     = numeric(),
    
    #' @field sv_s a positive scalar, the scale parameter of the gamma prior in 
    #' the hierarchical prior for the common stochastic volatility.
    sv_s     = numeric(),
    
    #' @description
    #' Create a new prior specification \code{PriorBVAR}.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the VAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @param stationary an \code{N} logical vector - its element set to 
    #' \code{FALSE} sets the prior mean for the autoregressive parameters of the 
    #' \code{N}th equation to the white noise process, otherwise to random walk.
    #' @param is_homoskedastic a logical scalar - if \code{TRUE} the model assumes 
    #' homoskedastic errors, otherwise it assumes stochastic volatility.
    #' @return A new prior specification \code{PriorBVAR}.
    #' @examples 
    #' # a prior for 3-variable example with one lag and stationary data
    #' prior = specify_prior_bvar$new(N = 3, p = 1, stationary = rep(TRUE, 3))
    #' prior$A # show autoregressive prior mean
    #' 
    initialize = function(N, p, d = 0, stationary = rep(FALSE, N), is_homoskedastic = TRUE){
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      stopifnot("Argument stationary must be a logical vector of length N." = length(stationary) == N & is.logical(stationary))
      
      K             = N * p + 1 + d
      self$A        = cbind(diag(as.numeric(!stationary)), matrix(0, N, K - N))
      self$S        = diag(N)
      self$nu       = N + 3       # as in Chan (2020, JBES)
      self$Psi      = diag(K)
      self$Gamma    = diag(K)
      self$lambda   = N + 1
      self$sv_a     = NA
      self$sv_s     = NA
      if (!is_homoskedastic) {
        self$sv_a   = 1
        self$sv_s   = 0.1
      }
    }, # END initialize
    
    #' @description
    #' Returns the elements of the prior specification \code{PriorBVAR} as 
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
        S        = self$S,
        nu       = self$nu,
        Psi      = self$Psi,
        Gamma    = self$Gamma,
        lambda   = self$lambda,
        sv_a     = self$sv_a,
        sv_s     = self$sv_s
      )
    } # END get_prior
    
  ) # END public
) # END specify_prior_bvar



#' R6 Class Representing \code{StartingValuesBVAR}
#'
#' @description
#' The class \code{StartingValuesBVAR} presents starting values for the BVAR model.
#' 
#' @examples 
#' # starting values for a 3-variable BVAR model.
#' sv = specify_starting_values_bvar$new(N = 3, p = 4, T = 100)
#' 
#' @export
specify_starting_values_bvar = R6::R6Class(
  "StartingValuesBVAR",
  
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
    
    #' @field h an \code{T}-vector with the starting values of the 
    #' log-volatility processes.
    h             = numeric(),
    
    #' @field rho a scalalr for the SV autoregressive parameter.
    rho           = numeric(),
    
    #' @field omega a scalar for the SV process conditional standard deviation.
    omega         = numeric(),
    
    #' @field sigma2v a scalar for SV process conditional variances.
    sigma2v       = numeric(),
    
    #' @field S a \code{T} integer vector with the auxiliary mixture 
    #' component indicator.
    S             = numeric(),
    
    #' @field sigma2_omega a scalar for the variance of the zero-mean 
    #' normal prior for \eqn{\omega}.
    sigma2_omega  = numeric(),
    
    #' @field s_ a positive scalar with the scale of the gamma prior of the 
    #' hierarchical prior for \eqn{\sigma^2_{\omega}}.
    s_            = numeric(),
    
    #' @field lambda a \code{T}-vetor of starting values for latent variable.
    lambda        = numeric(),
    
    #' @field df a scalar greater than 2 with the starting value 
    #' for the degrees of freedom parameter of the Student-t 
    #' conditional distribution of error term.
    df            = numeric(),
    
    #' @description
    #' Create new starting values \code{StartingValuesBVAR}.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the BVAR model.
    #' @param T a positive integer - the number of time periods in the data.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @param is_homoskedastic a logical scalar - if \code{TRUE} the model assumes 
    #' homoskedastic errors, otherwise it assumes stochastic volatility.
    #' @param is_normal a logical scalar - if \code{TRUE} the model assumes normal 
    #' error term, otherwise, it assumes Student-t errors.
    #' @param ar_sigma2 a positive \code{N}-vector with the autoregressive variance
    #' estimates for each variable to be used in the Minnesota prior for the autoregressive
    #' parameters.
    #' @param kappa a positive \code{2}-vector with the hyperparameters of
    #' the Minnesota prior for the autoregressive parameters - the first element 
    #' is the overall tightness hyperparameter, while the second element is the 
    #' tightness of the prior on the constant and exogenous variable coefficients.
    #' @return Starting values \code{StartingValuesBVAR}.
    #' @examples 
    #' # starting values for a 3-variable BVAR model
    #' sv = specify_starting_values_bvar$new(N = 3, p = 4, T = 100)
    #' 
    initialize = function(
      N, 
      p, 
      T,
      d = 0, 
      is_homoskedastic = TRUE, 
      is_normal = TRUE,
      ar_sigma2 = rep(1, N),
      kappa = c(0.2^2, 10^2)
    ){
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument T must be a positive integer number." = T > 0 & T %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      
      K                   = N * p + 1 + d
      self$A              = cbind(diag(runif(N)), matrix(0, N, K - N))
      self$Sigma          = diag(rgamma(N, 1))
      self$V              = diag(c(kappa[1] / (kronecker(rep(1, p), ar_sigma2) * kronecker((1:p)^2, rep(1, N))), rep(kappa[2], 1 + d)))
      
      self$h              = rnorm(T, sd = .01)
      self$rho            = .5
      self$omega          = .1
      self$sigma2v        = .1^2
      self$S              = rep(1, T)
      self$sigma2_omega   = 1
      self$s_             = 0.05
      
      self$lambda         = rep(1, T)
      self$df             = 30
    }, # END initialize
    
    #' @description
    #' Returns the elements of the starting values \code{StartingValuesBVAR} as a \code{list}.
    #' 
    #' @examples 
    #' # starting values for a 3-variable BVAR model
    #' sv = specify_starting_values_bvar$new(N = 3, p = 4, T = 100)
    #' sv$get_starting_values()   # show starting values as list
    #' 
    get_starting_values   = function(){
      list(
        A                 = self$A,
        Sigma             = self$Sigma,
        V                 = self$V,
        h                 = self$h,
        rho               = self$rho,
        omega             = self$omega,
        sigma2v           = self$sigma2v,
        S                 = self$S,
        sigma2_omega      = self$sigma2_omega,
        s_                = self$s_,
        lambda            = self$lambda,
        df                = self$df
      )
    }, # END get_starting_values
    
    #' @description
    #' Sets the elements of the starting values \code{StartingValuesBVAR} to 
    #' provided values.
    #' @param last_draw a list containing the last draw of elements \code{A} - 
    #' a \code{KxN} matrix, \code{Sigma} - an \code{NxN} matrix, and \code{V} - 
    #' a \code{KxK} matrix.
    #' @return An object of class \code{StartingValuesBVAR} including the 
    #' last draw of the current MCMC as the starting value to be passed to the 
    #' continuation of the MCMC estimation using \code{estimate()}.
    #' 
    #' @examples 
    #' # starting values for a 3-variable BVAR model
    #' sv = specify_starting_values_bvar$new(N = 3, p = 4, T = 100)
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
      self$h            = last_draw$h
      self$rho          = last_draw$rho
      self$omega        = last_draw$omega
      self$sigma2v      = last_draw$sigma2v
      self$S            = last_draw$S
      self$sigma2_omega = last_draw$sigma2_omega
      self$s_           = last_draw$s_
      self$lambda       = last_draw$lambda
      self$df           = last_draw$df
    } # END set_starting_values
  ) # END public
) # END specify_starting_values_bvar



#' R6 Class representing the specification of the \code{BVAR} model
#'
#' @description
#' The class \code{BVAR} presents complete specification for the BVAR model.
#' 
#' @examples 
#' spec = specify_bvar$new(us_macro_chan)
#' 
#' @export
specify_bvar = R6::R6Class(
  "BVAR",
  
  private = list(
    normal        = TRUE,
    homoskedastic = TRUE,
    centred_sv    = FALSE
  ), # END private
  
  public = list(
    
    #' @field p a non-negative integer specifying the autoregressive lag order of the model. 
    p                      = numeric(),
    
    #' @field prior an object \code{PriorBVAR} with the prior specification. 
    prior                  = list(),
    
    #' @field data_matrices an object \code{DataMatricesBSVAR} with the data matrices.
    data_matrices          = list(),
    
    #' @field starting_values an object \code{StartingValuesBVAR} with the starting values.
    starting_values        = list(),
    
    #' @description
    #' Create a new specification of the \code{BVAR} model.
    #' @param data a \code{(T+p)xN} matrix with time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @param exogenous a \code{(T+p)xd} matrix of exogenous variables. 
    #' @param common_volatility a character string specifying the common volatility 
    #' component of the error term covariance matrix. It can take three values: 
    #' \code{homoskedastic} - the model assumes homoskedastic errors, 
    #' \code{ncSV} - the model assumes non-centred stochastic volatility, and 
    #' \code{cSV} - the model assumes centred stochastic volatility.
    #' @param distribution a character string specifying the conditional distribution 
    #' of structural shocks. Value \code{"norm"} sets it to the normal distribution, 
    #' while value \code{"t"} sets the Student-t distribution.
    #' @param stationary an \code{N} logical vector - its element set to
    #' \code{FALSE} sets the prior mean for the autoregressive parameters of the 
    #' \code{N}th equation to the white noise process, otherwise to random walk.
    #' @return A new complete specification for the \code{BVAR} model.
    initialize = function(
    data,
    p = 1L,
    exogenous = NULL,
    common_volatility = c("homoskedastic", "ncSV", "cSV"),
    distribution = c("norm","t"),
    stationary = rep(FALSE, ncol(data))
    ) {
      stopifnot("Argument p has to be a positive integer." = ((p %% 1) == 0 & p > 0))
      self$p     = p
      
      common_volatility = match.arg(common_volatility)
      private$homoskedastic  = common_volatility == "homoskedastic"
      if (common_volatility == "cSV") {
        private$centred_sv = TRUE
      }
      
      distribution      = match.arg(distribution)
      private$normal    = distribution == "norm"
      TT            = nrow(data)
      T             = TT - self$p
      N             = ncol(data)
      d             = 0
      if (!is.null(exogenous)) {
        d           = ncol(exogenous)
      }
      K             = N * p + 1 + d
      
      ar_sigma2     = apply(data, 2, function(x){sum(ar(x, aic = FALSE, order.max = 4)$resid^2, na.rm=TRUE) / (dim(data)[1] - 5)})
      
      self$data_matrices   = bsvars::specify_data_matrices$new(data, p, exogenous)
      self$prior           = specify_prior_bvar$new(N, p, d, stationary, private$homoskedastic)
      self$starting_values = specify_starting_values_bvar$new(N, self$p, T, d, private$homoskedastic, private$normal, ar_sigma2)
      
    }, # END initialize
    
    #' @description
    #' Returns the logical value of whether the conditional shock distribution is normal.
    #' 
    #' @examples 
    #' spec = specify_bvar$new(us_macro_chan)
    #' spec$get_normal()
    #' 
    get_normal = function() {
      private$normal
    }, # END get_normal
    
    #' @description
    #' Returns the logical value of whether the common volatility is homoskedastic.
    #' 
    #' @examples 
    #' spec = specify_bvar$new(us_macro_chan)
    #' spec$get_homoskedastic()
    #' 
    get_homoskedastic = function() {
      private$homoskedastic
    }, # END get_homoskedastic
    
    #' @description
    #' Returns the logical value of whether the common volatility is centred 
    #' Stochastic Volatility
    #' 
    #' @examples 
    #' spec = specify_bvar$new(us_macro_chan)
    #' spec$get_centred_sv()
    #' 
    get_centred_sv = function() {
      private$centred_sv
    }, # END get_centred_sv
    
    #' @description
    #' Returns the data matrices as the \code{DataMatricesBSVAR} object.
    #' 
    #' @examples 
    #' spec = specify_bvar$new(us_macro_chan)
    #' spec$get_data_matrices()
    #' 
    get_data_matrices = function() {
      self$data_matrices$get_data_matrices()
    }, # END get_data_matrices
    
    #' @description
    #' Returns the prior specification as the \code{PriorBVAR} object.
    #' 
    #' @examples 
    #' spec = specify_bvar$new(us_macro_chan)
    #' spec$get_prior()
    #' 
    get_prior = function() {
      self$prior$clone()
    }, # END get_prior
    
    #' @description
    #' Returns the starting values as the \code{StartingValuesBVAR} object.
    #' 
    #' @examples 
    #' spec = specify_bvar$new(us_macro_chan)
    #' spec$get_starting_values()
    #' 
    get_starting_values = function() {
      self$starting_values$clone()
    } # END get_starting_values
  ) # END public
) # END specify_bvar








#' R6 Class Representing \code{PosteriorBVAR}
#'
#' @description
#' The class \code{PosteriorBVAR} contains posterior output and the 
#' specification including the last MCMC draw for the BVAR model. 
#' Note that due to the thinning of the MCMC output the starting value in element 
#' \code{last_draw} might not be equal to the last draw provided in element 
#' \code{posterior}.
#' 
#' @examples 
#' # This is a function that is used within estimate()
#' spec = specify_bvar$new(us_macro_chan)
#' post = estimate(spec, 5)
#' class(post)
#' 
#' @export
specify_posterior_bvar = R6::R6Class(
  "PosteriorBVAR",
  
  public = list(
    
    #' @field last_draw an object of class \code{BVAR} with the last draw of 
    #' the current MCMC run as the starting value to be passed to the 
    #' continuation of the MCMC estimation using \code{estimate()}. 
    last_draw = list(),
    
    #' @field posterior a list containing Bayesian estimation output collected 
    #' in elements \code{A},  \code{Sigma}, and \code{V}.
    posterior = list(),
    
    #' @description
    #' Create a new posterior output \code{PosteriorBVAR}.
    #' @param specification_bvar an object of class \code{BVAR} with the 
    #' last draw of the current MCMC run as the starting value.
    #' @param posterior_bvar a list containing Bayesian estimation output 
    #' collected in elements \code{A}, \code{Sigma}, and \code{V}.
    #' @return A posterior output \code{PosteriorBVAR}.
    initialize = function(specification_bvar, posterior_bvar) {
      
      stopifnot("Argument specification_bsvar must be of class BVAR." = any(class(specification_bvar) == "BVAR"))
      stopifnot("Argument posterior_bsvar must must contain MCMC output." = is.list(posterior_bvar) & is.array(posterior_bvar$A) & is.array(posterior_bvar$Sigma) & is.array(posterior_bvar$V))
      
      self$last_draw    = specification_bvar
      self$posterior    = posterior_bvar
    }, # END initialize
    
    #' @description
    #' Returns a list containing Bayesian estimation output collected in elements 
    #' \code{A}, \code{Sigma}, and \code{V}.
    #' 
    #' @examples 
    #' spec = specify_bvar$new(us_macro_chan)
    #' post = estimate(spec, 5)
    #' post$get_posterior()
    #' 
    get_posterior       = function(){
      self$posterior
    }, # END get_posterior
    
    #' @description
    #' Returns an object of class \code{BVAR} with the last draw of the 
    #' current MCMC run as the starting value to be passed to the continuation 
    #' of the MCMC estimation using \code{estimate()}.
    #' 
    #' @examples
    #' spec = specify_bvar$new(us_macro_chan)
    #' burn = estimate(spec, 5)
    #' post = estimate(burn, 5)
    get_last_draw      = function(){
      self$last_draw$clone()
    } # END get_last_draw
    
  ) # END public
) # END specify_posterior_bsvar
