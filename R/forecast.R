
#' @export
generics::forecast

#' @title Forecasting using Structural Vector Autoregression
#'
#' @description Samples from the joint predictive density of all of the dependent 
#' variables for the model by Chan (2020) <doi:10.1080/07350015.2018.1451336>, that is, a 
#' Bayesian Vector Autoregression with Minnesota priors and a flexible structure 
#' of the error term specification. The latter includes: conditional multivariate 
#' normal or Student’s t distributions, as well as homoskedastic or heteroskedastic 
#' specifications with a common volatility modelled by centred or non-centred 
#' Stochastic Volatility.
#' 
#' @method forecast PosteriorBVAR
#' 
#' @param object posterior estimation outcome - an object of class 
#' \code{PosteriorBVAR} obtained by running the \code{estimate} function.
#' @param horizon a positive integer, specifying the forecasting horizon.
#' @param exogenous_forecast a matrix of dimension \code{horizon x d} containing 
#' forecasted values of the exogenous variables. 
#' @param conditional_forecast a \code{horizon x N} matrix with forecasted values 
#' for selected variables. It should only contain \code{numeric} or \code{NA} 
#' values. The entries with \code{NA} values correspond to the values that are 
#' forecasted conditionally on the realisations provided as \code{numeric} values.
#' @param ... not used
#' 
#' @return A list of class \code{Forecasts} containing the
#' draws from the predictive density and data. The output list includes element:
#' 
#' \describe{
#'  \item{forecasts}{an \code{NxTxS} array with the draws from predictive density}
#'  \item{forecast_mean}{an \code{NxhorizonxS} array with the mean of the predictive density}
#'  \item{forecast_covariance}{an \code{NxNxhorizonxS} array with the covariance of the predictive density}
#'  \item{Y}{an \eqn{NxT} matrix with the data on dependent variables}
#' }
#' 
#' @author Rui Liu \email{rl3023@columbia.edu}, 
#' Andres Ramirez Hassan \email{aramir21@gmail.com} & 
#' Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' spec = specify_bvar$new(us_macro_chan)   # specify the model
#' burn = estimate(spec, 5)                      # run the burn-in
#' post = estimate(burn, 5)                     # estimate the model
#' pred = forecast(post, 4)                      # forecast 1 year ahead
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' us_macro_chan |>
#'   specify_bvar$new() |>
#'   estimate(S = 5) |> 
#'   estimate(S = 5) |> 
#'   forecast(horizon = 4) -> pred
#' 
#' @export
forecast.PosteriorBVAR = function(
    object, 
    horizon = 1, 
    exogenous_forecast = NULL,
    conditional_forecast = NULL,
    ...
) {
  
  posterior_Sigma = object$posterior$Sigma
  posterior_A     = object$posterior$A
  T               = ncol(object$last_draw$data_matrices$X)
  X_T             = object$last_draw$data_matrices$X[,T]
  Y               = object$last_draw$data_matrices$Y
  homoskedastic   = object$last_draw$get_homoskedastic()
  normal          = object$last_draw$get_normal()
  
  N               = dim(posterior_Sigma)[1]
  K               = length(X_T)
  d               = K - N * object$last_draw$p - 1
  S               = dim(posterior_Sigma)[3]
  
  # prepare forecasting with exogenous variables
  if (d == 0 ) {
    exogenous_forecast = matrix(NA, horizon, 1)
  } else {
    stopifnot("Forecasted values of exogenous variables are missing." = (d > 0) & !is.null(exogenous_forecast))
    stopifnot("The matrix of exogenous_forecast does not have a correct number of columns." = ncol(exogenous_forecast) == d)
    stopifnot("Provide exogenous_forecast for all forecast periods specified by argument horizon." = nrow(exogenous_forecast) == horizon)
    stopifnot("Argument exogenous has to be a matrix." = is.matrix(exogenous_forecast) & is.numeric(exogenous_forecast))
    stopifnot("Argument exogenous cannot include missing values." = sum(is.na(exogenous_forecast)) == 0 )
  }
  
  # prepare forecasting with conditional forecasts
  if ( is.null(conditional_forecast) ) {
    # this will not be used for forecasting, but needs to be provided
    conditional_forecast = matrix(NA, horizon, N)
  } else {
    stopifnot("Argument conditional_forecast must be a matrix with numeric values."
              = is.matrix(conditional_forecast) & is.numeric(conditional_forecast)
    )
    stopifnot("Argument conditional_forecast must have the number of rows equal to 
              the value of argument horizon."
              = nrow(conditional_forecast) == horizon
    )
    stopifnot("Argument conditional_forecast must have the number of columns 
              equal to the number of columns in the used data."
              = ncol(conditional_forecast) == N
    )
  }
  
  # forecast volatility
  forecast_sigma2   = matrix(1, horizon, S)
  if (!homoskedastic) {
    posterior_h_T   = object$posterior$h[T,]
    posterior_rho   = object$posterior$rho
    posterior_omega = object$posterior$omega
    
    forecast_sigma2 = .Call(`_bvars_forecast_sigma2_sv1`, 
                            posterior_h_T, posterior_rho, posterior_omega, horizon
                      ) # END .Call
  }
  
  # forecast Student-t
  forecast_lambda   = matrix(1, horizon, S)
  if (!normal) {
    posterior_df    = object$posterior$df
    forecast_lambda = .Call(`_bvars_forecast_lambda_t1`, 
                            posterior_df, horizon
                      ) # END .Call
  }
  forecast_sigma2 = forecast_sigma2 * forecast_lambda
      
  # perform forecasting
  fore        = .Call(`_bvars_forecast_bvarGIG`, 
                      posterior_Sigma,
                      posterior_A,
                      forecast_sigma2,    # (horizon, S)
                      X_T,
                      exogenous_forecast,
                      conditional_forecast,
                      horizon
  ) # END .Call
  
  fore = specify_forecasts$new(fore, Y)
  return(fore)
} # END forecast.PosteriorBVARGIG
