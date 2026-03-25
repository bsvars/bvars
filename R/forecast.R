
#' @title Forecasting using Structural Vector Autoregression
#'
#' @description Samples from the joint predictive density of all of the dependent 
#' variables for models from packages \pkg{bsvars}, \pkg{bsvarSIGNs} or 
#' \pkg{bvarPANELs} at forecast horizons from 1 to \code{horizon} specified as 
#' an argument of the function.
#' 
#' @method forecast PosteriorBVARGIG
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' \code{PosteriorBVARGIG} obtained by running the \code{estimate} function.
#' @param horizon a positive integer, specifying the forecasting horizon.
#' @param exogenous_forecast a matrix of dimension \code{horizon x d} containing 
#' forecasted values of the exogenous variables. 
#' @param conditional_forecast a \code{horizon x N} matrix with forecasted values 
#' for selected variables. It should only contain \code{numeric} or \code{NA} 
#' values. The entries with \code{NA} values correspond to the values that are 
#' forecasted conditionally on the realisations provided as \code{numeric} values.
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
#' # upload data
#' data(us_fiscal_lsuw)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' specification  = specify_bvarGIG$new(us_fiscal_lsuw, p = 1)
#' 
#' # run the burn-in
#' burn_in        = estimate(specification, 5)
#' 
#' # estimate the model
#' posterior      = estimate(burn_in, 10)
#' 
#' # sample from predictive density 1 year ahead
#' predictive     = forecast(posterior, 4)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' us_fiscal_lsuw |>
#'   specify_bvarGIG$new(p = 1) |>
#'   estimate(S = 5) |> 
#'   estimate(S = 10) |> 
#'   forecast(horizon = 4) -> predictive
#' 
#' @export
forecast.PosteriorBVARGIG = function(
    posterior, 
    horizon = 1, 
    exogenous_forecast = NULL,
    conditional_forecast = NULL
) {
  
  posterior_Sigma = posterior$posterior$Sigma
  posterior_A     = posterior$posterior$A
  T               = ncol(posterior$last_draw$data_matrices$X)
  X_T             = posterior$last_draw$data_matrices$X[,T]
  Y               = posterior$last_draw$data_matrices$Y
  
  N               = dim(posterior_Sigma)[1]
  K               = length(X_T)
  d               = K - N * posterior$last_draw$p - 1
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
  forecast_sigma2   = array(1, c(N, horizon, S))
  
  # perform forecasting
  fore        = .Call(`_bvars_forecast_bvarGIG`, 
                      posterior_Sigma,
                      posterior_A,
                      forecast_sigma2,    # (N, horizon, S)
                      X_T,
                      exogenous_forecast,
                      conditional_forecast,
                      horizon
  ) # END .Call
  
  SS                  = dim(fore$forecasts)[3]
  forecast_covariance = array(NA, c(N, N, horizon, SS))
  for (s in 1:SS) forecast_covariance[,,,s] = fore$forecast_cov[s,][[1]]
  fore$forecast_covariance = forecast_covariance
  
  fore$Y          = Y
  class(fore)     = "Forecasts"
  
  return(fore)
} # END forecast.PosteriorBVARGIG
