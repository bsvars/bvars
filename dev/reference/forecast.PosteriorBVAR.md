# Forecasting using Structural Vector Autoregression

Samples from the joint predictive density of all of the dependent
variables for the model by Chan (2020)
\<doi:10.1080/07350015.2018.1451336\>, that is, a Bayesian Vector
Autoregression with Minnesota priors and a flexible structure of the
error term specification. The latter includes: conditional multivariate
normal or Student’s t distributions, as well as homoskedastic or
heteroskedastic specifications with a common volatility modelled by
centred or non-centred Stochastic Volatility.

## Usage

``` r
# S3 method for class 'PosteriorBVAR'
forecast(
  object,
  horizon = 1,
  exogenous_forecast = NULL,
  conditional_forecast = NULL,
  ...
)
```

## Arguments

- object:

  posterior estimation outcome - an object of class `PosteriorBVAR`
  obtained by running the `estimate` function.

- horizon:

  a positive integer, specifying the forecasting horizon.

- exogenous_forecast:

  a matrix of dimension `horizon x d` containing forecasted values of
  the exogenous variables.

- conditional_forecast:

  a `horizon x N` matrix with forecasted values for selected variables.
  It should only contain `numeric` or `NA` values. The entries with `NA`
  values correspond to the values that are forecasted conditionally on
  the realisations provided as `numeric` values.

- ...:

  not used

## Value

A list of class `Forecasts` containing the draws from the predictive
density and data. The output list includes element:

- forecasts:

  an `NxTxS` array with the draws from predictive density

- forecast_mean:

  an `NxhorizonxS` array with the mean of the predictive density

- forecast_covariance:

  an `NxNxhorizonxS` array with the covariance of the predictive density

- Y:

  an \\NxT\\ matrix with the data on dependent variables

## Author

Rui Liu <rl3023@columbia.edu>, Andres Ramirez Hassan
<aramir21@gmail.com> & Tomasz Woźniak <wozniak.tom@pm.me>

## Examples

``` r
spec = specify_bvar$new(us_macro_chan)   # specify the model
burn = estimate(spec, 5)                      # run the burn-in
#> **************************************************|
#> bvars: Forecasting with Large                     |
#>        Bayesian Vector Autoregressions            |
#> **************************************************|
#>  Gibbs sampler for the BVAR model                 |
#> **************************************************|
#>  Progress of the MCMC simulation for 5 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|
post = estimate(burn, 5)                     # estimate the model
#> **************************************************|
#> bvars: Forecasting with Large                     |
#>        Bayesian Vector Autoregressions            |
#> **************************************************|
#>  Gibbs sampler for the BVAR model                 |
#> **************************************************|
#>  Progress of the MCMC simulation for 5 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|
pred = forecast(post, 4)                      # forecast 1 year ahead

# workflow with the pipe |>
############################################################
set.seed(123)
us_macro_chan |>
  specify_bvar$new() |>
  estimate(S = 5) |> 
  estimate(S = 5) |> 
  forecast(horizon = 4) -> pred
#> **************************************************|
#> bvars: Forecasting with Large                     |
#>        Bayesian Vector Autoregressions            |
#> **************************************************|
#>  Gibbs sampler for the BVAR model                 |
#> **************************************************|
#>  Progress of the MCMC simulation for 5 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|
#> **************************************************|
#> bvars: Forecasting with Large                     |
#>        Bayesian Vector Autoregressions            |
#> **************************************************|
#>  Gibbs sampler for the BVAR model                 |
#> **************************************************|
#>  Progress of the MCMC simulation for 5 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|
```
