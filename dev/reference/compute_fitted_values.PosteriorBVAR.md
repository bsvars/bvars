# Computes posterior draws from data predictive density

Each of the draws from the posterior estimation of the BVAR model is
transformed into a draw from the data predictive density.

## Usage

``` r
# S3 method for class 'PosteriorBVAR'
compute_fitted_values(posterior)
```

## Arguments

- posterior:

  posterior estimation outcome - an object of class `PosteriorBVAR`
  obtained by running the `estimate` function.

## Value

An object of class `PosteriorFitted`, that is, an `NxTxS` array with
attribute `PosteriorFitted` containing `S` draws from the data
predictive density.

## Author

Tomasz Woźniak <wozniak.tom@pm.me>

## Examples

``` r
spec = specify_bvar$new(us_macro_chan)        # specify the model
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
post = estimate(burn, 10)                     # estimate the model
#> **************************************************|
#> bvars: Forecasting with Large                     |
#>        Bayesian Vector Autoregressions            |
#> **************************************************|
#>  Gibbs sampler for the BVAR model                 |
#> **************************************************|
#>  Progress of the MCMC simulation for 10 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|
fitt = compute_fitted_values(post)            # fitted values

# workflow with the pipe |>
############################################################
us_macro_chan |>
  specify_bvar$new() |>
  estimate(S = 5) |> 
  estimate(S = 5) |> 
  compute_fitted_values() -> fitt
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
