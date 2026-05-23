# Computes posterior draws of shocks

Each of the draws from the posterior estimation of the BAVR model is
transformed into a draw from the posterior distribution of the shocks.

## Usage

``` r
# S3 method for class 'PosteriorBVAR'
compute_shocks(posterior)
```

## Arguments

- posterior:

  posterior estimation outcome - an object of class `PosteriorBVAR`
  obtained by running the `estimate` function.

## Value

An object of class `PosteriorShocks`, that is, an `NxTxS` array with
attribute `PosteriorShocks` containing `S` draws of the shocks.

## Author

Tomasz Woźniak <wozniak.tom@pm.me>

## Examples

``` r
# simple workflow
############################################################
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
shoc = compute_shocks(post)                   # compute shocks

# workflow with the pipe |>
############################################################
us_macro_chan |>
  specify_bvar$new() |>
  estimate(S = 5) |> 
  estimate(S = 10) |> 
  compute_shocks() -> shoc
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
#>  Progress of the MCMC simulation for 10 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|
  
```
