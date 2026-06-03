# Computes posterior draws of the forecast error variance decomposition

Each of the draws from the posterior estimation of the Vector
Autoregression is transformed into a draw from the posterior
distribution of the forecast error variance decomposition.

## Usage

``` r
# S3 method for class 'PosteriorBVAR'
compute_variance_decompositions(posterior, horizon)
```

## Arguments

- posterior:

  posterior estimation outcome obtained by running the `estimate`
  function.

- horizon:

  a positive integer number denoting the forecast horizon for the
  forecast error variance decomposition computations.

## Value

An object of class `PosteriorFEVD`, that is, an `NxNx(horizon+1)xS`
array with attribute `PosteriorFEVD` containing `S` draws of the
forecast error variance decomposition.

## References

Kilian, L., & Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In:
Structural vector autoregressive analysis. Cambridge University Press.

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
fevd = compute_variance_decompositions(post, horizon = 4)

# workflow with the pipe |>
############################################################
us_macro_chan |>
  specify_bvar$new() |>
  estimate(S = 5) |> 
  estimate(S = 10) |> 
  compute_variance_decompositions(horizon = 4) -> fevd
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
