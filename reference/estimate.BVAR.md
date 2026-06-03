# Bayesian Estimation via Gibbs sampler of a Bayesian VAR with a Flexible Error Term Specification

Estimates the model by Chan (2020)
\<doi:10.1080/07350015.2018.1451336\>, that is, a Bayesian Vector
Autoregression with Minnesota priors and a flexible structure of the
error term specification. The latter includes: conditional multivariate
normal or Student’s t distributions, as well as homoskedastic or
heteroskedastic specifications with a common volatility modelled by
centred or non-centred Stochastic Volatility. The estimation is
conducted via an efficient Gibbs sampler employing frontier numerical
techniques and algorithms written in C++ for excellent computational
speed.

## Usage

``` r
# S3 method for class 'BVAR'
estimate(specification, S, thin = 1, show_progress = TRUE)
```

## Arguments

- specification:

  an object of class `BVAR` generated using the `specify_bvar$new()`
  function.

- S:

  a positive integer, the number of posterior draws to be generated

- thin:

  a positive integer, specifying the frequency of MCMC output thinning

- show_progress:

  a logical value, if `TRUE` the estimation progress bar is visible

## Value

An object of class `PosteriorBVAR` containing the Bayesian estimation
output and containing two elements:

`posterior` a list with a collection of `S` draws from the posterior
distribution generated via Gibbs sampler containing:

- A:

  an `NxKxS` array with the posterior draws for matrix \\A\\

- Sigma:

  an `NxNxS` array with the posterior draws for matrix \\\Sigma\\

- V:

  a `KxKxS` array with the posterior draws for the hyper-parameter
  matrix \\\Sigma\\

`last_draw` an object of class `BVAR` with the last draw of the current
MCMC run as the starting value to be passed to the continuation of the
MCMC estimation using
[`estimate()`](https://bsvars.org/bsvars/reference/estimate.html).

## References

Barndorff-Nielsen, Blaesild, Jensen, Jorgensen (1982) Exponential
Transformation Models, Proceedings of the Royal Society of London. A.
Mathematical and Physical Sciences, 379, 41–-65,
\<doi:10.1098/rspa.1982.0004\>.

Chan (2020) Large Bayesian VARs: A Flexible Kronecker Error Covariance
Structure, Journal of Business and Economic Statistics, 38(1), 68–79,
\<doi:10.1080/07350015.2018.1451336\>.

Hamura, Irie, Sugasawa (2024) Gibbs Sampler for Matrix Generalized
Inverse Gaussian Distributions, Journal of Computational and Graphical
Statistics, 33(2), 331–340, \<doi:10.1080/10618600.2023.2258186\>.

Thabane, Safiul Haq (2004) On the Matrix-Variate Generalized Hyperbolic
Distribution and Its Bayesian Applications, Statistics: A Journal of
Theoretical and Applied Statistics, 38(6), 511–526,
\<doi:10.1080/02331880412331319279\>.

Woźniak (2016) Bayesian Vector Autoregressions, Australian Economic
Review, 49(3), 365–380, \<doi:10.1111/1467-8462.12179\>.

## See also

[`specify_bvar`](http://bsvars.org/bvars/reference/specify_bvar.md),
[`specify_posterior_bvar`](http://bsvars.org/bvars/reference/specify_posterior_bvar.md)

## Author

Rui Liu <rl3023@columbia.edu>, Andres Ramirez Hassan
<aramir21@gmail.com> & Tomasz Woźniak <wozniak.tom@pm.me>

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

# workflow with the pipe |>
############################################################
us_macro_chan |>
  specify_bvar$new() |>
  estimate(S = 5) |> 
  estimate(S = 10) -> post
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
