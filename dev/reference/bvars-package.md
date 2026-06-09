# Bayesian Forecasting with Large Vector Autoregressions

Provides fast and efficient procedures for Bayesian estimation and
forecasting using state-of-the-art Vector Autoregressions. This package
includes the model proposed by Chan (2020)
\<doi:10.1080/07350015.2018.1451336\>, that is, a Bayesian Vector
Autoregression with Minnesota priors and a flexible structure of the
error term specification. The latter includes: conditional multivariate
normal or Student’s t distributions, as well as homoskedastic or
heteroskedastic specifications with a common volatility modelled by
centred or non-centred Stochastic Volatility. Additionally, the package
facilitates predictive analyses using density forecasting and
forecast-error variance decompositions. All this is complemented by
simple workflows, useful plots and summary functions, and comprehensive
documentation. The 'bvars' package aligns with R packages 'bsvars' by
Woźniak (2024) \<doi:10.32614/CRAN.package.bsvars\>, 'bsvarSIGNs' by
Wang & Woźniak (2025) \<doi:10.32614/CRAN.package.bsvarSIGNs\>, and
'bpvars' by Woźniak (2025) \<doi:10.32614/CRAN.package.bpvars\>
regarding objects, workflows, and code structure, and they constitute an
integrated toolset.

## Details

**Models.** All the BVAR models in this package are specified by two
equations, including the reduced form equation: \$\$Y = AX + E\$\$ where
\\Y\\ is an `NxT` matrix of dependent variables, \\X\\ is a `KxT` matrix
of explanatory variables, \\E\\ is an `NxT` matrix of reduced form error
terms, and \\A\\ is an `NxK` matrix of autoregressive slope coefficients
and parameters on deterministic terms in \\X\\.

This package assumes that the error matrix follows a matrix normal
distribution: \$\$E \mid X \sim \mathcal{MN}\_{N \times T}(\mathbf{0},
\mathbf{\Sigma}, \mathbf{\Omega})\$\$ where \\\Sigma\\ is the `NxN`
covariance matrix of the error term at time \\t\\, and \\\Omega\\ is a
`TxT` diagonal matrix.

The diagonal elements of \\\Omega\\ determine the specification of the
error term covariance structure. Specifically, the error term at time
\\t\\ follows the multivariate normal distribution \$\$e_t \sim
\mathcal{N}\_N(\mathbf{0}, \sigma_t^2\lambda_t \mathbf{\Sigma})\$\$
where the scalar processes \\\sigma_t^2\\ and \\\lambda_t\\ determine
the diagonal elements of \\\Omega\\. The process \\\sigma_t^2\\
specifies conditional variance and includes three options:

- \\\sigma_t^2 = 1\\:

  homoskedastic error term

- \\\sigma_t^2\\:

  estimated and following non-centred stochastic volatility

- \\\sigma_t^2\\:

  estimated and following centred stochastic volatility

The process \\\lambda_t\\ specifies the conditional distribution of the
error term and includes two options:

- \\\lambda_t = 1\\:

  Gaussian error term specification

- \\\lambda_t\\:

  estimated and following a priori an inverse gamma 2 distribution
  \\\mathcal{IG}2(\nu - 2, \nu)\\, where \\\nu \> 2\\ is a degrees of
  freedom parameter

**Prior distributions.** The autoregressive matrix \\A\\ is assigned
matrix-variate normal distribution: \$\$ \mathbf{A} \mid
\underline{\mathbf{A}}, \mathbf{V}, \boldsymbol{\Sigma} \sim
\mathcal{MN}\_{N \times K}(\underline{\mathbf{A}}, \boldsymbol{\Sigma},
\mathbf{V}) \$\$ with the mean matrix \\\underline{\mathbf{A}}\\, and
covariance matrices \\\boldsymbol{\Sigma}\_{N\times N}\\ and
\\\mathbf{V}\_{K\times K}\\ defining the row- and column-covariance
structures.

This is complemented by the inverse Wishart prior for the error term
covariance \\\boldsymbol{\Sigma}\\: \$\$ \boldsymbol{\Sigma} \mid
\underline{\mathbf{S}}, \underline{\nu} \sim
\mathcal{IW}(\underline{\mathbf{S}}, \underline{\nu}) \$\$ with the
scale matrix \\\underline{\mathbf{S}}\\ and degrees of freedom
\\\underline{\nu}\\.

## Note

This package is currently in active development. Your comments,
suggestions and requests are warmly welcome!

## References

Chan (2020) Large Bayesian VARs: A Flexible Kronecker Error Covariance
Structure, Journal of Business and Economic Statistics, 38(1), 68–79,
\<doi:10.1080/07350015.2018.1451336\>.

## See also

Useful links:

- <https://bsvars.org/bvars/>

- Report bugs at <https://github.com/bsvars/bvars/issues>

## Author

Rui Liu <rl3023@columbia.edu>, Andres Ramirez Hassan
<aramir21@gmail.com>, Tomasz Woźniak <wozniak.tom@pm.me>

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
