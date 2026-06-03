# R6 Class Representing `PriorBVAR`

The class `PriorBVAR` presents a prior specification for the BVAR model.

\#' **The Model.** All the BVAR models in this package are specified by
two equations, including the reduced form equation: \$\$Y = AX + E\$\$
where \\Y\\ is an `NxT` matrix of dependent variables, \\X\\ is a `KxT`
matrix of explanatory variables, \\E\\ is an `NxT` matrix of reduced
form error terms, and \\A\\ is an `NxK` matrix of autoregressive slope
coefficients and parameters on deterministic terms in \\X\\.

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

## Public fields

- `A`:

  a real-valued `NxK` matrix, the mean matrix \\A_0\\ of the
  matrix-variate normal prior distribution for the parameter matrix
  \\A\\.

- `S`:

  a `NxN` positive definite scale matrix \\S_0\\ of the Inverse Wishart
  prior distribution for the error term covariance matrix \\\Sigma\\.

- `nu`:

  a positive scalar, shape parameter \\\nu_0\\ of the Inverse Wishart
  prior distribution for the error term covariance matrix \\\Sigma\\.

- `Psi`:

  a `KxK` scale matrix \\\Psi_0\\ of the matrix generalized inverse
  Gaussian distribution for the equation-specific prior covariance \\V\\

- `Gamma`:

  a `KxK` scale matrix \\\Gamma_0\\ of the matrix generalized inverse
  Gaussian distribution for the equation-specific prior covariance \\V\\

- `lambda`:

  a positive scalar shape parameter \\\lambda_0\\ of the matrix
  generalized inverse Gaussian distribution for the equation-specific
  prior covariance \\V\\

- `sv_a`:

  a positive scalar, the shape parameter of the gamma prior in the
  hierarchical prior for the common stochastic volatility.

- `sv_s`:

  a positive scalar, the scale parameter of the gamma prior in the
  hierarchical prior for the common stochastic volatility.

## Methods

### Public methods

- [`specify_prior_bvar$new()`](#method-PriorBVAR-new)

- [`specify_prior_bvar$get_prior()`](#method-PriorBVAR-get_prior)

- [`specify_prior_bvar$clone()`](#method-PriorBVAR-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new prior specification `PriorBVAR`.

#### Usage

    specify_prior_bvar$new(
      N,
      p,
      d = 0,
      stationary = rep(FALSE, N),
      is_homoskedastic = TRUE
    )

#### Arguments

- `N`:

  a positive integer - the number of dependent variables in the model.

- `p`:

  a positive integer - the autoregressive lag order of the VAR model.

- `d`:

  a positive integer - the number of `exogenous` variables in the model.

- `stationary`:

  an `N` logical vector - its element set to `FALSE` sets the prior mean
  for the autoregressive parameters of the `N`th equation to the white
  noise process, otherwise to random walk.

- `is_homoskedastic`:

  a logical scalar - if `TRUE` the model assumes homoskedastic errors,
  otherwise it assumes stochastic volatility.

#### Returns

A new prior specification `PriorBVAR`.

#### Examples

    # a prior for 3-variable example with one lag and stationary data
    prior = specify_prior_bvar$new(N = 3, p = 1, stationary = rep(TRUE, 3))
    prior$A # show autoregressive prior mean

------------------------------------------------------------------------

### Method `get_prior()`

Returns the elements of the prior specification `PriorBVAR` as a `list`.

#### Usage

    specify_prior_bvar$get_prior()

#### Examples

    # a prior for 3-variable example with four lags
    prior = specify_prior_bvar$new(N = 3, p = 4)
    prior$get_prior() # show the prior as list

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    specify_prior_bvar$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
prior = specify_prior_bvar$new(N = 3, p = 1)  # a prior for 3-variable example with one lag
prior$A                                        # show autoregressive prior mean
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0
#> [2,]    0    1    0    0
#> [3,]    0    0    1    0


## ------------------------------------------------
## Method `specify_prior_bvar$new`
## ------------------------------------------------

# a prior for 3-variable example with one lag and stationary data
prior = specify_prior_bvar$new(N = 3, p = 1, stationary = rep(TRUE, 3))
prior$A # show autoregressive prior mean
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    0    0    0    0
#> [3,]    0    0    0    0


## ------------------------------------------------
## Method `specify_prior_bvar$get_prior`
## ------------------------------------------------

# a prior for 3-variable example with four lags
prior = specify_prior_bvar$new(N = 3, p = 4)
prior$get_prior() # show the prior as list
#> $A
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#> [1,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#> [2,]    0    1    0    0    0    0    0    0    0     0     0     0     0
#> [3,]    0    0    1    0    0    0    0    0    0     0     0     0     0
#> 
#> $S
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
#> $nu
#> [1] 6
#> 
#> $Psi
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#>  [2,]    0    1    0    0    0    0    0    0    0     0     0     0     0
#>  [3,]    0    0    1    0    0    0    0    0    0     0     0     0     0
#>  [4,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#>  [5,]    0    0    0    0    1    0    0    0    0     0     0     0     0
#>  [6,]    0    0    0    0    0    1    0    0    0     0     0     0     0
#>  [7,]    0    0    0    0    0    0    1    0    0     0     0     0     0
#>  [8,]    0    0    0    0    0    0    0    1    0     0     0     0     0
#>  [9,]    0    0    0    0    0    0    0    0    1     0     0     0     0
#> [10,]    0    0    0    0    0    0    0    0    0     1     0     0     0
#> [11,]    0    0    0    0    0    0    0    0    0     0     1     0     0
#> [12,]    0    0    0    0    0    0    0    0    0     0     0     1     0
#> [13,]    0    0    0    0    0    0    0    0    0     0     0     0     1
#> 
#> $Gamma
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#>  [2,]    0    1    0    0    0    0    0    0    0     0     0     0     0
#>  [3,]    0    0    1    0    0    0    0    0    0     0     0     0     0
#>  [4,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#>  [5,]    0    0    0    0    1    0    0    0    0     0     0     0     0
#>  [6,]    0    0    0    0    0    1    0    0    0     0     0     0     0
#>  [7,]    0    0    0    0    0    0    1    0    0     0     0     0     0
#>  [8,]    0    0    0    0    0    0    0    1    0     0     0     0     0
#>  [9,]    0    0    0    0    0    0    0    0    1     0     0     0     0
#> [10,]    0    0    0    0    0    0    0    0    0     1     0     0     0
#> [11,]    0    0    0    0    0    0    0    0    0     0     1     0     0
#> [12,]    0    0    0    0    0    0    0    0    0     0     0     1     0
#> [13,]    0    0    0    0    0    0    0    0    0     0     0     0     1
#> 
#> $lambda
#> [1] 4
#> 
#> $sv_a
#> [1] NA
#> 
#> $sv_s
#> [1] NA
#> 
```
