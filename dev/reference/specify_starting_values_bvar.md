# R6 Class Representing `StartingValuesBVAR`

The class `StartingValuesBVAR` presents starting values for the BVAR
model.

## Public fields

- `A`:

  an `NxK` matrix of starting values for the autoregressive matrix
  \\A\\.

- `Sigma`:

  an `NxN` matrix of starting values for the error term covariance
  \\\Sigma\\.

- `V`:

  a `KxK` matrix of starting values for the prior equation-specific
  covariance \\V\\ of the hierarchical prior distribution for matrix
  \\A\\.

- `h`:

  an `T`-vector with the starting values of the log-volatility
  processes.

- `rho`:

  a scalalr for the SV autoregressive parameter.

- `omega`:

  a scalar for the SV process conditional standard deviation.

- `sigma2v`:

  a scalar for SV process conditional variances.

- `S`:

  a `T` integer vector with the auxiliary mixture component indicator.

- `sigma2_omega`:

  a scalar for the variance of the zero-mean normal prior for
  \\\omega\\.

- `s_`:

  a positive scalar with the scale of the gamma prior of the
  hierarchical prior for \\\sigma^2\_{\omega}\\.

- `lambda`:

  a `T`-vetor of starting values for latent variable.

- `df`:

  a scalar greater than 2 with the starting value for the degrees of
  freedom parameter of the Student-t conditional distribution of error
  term.

## Methods

### Public methods

- [`specify_starting_values_bvar$new()`](#method-StartingValuesBVAR-new)

- [`specify_starting_values_bvar$get_starting_values()`](#method-StartingValuesBVAR-get_starting_values)

- [`specify_starting_values_bvar$set_starting_values()`](#method-StartingValuesBVAR-set_starting_values)

- [`specify_starting_values_bvar$clone()`](#method-StartingValuesBVAR-clone)

------------------------------------------------------------------------

### Method `new()`

Create new starting values `StartingValuesBVAR`.

#### Usage

    specify_starting_values_bvar$new(
      N,
      p,
      T,
      d = 0,
      is_homoskedastic = TRUE,
      is_normal = TRUE,
      ar_sigma2 = rep(1, N),
      kappa = c(0.2^2, 10^2)
    )

#### Arguments

- `N`:

  a positive integer - the number of dependent variables in the model.

- `p`:

  a positive integer - the autoregressive lag order of the BVAR model.

- `T`:

  a positive integer - the number of time periods in the data.

- `d`:

  a positive integer - the number of `exogenous` variables in the model.

- `is_homoskedastic`:

  a logical scalar - if `TRUE` the model assumes homoskedastic errors,
  otherwise it assumes stochastic volatility.

- `is_normal`:

  a logical scalar - if `TRUE` the model assumes normal error term,
  otherwise, it assumes Student-t errors.

- `ar_sigma2`:

  a positive `N`-vector with the autoregressive variance estimates for
  each variable to be used in the Minnesota prior for the autoregressive
  parameters.

- `kappa`:

  a positive `2`-vector with the hyperparameters of the Minnesota prior
  for the autoregressive parameters - the first element is the overall
  tightness hyperparameter, while the second element is the tightness of
  the prior on the constant and exogenous variable coefficients.

#### Returns

Starting values `StartingValuesBVAR`.

#### Examples

    # starting values for a 3-variable BVAR model
    sv = specify_starting_values_bvar$new(N = 3, p = 4, T = 100)

------------------------------------------------------------------------

### Method `get_starting_values()`

Returns the elements of the starting values `StartingValuesBVAR` as a
`list`.

#### Usage

    specify_starting_values_bvar$get_starting_values()

#### Examples

    # starting values for a 3-variable BVAR model
    sv = specify_starting_values_bvar$new(N = 3, p = 4, T = 100)
    sv$get_starting_values()   # show starting values as list

------------------------------------------------------------------------

### Method `set_starting_values()`

Sets the elements of the starting values `StartingValuesBVAR` to
provided values.

#### Usage

    specify_starting_values_bvar$set_starting_values(last_draw)

#### Arguments

- `last_draw`:

  a list containing the last draw of elements `A` - a `KxN` matrix,
  `Sigma` - an `NxN` matrix, and `V` - a `KxK` matrix.

#### Returns

An object of class `StartingValuesBVAR` including the last draw of the
current MCMC as the starting value to be passed to the continuation of
the MCMC estimation using
[`estimate()`](https://bsvars.org/bsvars/reference/estimate.html).

#### Examples

    # starting values for a 3-variable BVAR model
    sv = specify_starting_values_bvar$new(N = 3, p = 4, T = 100)

    # Modify the starting values by:
    sv_list = sv$get_starting_values()   # getting them as list
    sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
    sv$set_starting_values(sv_list)      # providing to the class object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    specify_starting_values_bvar$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# starting values for a 3-variable BVAR model.
sv = specify_starting_values_bvar$new(N = 3, p = 4, T = 100)


## ------------------------------------------------
## Method `specify_starting_values_bvar$new`
## ------------------------------------------------

# starting values for a 3-variable BVAR model
sv = specify_starting_values_bvar$new(N = 3, p = 4, T = 100)


## ------------------------------------------------
## Method `specify_starting_values_bvar$get_starting_values`
## ------------------------------------------------

# starting values for a 3-variable BVAR model
sv = specify_starting_values_bvar$new(N = 3, p = 4, T = 100)
sv$get_starting_values()   # show starting values as list
#> $A
#>           [,1]      [,2]        [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#> [1,] 0.5741127 0.0000000 0.000000000    0    0    0    0    0    0     0     0
#> [2,] 0.0000000 0.2358854 0.000000000    0    0    0    0    0    0     0     0
#> [3,] 0.0000000 0.0000000 0.007971651    0    0    0    0    0    0     0     0
#>      [,12] [,13]
#> [1,]     0     0
#> [2,]     0     0
#> [3,]     0     0
#> 
#> $Sigma
#>           [,1]       [,2]      [,3]
#> [1,] 0.2874713 0.00000000 0.0000000
#> [2,] 0.0000000 0.02791551 0.0000000
#> [3,] 0.0000000 0.00000000 0.2883005
#> 
#> $V
#>       [,1] [,2] [,3] [,4] [,5] [,6]        [,7]        [,8]        [,9]  [,10]
#>  [1,] 0.04 0.00 0.00 0.00 0.00 0.00 0.000000000 0.000000000 0.000000000 0.0000
#>  [2,] 0.00 0.04 0.00 0.00 0.00 0.00 0.000000000 0.000000000 0.000000000 0.0000
#>  [3,] 0.00 0.00 0.04 0.00 0.00 0.00 0.000000000 0.000000000 0.000000000 0.0000
#>  [4,] 0.00 0.00 0.00 0.01 0.00 0.00 0.000000000 0.000000000 0.000000000 0.0000
#>  [5,] 0.00 0.00 0.00 0.00 0.01 0.00 0.000000000 0.000000000 0.000000000 0.0000
#>  [6,] 0.00 0.00 0.00 0.00 0.00 0.01 0.000000000 0.000000000 0.000000000 0.0000
#>  [7,] 0.00 0.00 0.00 0.00 0.00 0.00 0.004444444 0.000000000 0.000000000 0.0000
#>  [8,] 0.00 0.00 0.00 0.00 0.00 0.00 0.000000000 0.004444444 0.000000000 0.0000
#>  [9,] 0.00 0.00 0.00 0.00 0.00 0.00 0.000000000 0.000000000 0.004444444 0.0000
#> [10,] 0.00 0.00 0.00 0.00 0.00 0.00 0.000000000 0.000000000 0.000000000 0.0025
#> [11,] 0.00 0.00 0.00 0.00 0.00 0.00 0.000000000 0.000000000 0.000000000 0.0000
#> [12,] 0.00 0.00 0.00 0.00 0.00 0.00 0.000000000 0.000000000 0.000000000 0.0000
#> [13,] 0.00 0.00 0.00 0.00 0.00 0.00 0.000000000 0.000000000 0.000000000 0.0000
#>        [,11]  [,12] [,13]
#>  [1,] 0.0000 0.0000     0
#>  [2,] 0.0000 0.0000     0
#>  [3,] 0.0000 0.0000     0
#>  [4,] 0.0000 0.0000     0
#>  [5,] 0.0000 0.0000     0
#>  [6,] 0.0000 0.0000     0
#>  [7,] 0.0000 0.0000     0
#>  [8,] 0.0000 0.0000     0
#>  [9,] 0.0000 0.0000     0
#> [10,] 0.0000 0.0000     0
#> [11,] 0.0025 0.0000     0
#> [12,] 0.0000 0.0025     0
#> [13,] 0.0000 0.0000   100
#> 
#> $h
#>   [1]  1.696982e-03  3.001271e-03  1.109257e-02 -7.198301e-03  1.879593e-02
#>   [6] -4.177957e-03 -3.569055e-03 -1.859122e-03  9.036993e-03  1.007237e-02
#>  [11] -2.823084e-03  4.330427e-03  9.793907e-03 -2.217763e-03  1.629873e-04
#>  [16] -4.195563e-03  4.768723e-03 -5.340778e-03  1.492208e-02 -2.096973e-03
#>  [21] -1.526562e-03  9.234436e-03  2.052924e-02  4.052080e-03  2.081010e-02
#>  [26] -4.068178e-03  4.924050e-03 -2.678083e-04  2.221764e-02 -2.791158e-03
#>  [31]  2.705488e-02  4.442011e-03  4.657068e-04 -4.619567e-03 -7.104474e-03
#>  [36] -2.275246e-04 -5.306322e-03  5.466885e-03  1.074639e-02 -2.673047e-06
#>  [41]  2.035244e-03  4.382628e-03 -1.304181e-02  3.076142e-03  2.019408e-03
#>  [46]  9.527296e-03  5.320749e-03 -5.603359e-03  1.005445e-02 -3.508077e-03
#>  [51] -4.302493e-03 -1.206231e-02  1.394537e-02  6.091125e-03  9.837046e-03
#>  [56]  8.522565e-04 -5.298947e-03  3.313622e-02 -1.016320e-02 -1.201424e-02
#>  [61]  1.277549e-02 -2.384533e-03  8.251479e-03 -1.434804e-02  3.887690e-03
#>  [66] -1.706339e-02 -4.080955e-03 -8.146259e-03  3.395826e-03 -7.657576e-04
#>  [71] -1.406630e-03 -4.892756e-03  2.250724e-03  1.017381e-02  8.156466e-03
#>  [76]  4.473520e-05  1.244433e-02 -4.920910e-03  8.982512e-03  6.146918e-03
#>  [81] -4.003597e-03  9.680298e-05  6.290504e-03  1.133209e-02 -2.520703e-03
#>  [86]  2.116707e-03  5.049962e-03  3.973690e-03  7.809080e-03 -4.416303e-03
#>  [91] -5.216638e-03  1.013842e-02  5.594613e-04 -1.436022e-02 -5.100833e-04
#>  [96]  1.898591e-02 -1.569562e-02 -2.370156e-02 -3.606970e-03 -5.684785e-04
#> 
#> $rho
#> [1] 0.5
#> 
#> $omega
#> [1] 0.1
#> 
#> $sigma2v
#> [1] 0.01
#> 
#> $S
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> $sigma2_omega
#> [1] 1
#> 
#> $s_
#> [1] 0.05
#> 
#> $lambda
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> $df
#> [1] 30
#> 


## ------------------------------------------------
## Method `specify_starting_values_bvar$set_starting_values`
## ------------------------------------------------

# starting values for a 3-variable BVAR model
sv = specify_starting_values_bvar$new(N = 3, p = 4, T = 100)

# Modify the starting values by:
sv_list = sv$get_starting_values()   # getting them as list
sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
sv$set_starting_values(sv_list)      # providing to the class object
```
