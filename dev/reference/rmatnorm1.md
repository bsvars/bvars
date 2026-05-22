# Samples random numbers from the matrix-variate normal distribution

Samples random numbers from the matrix-variate normal distribution

## Usage

``` r
rmatnorm1(M, V, S)
```

## Arguments

- M:

  a real-valued matrix of the expected values

- V:

  a positive definite symmetric matrix of column-specific covariance

- S:

  a positive definite symmetric matrix of row-specific covariance

## Value

a matrix - a draw from the matrix-variate normal distribution

## Examples

``` r
rmatnorm1(matrix(0, 2, 3), diag(2), diag(3))
#>            [,1]       [,2]       [,3]
#> [1,]  0.6923777  0.4514032 -0.4523040
#> [2,] -0.5496063 -1.7444244  0.0337422
```
