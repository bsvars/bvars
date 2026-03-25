

#' @title Samples random numbers from the matrix-variate normal distribution
#' @param M a real-valued matrix of the expected values
#' @param V a positive definite symmetric matrix of column-specific covariance
#' @param S a positive definite symmetric matrix of row-specific covariance
#' @return  a matrix - a draw from the matrix-variate normal distribution
#' @examples
#' rmatnorm1(matrix(0, 2, 3), diag(2), diag(3))
#' @export
rmatnorm1 <- function(M, V, S) {
  
  # Check if the input matrices are positive definite
  stopifnot("M, V, and S must be matrices."
            = is.matrix(M) && is.matrix(V) && is.matrix(S))
  stopifnot("V, and S must be square matrices, and dimaensions of M must match."
            = nrow(M) == nrow(V) && ncol(M) == nrow(S) && nrow(V) == ncol(V) && nrow(S) == ncol(S))
  stopifnot("V and S must be positive definite matrices."
            = any(eigen(V)$values > 0) && any(eigen(S)$values > 0))
  stopifnot("M, V, and S must not contain NA values."
            = !any(is.na(M)) && !any(is.na(V)) && !any(is.na(S)))
  
  out = .Call(`_bvars_do_rmatnorm1`, M, V, S)
  return(out)
}