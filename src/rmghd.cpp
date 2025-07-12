#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Declare  MGIG sampler defined in rmgig.cpp
arma::mat do_rmgig1(
    arma::mat& aux_Sigma,       
    const double lambda,
    const arma::mat Psi,
    const arma::mat Gamma);

// [[Rcpp::export]]
arma::mat sample_mgh(const arma::mat& M,            // (n x p) mean
                     const arma::mat& U,            // (n x n) row covariance
                     const arma::mat& aux_Sigma,    // (p x p) SPD matrix
                     const double lambda,           // MGIG shape
                     const arma::mat& Psi,          // MGIG Psi
                     const arma::mat& Gamma) {      // MGIG Gamma
  
  int n = U.n_rows;
  int p = Psi.n_rows;
  
  arma::mat Sigma = do_rmgig1(const_cast<arma::mat&>(aux_Sigma), lambda, Psi, Gamma);
  
  // Sample A | Sigma ~ MN(M, U, Sigma)
  arma::mat Z = arma::randn(n, p); //an 𝑛x𝑝 matrix of independent standard normal random variables
  arma::mat A = M + chol(U, "lower") * Z * chol(Sigma, "upper").t();
  return A;
}
