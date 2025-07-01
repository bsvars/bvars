#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Forward declaration if needed
arma::mat do_rmatnorm1(const arma::mat& M, const arma::mat& V, const arma::mat& S);

// [[Rcpp::export]]
arma::mat do_rmatt(const arma::mat& M, const arma::mat& V, const arma::mat& S, const int df) {
  // Generate matrix-normal random matrix with zero mean
  arma::mat Z = do_rmatnorm1(arma::zeros(size(M)), V, S);
  
  // Chi-squared random variable
  double G = R::rchisq(df);
  
  // Matrix t draw
  return M + Z / std::sqrt(G / df);
}
