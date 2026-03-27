#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat do_rinvwishart(const arma::mat& Psi, int df) {
  
  // Sample W ~ Wishart(df, inv(Psi))
  arma::mat W = wishrnd(inv_sympd(Psi), df);
  
  // Return inverse: W^{-1} ~ InvWishart(df, Psi)
  return inv_sympd(W);
}