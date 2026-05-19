#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat do_rinvwishart(const arma::mat& Psi, int df) {
  return iwishrnd(Psi, df);
}