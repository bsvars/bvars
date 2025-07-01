#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat do_rwishart(const arma::mat& Sigma, int df) {
  return wishrnd(Sigma, df);
}