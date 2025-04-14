#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::vec nicerig2 (
  const int     n,    // a positive integer - number of draws
  const double  s,    // a positive scale parameter
  const double  nu    // a positive shape parameter
) {
  arma::vec rig2 = s / chi2rnd( nu, n );
  return rig2;
}