
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// the function below is the rmatnorm_cpp function copied from bsvarSIGNs package on 2025-05-20 and modified subsequently
// [[Rcpp:interface(cpp,r)]]
// [[Rcpp::export]]
arma::mat do_rmatnorm1(     // vec(X) ~N(vec(M), S %x% V)
    const arma::mat& M,     // mean
    const arma::mat& V,     // column-specific covariance
    const arma::mat& S      // row-specific covariance
) {
  
  mat X = mat(size(M), fill::randn);
  return M + chol(V).t() * X * chol(S);
} // END do_rmatnorm1
