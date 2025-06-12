
#include <RcppArmadillo.h>
#include "rmgig.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_V(
    arma::mat&      aux_V,          // (K,K) matrix
    const arma::mat&      aux_A,          // (N,K) matrix
    const arma::mat&      aux_Sigma_inv,  // (N,N) matrix
    const Rcpp::List&     prior           // a list of prior parameters
) {
  
  const int N         = aux_A.n_rows;
  
  const mat     prior_A       = as<mat>(prior["A"]);
  const mat     prior_Psi     = as<mat>(prior["Psi"]);
  const mat     prior_Gamma   = as<mat>(prior["Gamma"]);
  const double  prior_lambda  = as<double>(prior["lambda"]);
  
  double  full_lambda = prior_lambda - 0.5 * N;
  mat     AA_tmp      = trans(aux_A - prior_A);
  mat     full_Gamma  = prior_Gamma + AA_tmp * aux_Sigma_inv * AA_tmp.t();
  mat     draw        = do_rmgig1(aux_V, full_lambda, prior_Psi, full_Gamma);
  
  return draw;
} // END sample_V