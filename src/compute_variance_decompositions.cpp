
#include <RcppArmadillo.h>
#include <bsvars.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube Sigma2B (
    arma::cube&   posterior_Sigma_c    // (N, N, S)
) {
  
  const int   S = posterior_Sigma_c.n_slices;
  const int   N = posterior_Sigma_c.n_rows;
  
  cube        posterior_B_c(N, N, S);
  mat L;
  
  for (int s = 0; s < S; s++) {
    L                       = chol(posterior_Sigma_c.slice(s), "lower");
    posterior_B_c.slice(s)  = inv(trimatl(L));
  } // END s loop
  
  return posterior_B_c;
} // END Sigma2B


// [[Rcpp::interfaces(cpp,r)]]
// [[Rcpp::export]]
arma::field<arma::cube> compute_variance_decompositions (
    arma::cube&   posterior_Sigma,    // (S)(N, N, S)
    arma::cube&   posterior_A,        // (S)(K, N, S)
    const int     horizon,
    const int     p
) {
  
  const int   S = posterior_Sigma.n_slices;
  
  cube  posterior_B           = Sigma2B( posterior_Sigma );
  field<cube> posterior_irfs  = bsvars::bsvars_ir( posterior_B, posterior_A, horizon, p);
  field<cube> fevds           = bsvars::bsvars_fevd_homosk( posterior_irfs );
  
  return fevds;
} // END compute_variance_decompositions
