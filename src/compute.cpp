
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


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube compute_fitted_values (
    arma::cube&     posterior_A,        // NxKxS
    arma::cube&     posterior_Sigma,    // NxNxS
    arma::mat&      posterior_sigma2,   // TxS
    arma::mat&      X                   // KxT
) {
  
  const int   N = posterior_A.n_rows;
  const int   S = posterior_A.n_slices;
  const int   T = X.n_cols;
  cube        fitted_values(N, T, S, fill::randn);
  
  for (int s=0; s<S; s++) {
    mat L = chol(posterior_Sigma.slice(s), "lower");
    fitted_values.slice(s).each_row() %= trans(sqrt(posterior_sigma2.col(s)));
    fitted_values.slice(s) = posterior_A.slice(s) * X + L * fitted_values.slice(s); 
  } // END s loop
  
  return fitted_values;
} // END compute_fitted_values


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube compute_shocks (
    arma::cube&     posterior_A,    // (N, K, S)
    arma::mat&      Y,              // NxT dependent variables
    arma::mat&      X               // KxT dependent variables
) {
  
  const int       N = Y.n_rows;
  const int       T = Y.n_cols;
  const int       S = posterior_A.n_slices;
  cube            shocks(N, T, S);
  
  for (int s=0; s<S; s++) {
    shocks.slice(s)    =  (Y - posterior_A.slice(s) * X);
  } // END s loop
  
  return shocks;
} // END compute_shocks
