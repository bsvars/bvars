
#include <RcppArmadillo.h>
#include <bsvars.h>
#include "rmatnorm.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::mat> sample_ASigma(
    const arma::mat&      Y,              // (N,T) matrix
    const arma::mat&      X,              // (K,T) matrix
    arma::mat&            aux_V_inv,      // (K,K) matrix
    arma::vec&            aux_Omega_diag_inv, // (T) matrix
    const Rcpp::List&     prior           // a list of prior parameters
) {
  
  const int T             = aux_Omega_diag_inv.n_elem;
  
  const mat     prior_A   = as<mat>(prior["A"]);
  const mat     prior_S   = as<mat>(prior["S"]);
  const double  prior_nu  = as<double>(prior["nu"]);
  vec   Omega_diag_sqrt_inv = sqrt(aux_Omega_diag_inv);
  
  mat   Ys            = Y.each_row() % Omega_diag_sqrt_inv.t();
  mat   Xs            = X.each_row() % Omega_diag_sqrt_inv.t();
  mat   A_bar_tmp     = Ys * Xs.t() + prior_A * aux_V_inv;
  mat   V_bar_inv     = Xs * Xs.t() + aux_V_inv;
  V_bar_inv           = 0.5 * (V_bar_inv + V_bar_inv.t());
  mat   V_bar         = inv_sympd(V_bar_inv);
  mat   A_bar         = A_bar_tmp * V_bar;
  
  mat   S_bar         = prior_S;
  S_bar              += Ys * Ys.t();
  S_bar              += prior_A * aux_V_inv * prior_A.t();
  S_bar              -= A_bar * V_bar_inv * A_bar.t();
  S_bar               = 0.5 * (S_bar + S_bar.t());
  double nu_bar       = prior_nu + T;
  
  mat draw_Sigma      = iwishrnd(S_bar, nu_bar);
  mat draw_A          = do_rmatnorm1(A_bar, draw_Sigma, V_bar);
  
  field<mat> out(2);
  out(0) = draw_A;
  out(1) = draw_Sigma;
  
  return out;
} // END sample_ASigma



