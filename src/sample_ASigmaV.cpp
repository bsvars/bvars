
#include <RcppArmadillo.h>
#include <bsvars.h>
#include "rmgig.h"
#include "rmatnorm.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_V_mgig(
    arma::mat&            aux_V,          // (K,K) matrix
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
} // END sample_V_mgig



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec sample_V_gig(
    const arma::mat&      aux_A,          // (N,K) matrix
    const arma::mat&      aux_Sigma_inv,  // (N,N) matrix
    const Rcpp::List&     prior           // a list of prior parameters
) {
  
  const int N         = aux_A.n_rows;
  const int K         = aux_A.n_cols;
  
  const mat     prior_A       = as<mat>(prior["A"]);
  const double  prior_psi     = as<double>(prior["psi"]);
  const double  prior_chi     = as<double>(prior["chi"]);
  const double  prior_lambda  = as<double>(prior["lambda"]);
  
  double  full_lambda = prior_lambda - 0.5 * N;
  mat     AA_tmp      = trans(aux_A - prior_A);
  vec     draw(K);
  double  full_chi    = 0;
  
  for (int i=0; i<K; i++) {
    full_chi          = prior_chi + as_scalar(AA_tmp.row(i) * aux_Sigma_inv * trans(AA_tmp.row(i)));
    draw(i)           = bsvars::do_rgig1(full_lambda, full_chi, prior_psi);
  }
  
  return draw;
} // END sample_V_gig



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
  mat   Omega_inv         = diagmat(aux_Omega_diag_inv);
  
  mat   A_bar_tmp         = Y * Omega_inv * X.t() + prior_A * aux_V_inv;
  mat   V_bar_inv         = X * Omega_inv * X.t() + aux_V_inv;
  mat   V_bar             = inv_sympd(V_bar_inv);
  mat   A_bar             = A_bar_tmp * V_bar;
  
  mat   S_bar         = prior_S;
  S_bar              += Y * Omega_inv * Y.t();
  S_bar              += prior_A * aux_V_inv * prior_A.t();
  S_bar              += A_bar * V_bar_inv * A_bar.t();
  double nu_bar       = prior_nu + T;
  
  mat draw_Sigma      = iwishrnd(S_bar, nu_bar);
  mat draw_A          = do_rmatnorm1(A_bar, draw_Sigma, V_bar);
  
  field<mat> out(2);
  out(0) = draw_A;
  out(1) = draw_Sigma;
  
  return out;
} // END sample_ASigma




