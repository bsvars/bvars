
#include <RcppArmadillo.h>
#include <bsvars.h>
#include "progress.hpp"

#include "sample_ASigmaV.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bvar_gig_cpp(
    const int&  S,                        // number of draws from the posterior
    const arma::mat&  Y,                  // NxT dependent variables
    const arma::mat&  X,                  // KxT dependent variables
    const Rcpp::List& prior,              // a list of priors
    const Rcpp::List& starting_values,    // a list of starting values
    const int         thin = 100,         // introduce thinning
    const bool        show_progress = true
) {
  
  std::string oo = "";
  if ( thin != 1 ) {
    oo      = bsvars::ordinal(thin) + " ";
  }
  
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, S, 50));
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << "bvarGIGs: Bayesian Vector Autoregressions         |" << endl;
    Rcout << "          with Flexible Prior Distributions       |" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Gibbs sampler for the BVAR model                 |" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Progress of the MCMC simulation for " << S << " draws" << endl;
    Rcout << "    Every " << oo << "draw is saved via MCMC thinning" << endl;
    Rcout << " Press Esc to interrupt the computations" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress p(50, show_progress);
  
  const int N         = Y.n_rows;
  const int K         = X.n_rows;
  const int T         = Y.n_cols;
  
  mat   aux_A         = as<mat>(starting_values["A"]);
  mat   aux_Sigma     = as<mat>(starting_values["Sigma"]);
  mat   aux_Sigma_inv = inv_sympd(aux_Sigma);
  mat   aux_V         = as<mat>(starting_values["V"]);
  mat   aux_V_inv     = inv_sympd(aux_V);
  vec   aux_Omega_diag_inv(T, fill::ones);
    
  const int   SS      = floor(S / thin);
  
  cube  posterior_A(N, K, SS);
  cube  posterior_Sigma(N, N, SS);
  cube  posterior_V(K, K, SS);
  
  int   ss = 0;
  
  for (int s=0; s<S; s++) {
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    aux_V                 = sample_V_mgig( aux_V, aux_A, aux_Sigma_inv, prior );
    aux_V_inv             = inv_sympd(aux_V);
    
    field<mat> aux_ASigma = sample_ASigma( Y, X, aux_V_inv, aux_Omega_diag_inv, prior );
    aux_A                 = aux_ASigma(0);
    aux_Sigma             = aux_ASigma(1);
    
    if (s % thin == 0) {
      posterior_A.slice(ss)     = aux_A;
      posterior_Sigma.slice(ss) = aux_Sigma;
      posterior_V.slice(ss)     = aux_V;
      ss++;
    }
  } // END s loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["A"]        = aux_A,
      _["Sigma"]    = aux_Sigma,
      _["V"]        = aux_V
    ),
    _["posterior"]  = List::create(
      _["A"]        = posterior_A,
      _["Sigma"]    = posterior_Sigma,
      _["V"]        = posterior_V
    )
  );
} // END bvar_gig_cpp
