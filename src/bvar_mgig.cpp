
#include <RcppArmadillo.h>
#include <bsvars.h>
#include "progress.hpp"

#include "sample_ASigmaV.h"
#include "sample_Omega.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bvar_mgig_cpp(
    const int&        S,                  // number of draws from the posterior
    const arma::mat&  Y,                  // NxT dependent variables
    const arma::mat&  X,                  // KxT dependent variables
    const Rcpp::List& prior,              // a list of priors
    const Rcpp::List& starting_values,    // a list of starting values
    const arma::mat&  sv_aux_mix,         // provide the selected auxiliary mixture components
    const bool        homoskedastic = true, 
    const bool        centred_sv = false, // otherwise non-centred stochastic volatility
    const bool        normal = true,      // otherwise Student-t
    const int         thin = 100,         // introduce thinning
    const bool        show_progress = true
) {
  
  bool debug = false;
  
  std::string oo = "";
  if ( thin != 1 ) {
    oo      = bsvars::ordinal(thin) + " ";
  }
  
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, S, 50));
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << "bvars: Forecasting with Large                     |" << endl;
    Rcout << "       Bayesian Vector Autoregressions            |" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Gibbs sampler for the BVAR model                 |" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Progress of the MCMC simulation for " << S << " draws" << endl;
    Rcout << "    Every " << oo << "draw is saved via MCMC thinning" << endl;
    Rcout << " Press Esc to interrupt the computations" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress p(50, show_progress);
  
  if (debug) Rcout << " initialisation" << endl;
  const int N         = Y.n_rows;
  const int K         = X.n_rows;
  const int T         = Y.n_cols;
  const double ccc    = 0.000000001;      // a constant to make log((u+ccc)^2) feasible
  
  mat     aux_A         = as<mat>(starting_values["A"]);
  mat     aux_Sigma     = as<mat>(starting_values["Sigma"]);
  mat     aux_Sigma_inv = inv_sympd(aux_Sigma);
  mat     aux_V         = as<mat>(starting_values["V"]);
  mat     aux_V_inv     = inv_sympd(aux_V);
  vec     aux_h         = as<vec>(starting_values["h"]);
  double  aux_rho       = as<double>(starting_values["rho"]);
  double  aux_omega     = as<double>(starting_values["omega"]);
  double  aux_sigma2v   = as<double>(starting_values["sigma2v"]);
  uvec    aux_S         = as<uvec>(starting_values["S"]);
  double  aux_sigma2_omega = as<double>(starting_values["sigma2_omega"]);
  double  aux_s_        = as<double>(starting_values["s_"]);
  vec     aux_lambda    = as<vec>(starting_values["lambda"]);
  double  aux_df        = as<double>(starting_values["df"]);
  vec     aux_sigma2(T, fill::ones);
  vec     aux_hetero_inv(T, fill::ones);
  mat     U = Y - aux_A * X;
  mat     aux_L = chol(aux_Sigma, "lower");
  mat     U_std = solve(trimatl(aux_L), U);
  vec     u = trans(sum(U_std)) / N;
  
  if ( !homoskedastic & centred_sv ) {
    aux_sigma2 = exp(aux_h);
    aux_hetero_inv = 1 / (aux_sigma2 % aux_lambda);
  } else if ( !homoskedastic & !centred_sv ) {
    aux_sigma2 = exp(aux_omega * aux_h);
    aux_hetero_inv = 1 / (aux_sigma2 % aux_lambda);
  }
  
  if (debug) Rcout << " post init" << endl;
  const int   SS      = floor(S / thin);
  
  cube  posterior_A(N, K, SS);
  cube  posterior_Sigma(N, N, SS);
  cube  posterior_V(K, K, SS);
  mat   posterior_h(T, SS);
  vec   posterior_rho(SS);
  vec   posterior_omega(SS);
  vec   posterior_sigma2v(SS);
  umat  posterior_S(T, SS);
  vec   posterior_sigma2_omega(SS);
  vec   posterior_s_(SS);
  mat   posterior_sigma2(T, SS);
  mat   posterior_lambda(T, SS);
  vec   posterior_df(SS);
  
  // the initial value for the adaptive_scale is set to the negative inverse of 
  // Hessian for the posterior log_kenel for df evaluated at df = 30
  double    adaptive_scale      = abs(pow(0.25 * T * R::psigamma(15, 1) - T * 29 * pow(28, -2) - 2 * pow(29, -2), -1));
  const vec adptive_alpha_gamma = as<vec>(NumericVector::create(0.44, 0.6));
  int   ss = 0;
  
  // parameters of the SV auxiliary mixture
  if (debug) Rcout << " auxiliary mix" << endl;
  mat aux_mix         = sv_aux_mix;
  if (N > 100) {
    aux_mix           = sv_aux_mix_n (N);
  }
  
  for (int s=0; s<S; s++) {
    if (debug) Rcout << " s:" << s << endl;
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    if (debug) Rcout << " sample lambda" << endl;
    if ( !normal ) {
      List df_tmp     = sample_df ( aux_df, adaptive_scale, aux_lambda, s, adptive_alpha_gamma );
      aux_df          = as<double>(df_tmp["aux_df"]);
      adaptive_scale  = as<double>(df_tmp["adaptive_scale"]);
      
      U_std           = sum(square( solve(trimatl(aux_L), U) ));
      U_std          /= N * trans(aux_sigma2);
      u               = trans(U_std);
      aux_lambda      = sample_lambda ( aux_df, u , N);
      aux_hetero_inv  = 1 / (aux_sigma2 % aux_lambda);
    }
    
    if (debug) Rcout << " sample sv" << endl;
    List sv_n;
    if (!homoskedastic) {
      
      U_std           = solve(trimatl(aux_L), U);
      U_std.each_row() /= trans(sqrt(aux_lambda));
      u               = trans(sum( log(square(U_std) + ccc) )) / N;
      if ( centred_sv ) {
        sv_n          = svar_ce1( aux_h, aux_rho, aux_omega, aux_sigma2v, aux_sigma2_omega, aux_s_, aux_S, u, prior, aux_mix, true);
      } else {
        sv_n          = svar_nc1( aux_h, aux_rho, aux_omega, aux_sigma2v, aux_sigma2_omega, aux_s_, aux_S, u, prior, aux_mix, true, debug );
      }

      aux_h           = as<vec>(sv_n["aux_h"]);
      aux_S           = as<uvec>(sv_n["aux_S"]);
      aux_rho         = as<double>(sv_n["aux_rho"]);
      aux_omega       = as<double>(sv_n["aux_omega"]);
      aux_sigma2v     = as<double>(sv_n["aux_sigma2v"]);
      aux_sigma2_omega = as<double>(sv_n["aux_sigma2_omega"]);
      aux_s_          = as<double>(sv_n["aux_s_"]);

      if ( centred_sv ) {
        aux_sigma2    = exp(aux_h);
      } else {
        aux_sigma2    = exp(aux_omega * aux_h);
      }
      aux_hetero_inv  = 1 / (aux_sigma2 % aux_lambda);
    }
    
    if (debug) Rcout << " sample ASigma" << endl;
    aux_V                   = sample_V_mgig(aux_V, aux_A, aux_Sigma_inv, prior);
    aux_V_inv               = inv_sympd(aux_V);
    if (debug) Rcout << " aux_hetero_inv" << min(aux_hetero_inv) << endl;
    field<mat> aux_ASigma = sample_ASigma( Y, X, aux_V_inv, aux_hetero_inv, prior );
    aux_A                 = aux_ASigma(0);
    aux_Sigma             = aux_ASigma(1);
    aux_Sigma_inv         = inv_sympd(aux_Sigma);
    aux_L                 = chol(aux_Sigma, "lower");
    U                     = Y - aux_A * X;
    
    if (debug) Rcout << " store post " << endl;
    if (s % thin == 0) {
      posterior_A.slice(ss)     = aux_A;
      posterior_Sigma.slice(ss) = aux_Sigma;
      posterior_V.slice(ss)     = aux_V;
      posterior_sigma2.col(ss)  = aux_sigma2;
      posterior_h.col(ss)       = aux_h;
      posterior_S.col(ss)       = aux_S;
      posterior_rho(ss)         = aux_rho;
      posterior_omega(ss)       = aux_omega;
      posterior_sigma2v(ss)     = aux_sigma2v;
      posterior_sigma2_omega(ss) = aux_sigma2_omega;
      posterior_s_(ss)          = aux_s_;
      posterior_lambda.col(ss)  = aux_lambda;
      posterior_df(ss)          = aux_df;
      ss++;
    }
  } // END s loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["A"]        = aux_A,
      _["Sigma"]    = aux_Sigma,
      _["V"]        = aux_V,
      _["sigma2"]   = aux_sigma2,
      _["h"]        = aux_h,
      _["S"]        = aux_S,
      _["rho"]      = aux_rho,
      _["omega"]    = aux_omega,
      _["sigma2v"]  = aux_sigma2v,
      _["sigma2_omega"] = aux_sigma2_omega,
      _["s_"]       = aux_s_,
      _["lambda"]   = aux_lambda,
      _["df"]       = aux_df
    ),
    _["posterior"]  = List::create(
      _["A"]        = posterior_A,
      _["Sigma"]    = posterior_Sigma,
      _["V"]        = posterior_V,
      _["sigma2"]   = posterior_sigma2,
      _["h"]        = posterior_h,
      _["S"]        = posterior_S,
      _["rho"]      = posterior_rho,
      _["omega"]    = posterior_omega,
      _["sigma2v"]  = posterior_sigma2v,
      _["sigma2_omega"] = posterior_sigma2_omega,
      _["s_"]       = posterior_s_,
      _["lambda"]   = posterior_lambda,
      _["df"]       = posterior_df
    )
  );
} // END bvar_mgig_cpp
