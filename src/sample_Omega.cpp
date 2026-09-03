
#include <RcppArmadillo.h>
#include <bsvars.h>
#include <RcppTN.h>

using namespace Rcpp;
using namespace arma;



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec sample_lambda (
    double&     aux_df,
    arma::vec&  U,
    const int   N
) {
  
  const int T         = U.n_elem;
  U                  /= accu(U) / T;        // normalisation E[u] = 1
  double  nu_lambda   = aux_df + N;
  vec     s_lambda    = N * U + aux_df - 2;
  vec     aux_lambda  = chi2rnd(nu_lambda, T);
  aux_lambda          = s_lambda / aux_lambda;
  
  return aux_lambda;
} // END sample_lambda


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double log_kernel_df (
    const double&         aux_df,
    const arma::vec&   aux_lambda  // Tx1
) {
  
  const int T   = aux_lambda.n_elem;
  double lk_df  = 0;
  lk_df   -= T * lgamma(0.5 * aux_df);                        // lambda prior
  lk_df   += 0.5 * T * aux_df * log(0.5 * (aux_df - 2));      // lambda prior
  lk_df   -= 0.5 * (aux_df + 2) * accu(log(aux_lambda));      // lambda prior
  lk_df   -= 0.5 * (aux_df - 2) * accu(pow(aux_lambda, -1));  // lambda prior
  lk_df   -= 2 * log(aux_df - 1);                             // df prior
  
  return lk_df;
} // END log_kernel_df


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List sample_df (
    double&           aux_df,             // Nx1
    double&           adaptive_scale,     // Nx1
    const arma::vec&  aux_lambda,         // NxT
    const int&        s,                  // MCMC iteration
    const arma::vec&  adptive_alpha_gamma // 2x1 vector with target acceptance rate and step size
) {
  
  double aux_df_star;
  double alpha = 1;
  
  // by sampling from truncated normal it is assumed that the asymmetry from truncation
  // is negligible for alpha computation
  aux_df_star           = RcppTN::rtn1( aux_df, adaptive_scale, 2, R_PosInf );
  double lk_nu_star     = log_kernel_df(aux_df_star, aux_lambda);
  double lk_nu_old      = log_kernel_df(aux_df, aux_lambda);
  double cgd_ratio      = RcppTN::dtn1( aux_df_star, aux_df, adaptive_scale, 2, R_PosInf ) /
    RcppTN::dtn1( aux_df, aux_df_star, adaptive_scale, 2, R_PosInf );
  
  double kernel_ratio   = exp(lk_nu_star - lk_nu_old) * cgd_ratio;
  if ( kernel_ratio < 1 ) alpha = kernel_ratio;
  if ( R::runif(0, 1) < alpha ) {
    aux_df              = aux_df_star;
  }
  
  if (s > 1) {
    adaptive_scale      = exp( log(adaptive_scale) + 0.5 * log( 1 + pow(s, - adptive_alpha_gamma(1)) * (alpha - adptive_alpha_gamma(0))) );
  }
  
  return List::create(
    _["aux_df"] = aux_df,
    _["adaptive_scale"] = adaptive_scale
  );
} // END sample_df


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sv_aux_mix_n (
    const int N
) {
  
  // Computes the probabilities, means, and variances of the 10-component 
  // auxiliary mixture to approximate the log-chi-squared distribution with N 
  // degrees of freedom, adapting th one by Omori et al. (2007)
  double  NN = N;
  mat     data(1, 1e6);
  for (int n=0; n<N; n++) {
    rowvec nn(1e6, fill::randn);
    data += log(square(nn)) / NN;
  }
  
  gmm_diag model;
  model.learn(data, 10, maha_dist, random_subset, 50, 50, 1e-10, false);
  mat out(3, 10);
  out.row(0) = model.hefts;
  out.row(1) = model.means;
  out.row(2) = model.dcovs;
  return out;
} // END sv_aux_mix








/*______________________function find_mixture_indicator_cdf______________________*/
// utility function from file utils_latent_states.cc from the source code of package stochvol
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec find_mixture_indicator_cdf (
    const arma::vec&    datanorm,           // provide all that is conditionally normal
    const arma::mat&    aux_mix             // 3x10 matrix with the parameters of the auxiliary mixture rows: 1-probs, 2-means, 3-vars
) {
  
  const int T = datanorm.n_elem;
  vec mixprob(10 * T);
  for (int j = 0; j < T; j++) {  // TODO slow (10*T calls to exp)!
    const int first_index = 10*j;
    mixprob(first_index) = std::exp(aux_mix(0,0) - (datanorm(j) - aux_mix(1,0)) * (datanorm(j) - aux_mix(1,0)) / aux_mix(2,0) );
    for (int r = 1; r < 10; r++) {
      mixprob(first_index+r) = mixprob(first_index+r-1) + std::exp(aux_mix(0,r) - (datanorm(j) - aux_mix(1,r)) * (datanorm(j) - aux_mix(1,r)) / aux_mix(2,r) );
    }
  }
  return mixprob;
}



/*______________________function svar_nc1______________________*/
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List svar_nc1 (
    arma::vec&        aux_h,
    double&           aux_rho,
    double&           aux_omega,
    double&           aux_sigma2v,
    double&           aux_sigma2_omega, // omega prior hyper-parameter 
    double&           aux_s_,           // scale of IG2 prior for aux_sigma2_omega_n
    arma::uvec&       aux_S,
    const arma::vec&  u,
    const Rcpp::List& prior,
    const arma::mat&  aux_mix,          // 3x10 matrix with the parameters of the auxiliary mixture rows: 1-probs, 2-means, 3-vars
    bool              sample_s_ = true,
    bool              debug = false
) {
  
  if (debug) Rcout << " sv init" << endl;
  // sampler for the non-centred parameterisation of the SV process
  // const double        ccc     = 0.000000001;      // a constant to make log((u+ccc)^2) feasible
  
  // sample h and omega of the non-centered SV including ASIS step
  if (debug) Rcout << " sample h and omega" << endl;
  const int     T = u.n_rows;
  const vec     U = u;
  
  const double  prior_sv_a_ = prior["sv_a"];
  const double  prior_sv_s_ = prior["sv_s"];
  
  mat           H_rho(T, T, fill::eye);
  H_rho.diag(-1)       -= aux_rho;
  mat           HH_rho  = H_rho.t() * H_rho;
  
  // sample auxiliary mixture states aux_S
  if (debug) Rcout << " sample S" << endl;
  const vec   mixprob   = find_mixture_indicator_cdf(U - aux_omega * aux_h, aux_mix);
  aux_S                 = bsvars::inverse_transform_sampling(mixprob, T);
  
  rowvec    alpha_S(T);
  rowvec    sigma_S_inv(T);
  for (int t=0; t<T; t++) {
    alpha_S.col(t)      = aux_mix(1,aux_S(t));
    sigma_S_inv.col(t)  = 1/aux_mix(2,aux_S(t));
  }
  
  // sample aux_s_n
  if (debug) Rcout << " sample s_" << endl;
  if ( sample_s_ ) {
    aux_s_               = (prior_sv_s_ + 2 * aux_sigma2_omega)/chi2rnd(3 + 2 * prior_sv_a_);
  }
  
  // sample aux_sigma2_omega
  if (debug) Rcout << " sample sigma2_omega" << endl;
  aux_sigma2_omega      = bsvars::do_rgig1( prior_sv_a_-0.5, pow(aux_omega,2), 2/aux_s_ );
  
  // sample aux_rho
  if (debug) Rcout << " sample rho" << endl;
  vec       hm1         = aux_h.subvec(0,T-2);
  double    aux_rho_var = as_scalar(pow(hm1.t()*hm1, -1));
  double    aux_rho_mean = as_scalar(aux_rho_var * hm1.t() * aux_h.subvec(1,T-1));
  double    upper_bound = pow(1 - aux_sigma2_omega, 0.5);
  aux_rho               = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5),-upper_bound,upper_bound);
  
  mat       H_rho_new(T, T, fill::eye);
  H_rho_new.diag(-1)   -= aux_rho;
  H_rho                 = H_rho_new;
  HH_rho                = H_rho_new.t() * H_rho_new;
  
  // sample aux_omega
  if (debug) Rcout << " sample omega" << endl;
  double    V_omega_inv = 1/( as_scalar(aux_h.t() * diagmat(sigma_S_inv) * aux_h) + pow(aux_sigma2_omega, -1) );
  double    omega_bar   = as_scalar(aux_h.t() * diagmat(sigma_S_inv) * (U - alpha_S.t()));
  double    omega_aux   = randn( distr_param(V_omega_inv*omega_bar, sqrt(V_omega_inv) ));
  
  // sample aux_h
  if (debug) Rcout << " sample h" << endl;
  mat       V_h         = pow(omega_aux, 2) * diagmat(sigma_S_inv) + HH_rho;
  vec       h_bar       = omega_aux * diagmat(sigma_S_inv) * (U - alpha_S.t());
  vec       h_aux       = bsvars::precision_sampler_ar1( V_h.diag(), V_h(1, 0), h_bar);
  
  // ASIS
  if (debug) Rcout << " ASIS" << endl;
  vec       aux_h_tilde = omega_aux * h_aux;
  double    hHHh        = as_scalar( aux_h_tilde.t() * HH_rho * aux_h_tilde );
  aux_sigma2v           = bsvars::do_rgig1( -0.5*(T-1), hHHh, 1/aux_sigma2_omega );
  int       ss=1;
  if (R::runif(0,1)<0.5) ss *= -1;
  aux_omega             = ss * sqrt(aux_sigma2v);
  aux_h                 = aux_h_tilde / aux_omega;
  
  // ASIS: resample aux_rho
  hm1                   = aux_h.subvec(0,T-2);
  aux_rho_var           = as_scalar(pow(hm1.t() * hm1, -1));
  aux_rho_mean          = as_scalar(aux_rho_var * hm1.t() * aux_h.subvec(1,T-1));
  upper_bound           = pow(1 - aux_sigma2_omega, 0.5);
  aux_rho               = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5),-upper_bound,upper_bound);
  
  return List::create(
    _["aux_h"]              = aux_h,
    _["aux_rho"]            = aux_rho,
    _["aux_omega"]          = aux_omega,
    _["aux_sigma2v"]        = aux_sigma2v,
    _["aux_sigma2_omega"]   = aux_sigma2_omega,
    _["aux_s_"]             = aux_s_,
    _["aux_S"]              = aux_S
  );
} // END sv_nc1



/*______________________function svar_ce1______________________*/
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List svar_ce1 (
    arma::vec&          aux_h,
    double&             aux_rho,
    double&             aux_omega,
    double&             aux_sigma2v,
    double&             aux_sigma2_omega,   // omega prior hyper-parameter 
    double&             aux_s_,             // scale of IG2 prior for aux_sigma2_omega_n
    arma::uvec&         aux_S,
    const arma::vec&    u,
    const Rcpp::List&   prior,
    const arma::mat&    aux_mix,            // 3x10 matrix with the parameters of the auxiliary mixture rows: 1-probs, 2-means, 3-vars
    bool                sample_s_ = true
) {
  // sampler for the centred parameterisation of the SV process
  // const double        ccc     = 0.000000001;      // a constant to make log((u+ccc)^2) feasible
  
  // sample h and omega of the non-centered SV including ASIS step
  const int     T = u.n_elem;
  const vec     U = u;
  
  const double  prior_sv_a_ = prior["sv_a"];
  const double  prior_sv_s_ = prior["sv_s"];
  
  mat           H_rho(T, T, fill::eye);
  H_rho.diag(-1)       -= aux_rho;
  mat           HH_rho  = H_rho.t() * H_rho;
  
  // sample auxiliary mixture states aux_S
  const vec   mixprob   = find_mixture_indicator_cdf(U - aux_omega * aux_h, aux_mix);
  aux_S                 = bsvars::inverse_transform_sampling(mixprob, T);
  
  rowvec    alpha_S(T);
  rowvec    sigma_S_inv(T);
  for (int t=0; t<T; t++) {
    alpha_S.col(t)      = aux_mix(1,aux_S(t));
    sigma_S_inv.col(t)  = 1/aux_mix(2,aux_S(t));
  }
  
  // sample aux_s_n
  if ( sample_s_ ) {
    aux_s_              = (1 + 2 * aux_sigma2_omega) / chi2rnd(3 + 2 * prior_sv_a_);
  }
  
  // sample aux_sigma2_omega
  aux_sigma2_omega      = randg( distr_param(1 + 0.5 * prior_sv_a_, pow(pow(prior_sv_s_,-1) + pow(2 * aux_sigma2v,-1), -1)  ) );
  
  // sample aux_rho
  vec    hm1            = aux_h.subvec(0,T-2);
  double    aux_rho_var = as_scalar(pow( hm1.t() * hm1 / aux_sigma2v, -1));
  double    aux_rho_mean = as_scalar(aux_rho_var * (hm1.t() * aux_h.subvec(1,T-1) / aux_sigma2v) );
  aux_rho               = RcppTN::rtn1(aux_rho_mean, pow(aux_rho_var, 0.5),-1,1);
  
  mat       H_rho_new(T, T, fill::eye);
  H_rho_new.diag(-1)   -= aux_rho;
  H_rho                 = H_rho_new;
  HH_rho                = H_rho_new.t() * H_rho_new;
  
  // sample aux_sigma2v
  aux_sigma2v           = (aux_sigma2_omega + as_scalar(aux_h.t() * HH_rho * aux_h)) / chi2rnd( 3 + T );
  aux_omega             = pow(aux_sigma2v, 0.5);
  
  // sample aux_h
  mat       V_h         = diagmat(sigma_S_inv) + (HH_rho / aux_sigma2v);
  vec       h_bar       = diagmat(sigma_S_inv) * (U - alpha_S.t());
  aux_h                 = bsvars::precision_sampler_ar1( V_h.diag(), V_h(1, 0), h_bar);
  
  return List::create(
    _["aux_h"]              = aux_h,
    _["aux_rho"]            = aux_rho,
    _["aux_omega"]          = aux_omega,
    _["aux_sigma2v"]        = aux_sigma2v,
    _["aux_sigma2_omega"]   = aux_sigma2_omega,
    _["aux_s_"]             = aux_s_,
    _["aux_S"]              = aux_S
  );
} // END svar_ce1
