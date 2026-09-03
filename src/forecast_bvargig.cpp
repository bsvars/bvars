
#include <RcppArmadillo.h>
#include <bsvars.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List forecast_bvarGIG (
    arma::cube&   posterior_Sigma,    // (N, N, S)
    arma::cube&   posterior_A,        // (N, K, S)
    arma::mat&    forecast_sigma2,    // (horizon, S)
    arma::vec&    X_T,                // (K)
    arma::mat&    exogenous_forecast, // (horizon, d)
    arma::mat&    cond_forecast,     // (horizon, N)
    const int&    horizon
) {
  
  const int   N = posterior_Sigma.n_rows;
  const int   S = posterior_Sigma.n_slices;
  const int   K = posterior_A.n_cols;
  const int   d = exogenous_forecast.n_cols;
  
  bool        do_exog = exogenous_forecast.is_finite();
  vec         x_t;
  if ( do_exog ) {
    x_t       = X_T.rows(0, K - 1 - d);
  } else {
    x_t       = X_T.rows(0, K - 1);
  } // END if do_exog
  
  vec         Xt(K);
  cube        forecasts(N, horizon, S);
  cube        out_forecast_mean(N, horizon, S);
  cube        SigmaT(N, N, horizon);
  field<cube> out_forecast_cov(S);
  
  for (int s=0; s<S; s++) {
    
    if ( do_exog ) {
      Xt          = join_cols(x_t, trans(exogenous_forecast.row(0)));
    } else {
      Xt          = x_t;
    } // END if do_exog
    
    for (int h=0; h<horizon; h++) {
      
      vec   cond_forecast_h   = trans(cond_forecast.row(h));
      uvec  nonf_el           = find_nonfinite( cond_forecast_h );
      int   nonf_no           = nonf_el.n_elem;
      
      out_forecast_mean.slice(s).col(h) = posterior_A.slice(s) * Xt;
      SigmaT.slice(h)         = forecast_sigma2(h, s) * posterior_Sigma.slice(s);
      
      if ( nonf_no == N ) {
        forecasts.slice(s).col(h) = mvnrnd( out_forecast_mean.slice(s).col(h), SigmaT.slice(h) );
      } else {
        forecasts.slice(s).col(h) = bsvars::mvnrnd_cond( cond_forecast_h, out_forecast_mean.slice(s), SigmaT.slice(h) );   // does not work if cond_fc_h is all nan
      } // END if nonf_no
      
      if ( h != horizon - 1 ) {
        if ( do_exog ) {
          Xt                  = join_cols( forecasts.slice(s).col(h), Xt.subvec(N, K - 1 - d), trans(exogenous_forecast.row(h + 1)) );
        } else {
          Xt                  = join_cols( forecasts.slice(s).col(h), Xt.subvec(N, K - 1) );
        }
      } // END if h
      
      out_forecast_cov(s)     = SigmaT;
      
    } // END h loop
  } // END s loop
  
  return List::create(
    _["forecasts"]       = forecasts,
    _["forecast_mean"]  = out_forecast_mean,
    _["forecast_cov"]   = out_forecast_cov
  );
} // END forecast_bvarGIG


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat forecast_sigma2_sv1 (
    arma::vec&    posterior_h_T,      // S
    arma::vec&    posterior_rho,      // S
    arma::vec&    posterior_omega,    // S
    const int&    horizon
) {
  
  const int S = posterior_rho.n_elem;
  mat       forecasts_sigma2(horizon, S);
  vec       ht = posterior_h_T;
  
  for (int h=0; h<horizon; h++) {
    ht = posterior_rho % ht + posterior_omega % randn(S);
    forecasts_sigma2.row(h) = trans(exp(ht));
  } // END h loop
  
  return forecasts_sigma2;
} // END forecast_sigma2_sv1



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat forecast_lambda_t1 (
    arma::vec&    posterior_df,
    const int&    horizon
) {
  
  const int       S = posterior_df.n_elem;
  mat             forecasts_lambda(horizon, S);
  forecasts_lambda.each_row() += trans(posterior_df) - 2;
  for (int h=0; h<horizon; h++) {
    forecasts_lambda.row(h)   /= trans(chi2rnd( posterior_df ));
  }
  
  return forecasts_lambda;
} // END forecast_lambda_t1
