
#ifndef _FORECAST_BVARGIG_H_
#define _FORECAST_BVARGIG_H_

#include <RcppArmadillo.h>

Rcpp::List forecast_bvarGIG (
    arma::cube&   posterior_Sigma,    // (N, N, S)
    arma::cube&   posterior_A,        // (N, K, S)
    arma::cube&   forecast_sigma2,    // (N, horizon, S)
    arma::vec&    X_T,                // (K)
    arma::mat&    exogenous_forecast, // (horizon, d)
    arma::mat&    cond_forecast,     // (horizon, N)
    const int&    horizon
);

#endif  // _FORECAST_BVARGIG_H_
