#ifndef _BVAR_GIG_
#define _BVAR_GIG_

#include <RcppArmadillo.h>

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
);

#endif  // _BVAR_GIG_
