#ifndef _BVAR_GIG_

#define _BVAR_GIG_

#include <RcppArmadillo.h>

Rcpp::List bvar_gig_cpp(
    const int&  S,                        // number of draws from the posterior
    const arma::mat&  Y,                  // NxT dependent variables
    const arma::mat&  X,                  // KxT dependent variables
    const Rcpp::List& prior,              // a list of priors
    const Rcpp::List& starting_values,    // a list of starting values
    const int         thin = 100,         // introduce thinning
    const bool        show_progress = true
);

#endif  // _BVAR_GIG_
