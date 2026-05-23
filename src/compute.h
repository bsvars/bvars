
#ifndef _COMPUTE_H_
#define _COMPUTE_H_

#include <RcppArmadillo.h>


arma::cube Sigma2B (
    arma::cube&   posterior_Sigma_c    // (N, N, S)
);


arma::field<arma::cube> compute_variance_decompositions (
    arma::cube&   posterior_Sigma,    // (S)(N, N, S)
    arma::cube&   posterior_A,        // (S)(K, N, S)
    const int     horizon,
    const int     p
);


arma::cube compute_fitted_values (
    arma::cube&     posterior_A,        // NxKxS
    arma::cube&     posterior_Sigma,    // NxNxS
    arma::mat&      posterior_sigma2,   // TxS
    arma::mat&      X                   // KxT
);


arma::cube compute_shocks (
    arma::cube&     posterior_A,    // (N, K, S)
    arma::mat&      Y,              // NxT dependent variables
    arma::mat&      X               // KxT dependent variables
);


#endif