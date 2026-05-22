
#ifndef _COMPUTE_VARIANCE_DECOMPOSITIONS_H_
#define _COMPUTE_VARIANCE_DECOMPOSITIONS_H_

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

#endif