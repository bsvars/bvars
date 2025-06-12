
#ifndef _SAMPLE_ASIGMAV_H_
#define _SAMPLE_ASIGMAV_H_

#include <RcppArmadillo.h>


arma::mat sample_V(
    arma::mat&      aux_V,          // (K,K) matrix
    const arma::mat&      aux_A,          // (N,K) matrix
    const arma::mat&      aux_Sigma_inv,  // (N,N) matrix
    const Rcpp::List&     prior           // a list of prior parameters
);


#endif  // _SAMPLE_ASIGMAV_H_