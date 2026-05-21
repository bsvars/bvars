
#ifndef _SAMPLE_ASIGMAV_H_
#define _SAMPLE_ASIGMAV_H_

#include <RcppArmadillo.h>


arma::field<arma::mat> sample_ASigma(
    const arma::mat&      Y,              // (N,T) matrix
    const arma::mat&      X,              // (K,T) matrix
    arma::mat&            aux_V_inv,      // (K,K) matrix
    arma::vec&            aux_Omega_diag_inv, // (T) matrix
    const Rcpp::List&     prior           // a list of prior parameters
);

#endif  // _SAMPLE_ASIGMAV_H_