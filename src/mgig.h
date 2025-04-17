
#ifndef _MGIG_H_
#define _MGIG_H_

#include <RcppArmadillo.h>


arma::mat do_rmgig1(
    arma::mat&      aux_Sigma,  // (N x N) matrix
    const double    lambda,     // (scalar) scale parameter  
    const arma::mat Psi,        // (N x N) positive definite matrix
    const arma::mat Gamma       // (N x N) positive definite matrix
);


#endif  // _MGIG_H_