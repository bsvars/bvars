
#ifndef _RMGIG_H_
#define _RMGIG_H_

#include <RcppArmadillo.h>


arma::mat do_rmgig1(
    arma::mat&      aux_Sigma,  // (N x N) matrix
    const double    lambda,     // (scalar) scale parameter  
    const arma::mat Psi,        // (N x N) positive definite matrix
    const arma::mat Gamma       // (N x N) positive definite matrix
);


#endif  // _RMGIG_H_