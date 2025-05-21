
#ifndef _MATNORM_H_
#define _MATNORM_H_

#include <RcppArmadillo.h>


arma::mat do_rmatnorm1(     // vec(X) ~N(vec(M), S %x% V)
    const arma::mat& M,     // mean
    const arma::mat& V,     // column-specific covariance
    const arma::mat& S      // row-specific covariance
);


#endif  // _MATNORM_H_