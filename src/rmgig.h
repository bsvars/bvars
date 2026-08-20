#ifndef _RMGIG_H_
#define _RMGIG_H_

#include <RcppArmadillo.h>


arma::mat do_rmgig1(
    arma::mat&        current_value,
    const double      lambda,
    const arma::mat&  Psi,
    const arma::mat&  Gamma
);


#endif  // _RMGIG_H_
