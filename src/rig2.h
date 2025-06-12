#ifndef _RIG_

#define _RIG_

#include <RcppArmadillo.h>

arma::vec rig2 (
    const int     n,    // a positive integer - number of draws
    const double  s,    // a positive scale parameter
    const double  nu    // a positive shape parameter
);

#endif  // _RIG_
 