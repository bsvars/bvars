#ifndef _NICERIG_

#define _NICERIG_

#include <RcppArmadillo.h>

arma::vec nicerig2 (
    const int     n,    // a positive integer - number of draws
    const double  s,    // a positive scale parameter
    const double  nu    // a positive shape parameter
);
#endif  // _NICERIG_
 