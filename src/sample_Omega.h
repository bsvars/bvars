
#ifndef _SAMPLE_OMEGA_H_
#define _SAMPLE_OMEGA_H_

#include <RcppArmadillo.h>


arma::vec sample_lambda (
    double&     aux_df,
    arma::vec&  U,
    const int   N
);


double log_kernel_df (
    const double&         aux_df,
    const arma::vec&   aux_lambda  // Tx1
);


Rcpp::List sample_df (
    double&           aux_df,             // Nx1
    double&           adaptive_scale,     // Nx1
    const arma::vec&  aux_lambda,         // NxT
    const int&        s,                  // MCMC iteration
    const arma::vec&  adptive_alpha_gamma // 2x1 vector with target acceptance rate and step size
);


arma::mat sv_aux_mix_n (
    const int N
);


arma::vec find_mixture_indicator_cdf (
    const arma::vec&    datanorm,           // provide all that is conditionally normal
    const arma::mat&    aux_mix             // 3x10 matrix with the parameters of the auxiliary mixture rows: 1-probs, 2-means, 3-vars
);


Rcpp::List svar_nc1 (
    arma::vec&        aux_h,
    double&           aux_rho,
    double&           aux_omega,
    double&           aux_sigma2v,
    double&           aux_sigma2_omega, // omega prior hyper-parameter 
    double&           aux_s_,           // scale of IG2 prior for aux_sigma2_omega_n
    arma::uvec&       aux_S,
    const arma::vec&  u,
    const Rcpp::List& prior,
    const arma::mat&  aux_mix,          // 3x10 matrix with the parameters of the auxiliary mixture rows: 1-probs, 2-means, 3-vars
    bool              sample_s_ = true,
    bool              debug = false
);


Rcpp::List svar_ce1 (
    arma::vec&          aux_h,
    double&             aux_rho,
    double&             aux_omega,
    double&             aux_sigma2v,
    double&             aux_sigma2_omega,   // omega prior hyper-parameter 
    double&             aux_s_,             // scale of IG2 prior for aux_sigma2_omega_n
    arma::uvec&         aux_S,
    const arma::vec&    u,
    const Rcpp::List&   prior,
    const arma::mat&    aux_mix,            // 3x10 matrix with the parameters of the auxiliary mixture rows: 1-probs, 2-means, 3-vars
    bool                sample_s_ = true
);


#endif