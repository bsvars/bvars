
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat sv_aux_mix_n (
    const int N
) {
  // Computes the probabilities, means, and variances of the 10-component 
  // auxiliary mixture to approximate the log-chi-squared distribution with N 
  // degrees of freedom, adapting th one by Omori et al. (2007)
  
  mat nn(N, N * 1e6, fill::randn);
  mat data = mean(log(square(nn)));
  gmm_diag model;
  model.learn(data, 10, maha_dist, random_subset, 50, 50, 1e-10, false);
  mat out(3, 10);
  out.row(0) = model.hefts;
  out.row(1) = model.means;
  out.row(2) = model.dcovs;
  return out;
} // END sv_aux_mix


// [[Rcpp::export]]
arma::cube sv_aux_mix_go (
    const int N_max
) {
  cube out(3,10,N_max);
  for (int i=0; i<N_max; i++) {
    Rcout << " i: " << i << endl;
    out.slice(i) = sv_aux_mix_n(i + 1);
  }
  return out;
}


/*** R
# @description Parameters of the 10-component auxiliary mixture of normal 
# components model approximating the log-chi-squared distribution with \code{N} 
# degrees of freedom. These components are used for the common Stochastic 
# Volatility estimation and are applied to the log-square-linearised system 
# following a similar approach as in Omori, Chib, Shephard, Nakajima (2007).
#
# @references 
# Omori, Y., Chib, S., Shephard, N., Nakajima, J., (2007). Stochastic 
# volatility with leverage: Fast and efficient likelihood inference, Journal of 
# Econometrics, 140(2), 425-449, DOI: \doi{10.1016/j.jeconom.2006.07.008}.

sv_aux_mix = sv_aux_mix_go(100)
dimnames(sv_aux_mix) <- list(
  c("probability","mean","variance"),
  paste0("component",1:10),
  paste0("N=",1:dim(sv_aux_mix)[3])
)
usethis::use_data(sv_aux_mix, internal = TRUE)
*/

