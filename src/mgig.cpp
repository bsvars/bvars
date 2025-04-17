
#include <RcppArmadillo.h>
#include <bsvars.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::interfaces(cpp,r)]]
// [[Rcpp::export]]
arma::mat do_rmgig1(
    arma::mat&      aux_Sigma,  // (N x N) matrix
    const double    lambda,     // (scalar) scale parameter  
    const arma::mat Psi,        // (N x N) positive definite matrix
    const arma::mat Gamma       // (N x N) positive definite matrix
) {
  
  const int N         = Psi.n_rows;
  mat       B_tilde   = chol( aux_Sigma, "upper" );
  vec       a         = square( diagvec( B_tilde ) );
  mat       B         = trans( B_tilde.each_col() / sqrt( a ) );
  
  if ( N == 1 ) {
    mat     inv_b     = inv( trimatl( B ) );
    double  chi_a     = as_scalar( inv_b * Gamma * inv_b.t() );
    double  psi_a     = as_scalar( B.t() * Psi * B );
    double  lambda_a  = lambda + 1;
    a                 = bsvars::do_rgig1( lambda_a, chi_a, psi_a );
  } else {
    
    mat     inv_B     = inv( trimatl( B ) );
    vec     Chi_a     = diagvec( inv_B * Gamma * inv_B.t() );
    vec     Psi_a     = diagvec( B.t() * Psi * B );
    vec     La_tmp(N, fill::ones);
    vec     Lambda_a  = lambda + N - cumsum(La_tmp) + 1;
    for (int n=0; n<N; n++) {
      a(n)            = bsvars::do_rgig1( Lambda_a(n), Chi_a(n), Psi_a(n) );
    }
    
    // b (i=1)
    mat     P         = inv_B.t() * diagmat( pow(a, -1) ) * inv_B;
    mat     M1        = Psi;
    mat     R1        = B;  
    R1.submat(1, 0, N-1, 0) = zeros<vec>(N-1);
    R1                = R1.t();
    R1                = trans( R1.each_col() % sqrt(a) );
    mat     M2        = Gamma;
    mat     R2        = inv_B.t();
    R2.row(0)         = R2.row(0) + trans(B.submat(1, 0, N-1, 0)) * R2.rows(1, N-1);
    R2                = R2.t();
    R2                = trans( R2.each_col() / sqrt(a) );
    vec     vn        = (-1) * M1.submat(1, 0, N-1, N-1) * R1 * trans(R1.row(0));
    vn               += R2.rows(1, N-1) * R2.t() * trans(M2.row(0));
    mat     NN        = a(0) * Psi.submat(1, 1, N-1, N-1) + M2(0, 0) * P.submat(1, 1, N-1, N-1);
    mat     chol_N    = chol(NN, "upper");
    vec     z(N-1, fill::randn);
    z                += solve( trimatl(chol_N.t()), vn);
    B.submat(1, 0, N-1, 0) = solve( trimatu(chol_N), z);
    
    // b (i>1))
    if ( N > 2 ) {
      for (int i=1; i<N-1; i++) {
        M1.col(i-1)       += M1.cols(i, N-1) * B.submat(i, i-1, N-1, i - 1);
        M1.row(i-1)       += trans( B.submat(i, i-1, N-1, i-1) ) * M1.rows(i, N-1);
        R1.rows(i+1, N-1) -= B.submat(i+1, i, N-1, i) * R1.row(i);
        M2.cols(i, N-1)   -= M2.col(i-1) * trans(B.submat(i, i-1, N-1, i-1));
        M2.rows(i, N-1) -= B.submat(i, i-1, N-1, i-1) * M2.row(i);
        R2.row(i)         += trans(B.submat(i+1, i, N-1, i)) * R2.rows(i+1, N-1);
        vn                 = (-1) * M1.rows(i+1, N-1) * R1 * trans(R1.row(i));
        vn                += R2.rows(i+1, N-1) * R2.t() * trans(M2.row(i));
        NN                 = a(i) * Psi.submat(i+1, i+1, N-1, N-1) + M2(i, i) * P.submat(i+1, i+1, N-1, N-1);
        chol_N             = chol(NN, "upper");
        vec     zz(N-i-1, fill::randn);
        zz                += solve( trimatl(chol_N.t()), vn);
        B.submat(i+1, i, N-1, i) = solve( trimatu(chol_N), zz);
      } // END i loop
    } // END if
  } // END else
  
  mat Sigma = B * diagmat(a) * B.t();
  return Sigma;
} // END do_rmgig1
