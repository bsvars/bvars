#include <RcppArmadillo.h>
#include <bsvars.h>
#include "rmgig.h"

using namespace arma;


arma::mat do_rmgig1(
    arma::mat&        current_value,
    const double      lambda,
    const arma::mat&  Psi,
    const arma::mat&  Gamma
) {
  const int N = Psi.n_rows;
  mat B_tilde = chol(current_value, "upper");
  vec a       = square(diagvec(B_tilde));
  mat B       = trans(B_tilde.each_col() / sqrt(a));

  if (N == 1) {
    mat inv_B      = inv(trimatl(B));
    double chi     = as_scalar(inv_B * Gamma * inv_B.t());
    double psi     = as_scalar(B.t() * Psi * B);
    double shape   = lambda + 1;
    a              = bsvars::do_rgig1(shape, chi, psi);
  } else {
    mat inv_B      = inv(trimatl(B));
    vec chi        = diagvec(inv_B * Gamma * inv_B.t());
    vec psi        = diagvec(B.t() * Psi * B);
    vec shape      = lambda + N - regspace<vec>(1, N) + 1;

    for (int n = 0; n < N; n++) {
      a(n) = bsvars::do_rgig1(shape(n), chi(n), psi(n));
    }

    mat P  = inv_B.t() * diagmat(pow(a, -1)) * inv_B;
    mat M1 = Psi;
    mat R1 = B;
    R1.submat(1, 0, N - 1, 0).zeros();
    R1 = R1.t();
    R1 = trans(R1.each_col() % sqrt(a));

    mat M2 = Gamma;
    mat R2 = inv_B.t();
    R2.row(0) += trans(B.submat(1, 0, N - 1, 0)) * R2.rows(1, N - 1);
    R2 = R2.t();
    R2 = trans(R2.each_col() / sqrt(a));

    vec mean_n = -M1.submat(1, 0, N - 1, N - 1) * R1 * trans(R1.row(0));
    mean_n    += R2.rows(1, N - 1) * R2.t() * trans(M2.row(0));
    mat precision_n = a(0) * Psi.submat(1, 1, N - 1, N - 1)
      + M2(0, 0) * P.submat(1, 1, N - 1, N - 1);
    mat chol_precision = chol(precision_n, "upper");
    vec z(N - 1, fill::randn);
    z += solve(trimatl(chol_precision.t()), mean_n);
    B.submat(1, 0, N - 1, 0) = solve(trimatu(chol_precision), z);

    for (int i = 1; i < N - 1; i++) {
      M1.col(i - 1) += M1.cols(i, N - 1) * B.submat(i, i - 1, N - 1, i - 1);
      M1.row(i - 1) += trans(B.submat(i, i - 1, N - 1, i - 1)) * M1.rows(i, N - 1);
      R1.rows(i + 1, N - 1) -= B.submat(i + 1, i, N - 1, i) * R1.row(i);
      M2.cols(i, N - 1) -= M2.col(i - 1) * trans(B.submat(i, i - 1, N - 1, i - 1));
      M2.rows(i, N - 1) -= B.submat(i, i - 1, N - 1, i - 1) * M2.row(i - 1);
      R2.row(i) += trans(B.submat(i + 1, i, N - 1, i)) * R2.rows(i + 1, N - 1);

      mean_n = -M1.rows(i + 1, N - 1) * R1 * trans(R1.row(i));
      mean_n += R2.rows(i + 1, N - 1) * R2.t() * trans(M2.row(i));
      precision_n = a(i) * Psi.submat(i + 1, i + 1, N - 1, N - 1)
        + M2(i, i) * P.submat(i + 1, i + 1, N - 1, N - 1);
      chol_precision = chol(precision_n, "upper");
      vec zz(N - i - 1, fill::randn);
      zz += solve(trimatl(chol_precision.t()), mean_n);
      B.submat(i + 1, i, N - 1, i) = solve(trimatu(chol_precision), zz);
    }
  }

  return B * diagmat(a) * B.t();
} // END do_rmgig1
