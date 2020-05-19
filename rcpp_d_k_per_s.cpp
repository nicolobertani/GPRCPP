#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double d_k_per_s(const vec &x, const vec &y, const double l, const double s, const int p) {
  double df = norm(x - y, 2);
  double res = 2 * s * exp(- 2 * std::pow(sin(datum::pi * df / p), 2) / std::pow(l, 2));
  return res;
}


// [[Rcpp::export]]
mat rcpp_d_k_per_s(const mat &M, const mat &N, const double l, const double s, const int p, bool equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = d_k_per_s(M.row(r).t(), M.row(c).t(), l, s, p);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = 2 * s;
    }

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_per_s(M.row(r).t(), N.row(c).t(), l, s, p);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}
