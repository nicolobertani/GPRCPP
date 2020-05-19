#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double k_se_diff(double df, const double l, const double s) {
  double res = std::pow(s, 2) * exp(- df / (2 * std::pow(l, 2)));
  return res;
}

// [[Rcpp::export]]
mat rcpp_k_se_diff(const mat &A, const double l, const double s, bool symmetric) {

  mat K(A.n_rows, A.n_cols, fill::zeros);
  if (symmetric) {

    // Rcout << "I am executing TRUE." << "\n"; // progress message
    // fill upper triangular wo diag
    for (int r = 0; r < A.n_rows; r++) {
      for (int c = r + 1; c < A.n_cols; c++) {
        K(r,c) = k_se_diff(A(r,c), l, s);
      }
    }
    // fill lower triangle
    K = K + K.t();
    // fill diag
    for (int i = 0; i < A.n_rows; i++) {
      K(i,i) = k_se_diff(A(i,i), l, s);
    }

  } else {

    // Rcout << "I am executing FALSE." << "\n"; // progress message
    for (int r = 0; r < A.n_rows; r++) {
      for (int c = 0; c < A.n_cols; c++) {
        K(r,c) = k_se_diff(A(r,c), l, s);
      }
    }
  }
  // Rcout << "I ACTUALLY GOT TO THE END!\n"; // progress message
  return K;
}
