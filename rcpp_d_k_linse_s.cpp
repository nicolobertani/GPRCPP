#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double d_k_lin_s(const vec &x, const vec &y, const double l, const double s) {
  double res = 2 * s * dot(x - l, y - l);
  return res;
}

double k_se(const vec &x, const vec &y, const double m) {
  double df = pow(norm(x - y, 2), 2);
  double res = exp(- df / (2 * pow(m, 2)));
  return res;
}

mat rcpp_d_k_lin_s(const mat &M, const mat &N, const double l, const double s, const bool equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = d_k_lin_s(M.row(r).t(), M.row(c).t(), l, s);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = d_k_lin_s(M.row(i).t(), M.row(i).t(), l, s);
    }

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_lin_s(M.row(r).t(), N.row(c).t(), l, s);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}

mat rcpp_k_se(const mat &M, const mat &N, const double m, const bool equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = k_se(M.row(r).t(), M.row(c).t(), m);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = k_se(M.row(i).t(), M.row(i).t(), m);
    }

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = k_se(M.row(r).t(), N.row(c).t(), m);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}


// [[Rcpp::export]]
mat rcpp_d_k_linse_s(const mat &M, const mat &N, const mat &P, const mat &Q, const double l, const double s, const double m, const bool equal_matrices) {
  mat d_K_lin_s, K_se, res;
  d_K_lin_s = rcpp_d_k_lin_s(M, N, l, s, equal_matrices);
  K_se = rcpp_k_se(P, Q, m, equal_matrices);
  res = d_K_lin_s % K_se;
  return res;
}
