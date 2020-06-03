#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double d_k_per_l(const vec &x, const vec &y, const double l, const double s, const int p) {
  double df = norm(x - y, 2);
  double res = std::pow(s, 2) * exp(- 2 * std::pow(sin(datum::pi * df / p), 2) / std::pow(l, 2)) * 4 * std::pow(sin(datum::pi * df / p), 2) / std::pow(l, 3);
  return res;
}

double d_k_per_s(const vec &x, const vec &y, const double l, const double s, const int p) {
  double df = norm(x - y, 2);
  double res = 2 * s * exp(- 2 * std::pow(sin(datum::pi * df / p), 2) / std::pow(l, 2));
  return res;
}

mat rcpp_d_k_per_l(const mat &M, const mat &N, const double l, const double s, const int p, bool equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = d_k_per_l(M.row(r).t(), M.row(c).t(), l, s, p);
      }
    }
    K = K + K.t();
    // fill diag
    /*for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = d_k_per_l(M.row(i).t(), M.row(i).t(), l, s, p);
    }*/

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_per_l(M.row(r).t(), N.row(c).t(), l, s, p);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}

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

mat d_sum_Lambda_fun(const mat &M1, const mat &M1D, const mat &M2, const mat &M2D, const double sigma_sq) {
  mat D;
  D = M1D * M2  * M2.t()  * M1 +
      M1  * M2D * M2.t()  * M1 + // problem here
      M1  * M2  * M2D.t() * M1 +
      M1  * M2  * M2.t()  * M1D;
  D = D / sigma_sq;
  return D;
}

// [[Rcpp::export]]
double rcpp_step_k_per_svi(const mat &M, const mat &N,
  const double l, const double s, const double sigma_sq, const int p,
  const mat &inv_K_mm, const mat &K_mn, const vec &y, const vec &m, const mat &S, const int par_id) {
    // Rcout << "I am running.\n"; // progress message
    mat d_K_mm(M.n_rows, M.n_rows, fill::zeros);
    mat d_inv_K_mm(M.n_rows, M.n_rows, fill::zeros);
    mat d_K_mn(M.n_rows, N.n_rows, fill::zeros);
    vec diag_d_K_nn(N.n_rows);
    double tr_d_K_nn;
    mat d_Q_nn(N.n_rows, N.n_rows, fill::zeros);
    mat d_sum_Lambda(M.n_rows, M.n_rows, fill::zeros);
    vec S_d_sum_Lambda(m.n_elem);
    double output;

    if (par_id == 1) {

      d_K_mm = rcpp_d_k_per_l(M, N, l, s, p, 1);
      d_K_mn = rcpp_d_k_per_l(M, N, l, s, p, 0);
      for (size_t i = 0; i < N.n_rows; i++) {
        diag_d_K_nn(i) = d_k_per_l(N.row(i).t(), N.row(i).t(), l, s, p);
      }
    } else {

      d_K_mm = rcpp_d_k_per_s(M, N, l, s, p, 1);
      d_K_mn = rcpp_d_k_per_s(M, N, l, s, p, 0);
      for (size_t i = 0; i < N.n_rows; i++) {
        diag_d_K_nn(i) = d_k_per_s(N.row(i).t(), N.row(i).t(), l, s, p);
      }
    }

    d_inv_K_mm = - inv_K_mm * d_K_mm * inv_K_mm;
    tr_d_K_nn = sum(diag_d_K_nn);
    d_Q_nn = (d_K_mn.t() * inv_K_mm * K_mn) + (K_mn.t() * d_inv_K_mm * K_mn) + (K_mn.t() * inv_K_mm * d_K_mn);
    d_sum_Lambda = d_sum_Lambda_fun(inv_K_mm, d_inv_K_mm, K_mn, d_K_mn, sigma_sq);
    for (size_t i = 0; i < m.n_elem; i++) {
      S_d_sum_Lambda(i) = dot(S.row(i), d_sum_Lambda.col(i));
    }

    output = as_scalar(
      (y - K_mn.t() * (inv_K_mm * m)).t() * (d_K_mn.t() * (inv_K_mm * m) + K_mn.t() * (d_inv_K_mm * m)) / sigma_sq -
      (tr_d_K_nn - trace(d_Q_nn)) / (2 * sigma_sq) -
      sum(S_d_sum_Lambda) / 2 +
      (trace((inv_K_mm * S * inv_K_mm - inv_K_mm) * d_K_mm) - m.t() * (d_inv_K_mm * m)) / 2
    );
    // Rcout << "Ok, I got to the end.\n"; // progress message
    return output;
}
