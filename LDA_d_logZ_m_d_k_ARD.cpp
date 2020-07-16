#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double k_se(const vec &x, const vec &y, const double &m) {
  double df = pow(norm(x - y, 2), 2);
  double res = exp(- df / (2 * pow(m, 2)));
  return res;
}

double d_k_se_m(const vec &x, const vec &y, const double &m) {
  double df = pow(norm(x - y, 2), 2);
  double res = df / pow(m, 3) * exp(- df / (2 * pow(m, 2)));
  return res;
}

mat K_se(const mat &M, const mat &N, const double &m, const bool &equal_matrices) {
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

mat d_K_se(const mat &M, const mat &N, const double &m, const bool &equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = d_k_se_m(M.row(r).t(), M.row(c).t(), m);
      }
    }
    K = K + K.t();

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_se_m(M.row(r).t(), N.row(c).t(), m);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}

cube rcpp_k_ARD(const mat &X, const mat &Y, const vec &p_vec, const bool &equal_mx) {
  // initial checks
  if (!equal_mx) {
    if (X.n_cols != Y.n_cols) Rcout << "Unequal number of columns in the matrices.\n";
  }
  int p = p_vec.n_elem;
  // Rcout << p << "\n";
  if (X.n_cols != (p - 1)) Rcout << "Unequal number of matrix columns and parameters.\n";
  // generate variables
  cube K_cube;
  cube d_cube;
  cube out;
  mat K_0;
  if (equal_mx == 1) {
    K_0.set_size(X.n_rows, X.n_rows);
    K_0.ones();
    K_cube.set_size(X.n_rows, X.n_rows, p - 1);
    d_cube.set_size(X.n_rows, X.n_rows, p - 1); // p-1 because only needed for characteristic length-scales
    out.set_size(X.n_rows, X.n_rows, p + 1); // 1 for covariance plus p for each derivative
  } else {
    K_0.set_size(X.n_rows, Y.n_rows);
    K_0.ones();
    K_cube.set_size(X.n_rows, Y.n_rows, p - 1);
    d_cube.set_size(X.n_rows, Y.n_rows, p - 1); // p-1 because only needed for characteristic length-scales
    out.set_size(X.n_rows, Y.n_rows, p + 1); // 1 for covariance plus p for each derivative
  }
  // COMPUTE COVARIANCE MATRIX
  // populate list of covariance matrices
  for (size_t i = 0; i < (p - 1); i++) {
    // Rcout << i << "\n";
    K_cube.slice(i) = K_se(X.col(i), Y.col(i), p_vec(i), equal_mx);
  }
  // multiply the matrices in the list
  for (size_t i = 0; i < (p - 1); i++) {
      K_0 = K_0 % K_cube.slice(i);
  }
  out.slice(0) = pow(p_vec(p - 1), 2) * K_0;
  // COMPUTE DERIVATIVE MATRICES
  // compute all individual derivative matrices, first for characteristic length-scales l
  // Rcout << "I got to the derivatives!\n";
  for (size_t i = 0; i < p - 1; i++) {
    d_cube.slice(i) = d_K_se(X.col(i), Y.col(i), p_vec(i), equal_mx);
  }
  // Rcout << "I got computed d_cube fine!\n";
  for (size_t i = 0; i < p - 1; i++) {
    mat temp;
    temp.copy_size(K_0);
    temp.ones();
    for (size_t j = 0; j < p - 1; j++) {
      if (i == j) {
        temp = temp % d_cube.slice(j);
      } else {
        temp = temp % K_cube.slice(j);
      }
    }
    out.slice(i + 1) = pow(p_vec(p - 1), 2) * temp;
  }
  // Last entry in the outpust list is for the derivative of the amplitude s
  out.slice(p) = 2 * p_vec(p - 1) * K_0;
  return out;
}

vec d_likelihood_m (const vec &y, const vec &f, const mat &M) {
  vec E = exp(M * f);
  int m_size = f.n_elem;
  vec v_out = zeros<vec>(m_size);
  for (size_t j = 0; j < m_size; j++) {
    v_out(j) = dot(y - E, M.col(j));
  }
  return(v_out);
}

mat d2_likelihood_m (const vec &y, const vec &f, const mat &M) {
  vec E = exp(M * f);
  int m_size = f.n_elem;
  mat mx_out = zeros<mat>(m_size, m_size);
  for (size_t j = 0; j < m_size; j++) {
    for (size_t k = 0; k < m_size; k++) {
      mx_out(j, k) = - dot(E, M.col(j) % M.col(k));
    }
  }
  return(mx_out);
}

vec d3_likelihood_m (const vec &y, const vec &f, const mat &M, const int i) {
  vec E = exp(M * f);
  int m_size = f.n_elem;
  vec v_out = zeros<vec>(m_size);
  for (size_t k = 0; k < m_size; k++) {
    v_out(k) = - dot(E, M.col(k) % M.col(k) % M.col(i));
  }
  return(v_out);
}

vec d_f_m (const vec &par_hat, const vec &f, const vec &y,
  const mat &K_mm, const mat &W, const mat &M, mat const &d_mx) {
    int m_size = f.n_elem;
    vec v_out = zeros<vec>(m_size);
    sp_mat I_m = zeros<sp_mat>(m_size, m_size);
    I_m.diag().ones();
    vec d_lik_vec = d_likelihood_m(y, f, M);
    v_out = solve(I_m + K_mm * W, d_mx * d_lik_vec);
    return(v_out);
  }


// [[Rcpp::export]]
mat ARD(const mat &X, const mat &Y, const vec &p_vec, const bool &equal_mx, const int &slice_n) {
  mat out;
  cube temp;
  temp = rcpp_k_ARD(X, Y, p_vec, equal_mx);
  out = temp.slice(slice_n);
  return(out);
}
