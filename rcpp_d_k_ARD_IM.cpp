#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double sq2pi = pow(2 * datum::pi, .5);

double k_se(const vec &x, const vec &y, const double m) {
  double df = pow(norm(x - y, 2), 2);
  double res = exp(- df / (2 * pow(m, 2)));
  return res;
}

double d_k_se_m(const vec &x, const vec &y, const double m) {
  double df = pow(norm(x - y, 2), 2);
  double res = df / pow(m, 3) * exp(- df / (2 * pow(m, 2)));
  return res;
}

mat K_se(const mat &M, const mat &N, const double m, const bool equal_matrices) {
  mat K;

  if (equal_matrices == 1) {

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

    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = k_se(M.row(r).t(), N.row(c).t(), m);
      }
    }
  }
  return K;
}

mat d_K_se(const mat &M, const mat &N, const double m, const bool equal_matrices) {
  mat K;

  if (equal_matrices == 1) {

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

    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_se_m(M.row(r).t(), N.row(c).t(), m);
      }
    }
  }
  return K;
}

double d_k_IM_partial(const double &SUK_i, const double &dSUK_i, const double &SUK_j, const double &dSUK_j,
                      const double &ell, const double &bandwidth, const double &nb) {
  double res = (pow((SUK_i - SUK_j), 2) / bandwidth - (SUK_i - SUK_j) * (dSUK_i - dSUK_j)) / pow(sq2pi * nb * ell, 2);
  return res;
}

mat d_k_IM_partial(const vec &SUK_vec_x, const vec &dSUK_vec_x, const vec &SUK_vec_y, const vec &dSUK_vec_y,
  const double &ell, const double &bandwidth, const int &n, const bool &equal_vec) {
  mat K;
  double nb = n * bandwidth;

  if (equal_vec == 1) {

    K.set_size(SUK_vec_x.n_elem, SUK_vec_x.n_elem);
    K.zeros();
    // fill upper triangular wo diag
    for (int i = 0; i < SUK_vec_x.n_elem; i++) {
      for (int j = i + 1; j < SUK_vec_x.n_elem; j++) {
        K(i, j) = d_k_IM_partial(SUK_vec_x(i), dSUK_vec_x(i), SUK_vec_x(j), dSUK_vec_x(j), ell, bandwidth, nb);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < SUK_vec_x.n_elem; i++) {
      K(i,i) = d_k_IM_partial(SUK_vec_x(i), dSUK_vec_x(i), SUK_vec_x(i), dSUK_vec_x(i), ell, bandwidth, nb);
    }

  } else {

    K.set_size(SUK_vec_x.n_elem, SUK_vec_y.n_elem);
    // fill everything
    for (int i = 0; i < SUK_vec_x.n_elem; i++) {
      for (int j = 0; j < SUK_vec_y.n_elem; j++) {
        K(i, j) = d_k_IM_partial(SUK_vec_x(i), dSUK_vec_x(i), SUK_vec_y(j), dSUK_vec_y(j), ell, bandwidth, nb);
      }
    }
  }
  return K;
}

mat SUK_to_IM(const mat &SUK, const vec &b_vec, const vec &n_vec) {
  mat IM;
  IM.copy_size(SUK);
  for (size_t i = 0; i < IM.n_rows; i++) {
    for (size_t j = 0; j < IM.n_cols; j++) {
      IM(i,j) = SUK(i,j) / (sq2pi * n_vec(j) * b_vec(j));
    }
  }
  return(IM);
}


// [[Rcpp::export]]
List rcpp_k_ARD_IM(const mat &SUK_X, const mat &dSUK_X, const mat &SUK_Y, const mat &dSUK_Y,
                   const vec &p_vec, const vec &b_vec, const vec &n_vec,
                   const bool equal_mx, const bool compute_ds) {
  // GENERATE X AND Y FROM SUK_X AND SUK_Y
  // These are the Intensity Measures
  mat X = SUK_to_IM(SUK_X, b_vec, n_vec);
  mat Y = SUK_to_IM(SUK_Y, b_vec, n_vec);
  // initial checks
  if (!equal_mx) {
    if (X.n_cols != Y.n_cols) Rcout << "Unequal number of columns in the matrices.\n";
  }
  int p = p_vec.n_elem;
  if (X.n_cols != (p - 1)) Rcout << "Unequal number of matrix columns and parameters.\n";
  // generate variables
  cube K_cube;
  cube d_cube;
  cube b_cube;
  List out(1 + p + b_vec.n_elem);
  mat K_0;
  if (equal_mx == 1) {
    K_0.set_size(X.n_rows, X.n_rows);
    K_0.ones();
    K_cube.set_size(X.n_rows, X.n_rows, p-1);
    d_cube.set_size(X.n_rows, X.n_rows, p-1); // p-1 because only needed for characteristic length-scales
    b_cube.set_size(X.n_rows, X.n_rows, p-1); // p-1 = #b
  } else {
    K_0.set_size(X.n_rows, Y.n_rows);
    K_0.ones();
    K_cube.set_size(X.n_rows, Y.n_rows, p-1);
    d_cube.set_size(X.n_rows, Y.n_rows, p-1); // p-1 because only needed for characteristic length-scales
    b_cube.set_size(X.n_rows, Y.n_rows, p-1); // p-1 = #b
  }
  // COMPUTE COVARIANCE MATRIX
  // populate list of covariance matrices
  for (size_t i = 0; i < (p - 1); i++) {
    K_cube.slice(i) = K_se(X.col(i), Y.col(i), p_vec(i), equal_mx);
  }
  // multiply the matrices in the list
  for (size_t i = 0; i < (p - 1); i++) {
      K_0 = K_0 % K_cube.slice(i);
  }
  out(0) = pow(p_vec(p - 1), 2) * K_0;
  // COMPUTE DERIVATIVE MATRICES
  // compute all individual derivative matrices, first for characteristic length-scales l
  if (compute_ds) {
    // COMPUTE DERIVATIVES FOR LENGTH-SCALES
    for (size_t i = 0; i < p - 1; i++) {
      d_cube.slice(i) = d_K_se(X.col(i), Y.col(i), p_vec(i), equal_mx);
    }
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
      out(i + 1) = pow(p_vec(p - 1), 2) * temp;
    }
    // Last entry in the output list is for the derivative of the amplitude s
    out(p) = 2 * p_vec(p - 1) * K_0;
    // COMPUTE DERIVATIVES FOR BANDWIDTHS
    for (size_t i = 0; i < p - 1; i++) {
      b_cube.slice(i) = d_k_IM_partial(SUK_X.col(i), dSUK_X.col(i), SUK_Y.col(i), dSUK_Y.col(i), p_vec(i), b_vec(i), n_vec(i), equal_mx);
    }
    for (size_t i = 0; i < p - 1; i++) {
      out(p + 1 + i) = pow(p_vec(p - 1), 2) * K_0 % b_cube.slice(i);
    }
  }
  return out;
}
