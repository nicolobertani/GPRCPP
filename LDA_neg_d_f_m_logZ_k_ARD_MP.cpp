#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

double k_se(const vec &x, const vec &y, const double &m) {
  double df = pow(norm(x - y, 2), 2);
  double res = exp(- df / (2 * pow(m, 2)));
  return res;
}

mat K_se(const mat &M, const mat &N, const double &m, const bool &equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;
  omp_set_dynamic(0);

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = k_se(M.row(r).t(), M.row(c).t(), m);
      }
    }
    K = K + K.t();
    // fill diag
    #pragma omp parallel for
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = k_se(M.row(i).t(), M.row(i).t(), m);
    }

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = k_se(M.row(r).t(), N.row(c).t(), m);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}

mat k_ARD(const mat &X, const mat &Y, const vec &p_vec, const bool &equal_mx) {
  // initial checks
  if (!equal_mx) {
    if (X.n_cols != Y.n_cols) Rcout << "Unequal number of columns in the matrices.\n";
  }
  int p = p_vec.n_elem;
  // Rcout << p << "\n";
  if (X.n_cols != (p - 1)) Rcout << "Unequal number of matrix columns and parameters.\n";
  // generate variables
  cube K_cube;
  mat K_0;
  mat out;
  if (equal_mx == 1) {
    K_0.set_size(X.n_rows, X.n_rows);
    K_0.ones();
    K_cube.set_size(X.n_rows, X.n_rows, p - 1);
  } else {
    K_0.set_size(X.n_rows, Y.n_rows);
    K_0.ones();
    K_cube.set_size(X.n_rows, Y.n_rows, p - 1);
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
  out = pow(p_vec(p - 1), 2) * K_0;
  return(out);
}

double log_d_pois(const vec &y, const vec &lambda) {
  Function f("dpois");
  NumericVector temp = f(y, lambda, 1);
  vec out(as<vec>(temp));
  return(sum(out));
}

vec d_likelihood_m (const vec &y, const vec &f, const mat &M, const vec &lambda, const int &m_size) {
  vec v_out = zeros<vec>(m_size);
  #pragma omp parallel for schedule(static)
  for (size_t j = 0; j < m_size; j++) {
    v_out(j) = dot(y - lambda, M.col(j));
  }
  return(v_out);
}


// [[Rcpp::export]]
List neg_d_f_m_logZ_m (const vec &par, const vec &y, const vec &f, const mat &X, const mat &Y, const double &jitter, const bool &compute_d) {
  // initialize values
  int m_size = f.n_elem;
  List out(2);
  double approx_logZ;
  // mm part
  sp_mat jitter_mx = zeros<sp_mat>(m_size, m_size);
  jitter_mx.diag().fill(jitter);
  mat K_mm = k_ARD(X, X, par, 1) + jitter_mx;
  mat inv_K_mm = inv_sympd(K_mm);
  // nm part
  mat K_nm = k_ARD(Y, X, par, 0);
  // other matrices and input
  mat M = K_nm * inv_K_mm;
  vec lambda = exp(M * f);
  vec inv_K_mm_f = inv_K_mm * f;
  // LOG-LIKELIHOOD
  approx_logZ = log_d_pois(y, lambda) - .5 * dot(f, inv_K_mm_f); //omitting terms independent of f_m
  out(0) = - approx_logZ;
  if(compute_d) {
    vec d_logZ_f_m;
    d_logZ_f_m = d_likelihood_m(y, f, M, lambda, m_size) - inv_K_mm_f;
    out(1) = - d_logZ_f_m;
  }
  return(out);
}
