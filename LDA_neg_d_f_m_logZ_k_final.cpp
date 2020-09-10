#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// COVARIANCE FUNCTIONS --------------------------------------------------------
// single entry functions
double k_per(const vec &x, const vec &y, const double &l, const double &s, const int &p) {
  double df = norm(x - y, 2);
  double res = pow(s, 2) * exp(- 2 * std::pow(sin(datum::pi * df / p), 2) / pow(l, 2));
  return res;
}

double k_se(const vec &x, const vec &y, const double &l, const double &s) {
  double df = std::pow(norm(x - y, 2), 2);
  double res = std::pow(s, 2) * exp(- df / (2 * std::pow(l, 2)));
  return res;
}

double k_se_m(const vec &x, const vec &y, const double &m) {
  double df = pow(norm(x - y, 2), 2);
  double res = exp(- df / (2 * pow(m, 2)));
  return res;
}

// matrix functions
mat K_per(const mat &M, const mat &N, const double &l, const double &s, const int &p, const bool &equal_matrices) {
  mat K;

  if (equal_matrices == 1) {

    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = k_per(M.row(r).t(), M.row(c).t(), l, s, p);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = k_per(M.row(i).t(), M.row(i).t(), l, s, p);
    }

  } else {

    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = k_per(M.row(r).t(), N.row(c).t(), l, s, p);
      }
    }
  }
  return K;
}

mat K_se(const mat &M, const mat &N, const double &l, const double &s, const bool &equal_matrices) {
  mat K;

  if (equal_matrices == 1) {

    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = k_se(M.row(r).t(), M.row(c).t(), l, s);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = k_se(M.row(i).t(), M.row(i).t(), l, s);
    }

  } else {

    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = k_se(M.row(r).t(), N.row(c).t(), l, s);
      }
    }
  }
  return K;
}

mat K_se_m(const mat &M, const mat &N, const double &m, const bool &equal_matrices) {
  mat K;

  if (equal_matrices == 1) {

    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = k_se_m(M.row(r).t(), M.row(c).t(), m);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = k_se_m(M.row(i).t(), M.row(i).t(), m);
    }

  } else {

    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = k_se_m(M.row(r).t(), N.row(c).t(), m);
      }
    }
  }
  return K;
}


// COMPLETE COVARIANCE MATRIX AND DERIVATIVES ----------------------------------
mat K_full(
  const mat &X_time, const mat &Y_time, const vec &par_time, const int &period,
  const mat &X_geo, const mat &Y_geo, const vec &par_geo,
  const mat &X, const mat &Y,
  const vec &l_IM_vec, const vec &l_other_vec,
  const bool &equal_mx
) {
  // create size variables
  int l_IM_size = l_IM_vec.n_elem;
  int l_other_size = l_other_vec.n_elem;
  int l_tot = l_IM_size + l_other_size; // amplitude is the last parameter
  // initial checks
  if (par_time.n_elem != 2) Rcout << "Wrong number of parameters for time.\n";
  if (par_geo.n_elem != 2) Rcout << "Wrong number of parameters for geography.\n";
  if (X.n_cols != (l_tot - 1)) Rcout << "Unequal number of matrix columns and parameters.\n";
  if (!equal_mx) {
    if (X.n_cols != Y.n_cols) Rcout << "Unequal number of columns in the matrices.\n";
  }
  // generate variables
  vec l_vec = join_cols(l_IM_vec, l_other_vec);
  // initialize covariance function terms
  mat K_time;
  mat K_geo;
  mat K_ARD_0;
  cube K_ARD_cube;
  mat out;
  // set sizes
  if (equal_mx == 1) {
    K_ARD_0.set_size(X.n_rows, X.n_rows);
    K_ARD_0.ones();
    K_ARD_cube.set_size(X.n_rows, X.n_rows, l_tot - 1);
  } else {
    K_ARD_0.set_size(X.n_rows, Y.n_rows);
    K_ARD_0.ones();
    K_ARD_cube.set_size(X.n_rows, Y.n_rows, l_tot - 1);
  }
  // COMPUTE COVARIANCE MATRIX
  // populate list of terms for the ARD
  for (size_t i = 0; i < (l_tot - 1); i++) {
    K_ARD_cube.slice(i) = K_se_m(X.col(i), Y.col(i), l_vec(i), equal_mx);
  }
  // multiply the matrices in the list
  for (size_t i = 0; i < (l_tot - 1); i++) {
      K_ARD_0 = K_ARD_0 % K_ARD_cube.slice(i);
  }
  // full covariance matrix
  out = K_per(X_time, Y_time, par_time(0), par_time(1), period, equal_mx) +
        K_se(X_geo, Y_geo, par_geo(0), par_geo(1), equal_mx) +
        pow(l_vec(l_tot - 1), 2) * K_ARD_0;
  return out;
}


// LIKELIHOOD COMPONENTS -------------------------------------------------------
double log_d_pois(const vec &y, const vec &lambda) {
  Function f("dpois");
  NumericVector temp = f(y, lambda, 1);
  vec out(as<vec>(temp));
  return(sum(out));
}

vec d_likelihood_m (const vec &y, const vec &f, const mat &M, const vec &lambda, const int &m_size) {
  vec v_out = zeros<vec>(m_size);
  for (size_t j = 0; j < m_size; j++) {
    v_out(j) = dot(y - lambda, M.col(j));
  }
  return(v_out);
}


// [[Rcpp::export]]
List neg_d_f_m_logZ_m_final (
  const vec &y, const vec &f,
  const mat &X_time, const mat &Y_time, const vec &par_time, const int &period,
  const mat &X_geo, const mat &Y_geo, const vec &par_geo,
  const mat &X, const mat &Y, const vec &l_IM_vec, const vec &l_other_vec,
  const double &jitter, const bool &compute_d
) {
  // initialize values
  int m_size = f.n_elem;
  List out(2);
  double approx_logZ;
  // mm part
  sp_mat jitter_mx = zeros<sp_mat>(m_size, m_size);
  jitter_mx.diag().fill(jitter);
  mat K_mm = K_full(
    X_time, X_time, par_time, period,
    X_geo, X_geo, par_geo,
    X, X, l_IM_vec, l_other_vec,
    1) + jitter_mx;
  mat inv_K_mm = inv_sympd(K_mm);
  // nm part
  mat K_nm = K_full(
    Y_time, X_time, par_time, period,
    Y_geo, X_geo, par_geo,
    Y, X, l_IM_vec, l_other_vec,
    0);
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
