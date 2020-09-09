#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double sq2pi = pow(2 * datum::pi, .5);


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


// DERIVATIVES OF THE COVARIANCE FUNCTION --------------------------------------
// single entry functions
double d_k_per_l(const vec &x, const vec &y, const double &l, const double &s, const int &p) {
  double df = norm(x - y, 2);
  double res = std::pow(s, 2) * exp(- 2 * std::pow(sin(datum::pi * df / p), 2) / std::pow(l, 2)) * 4 * std::pow(sin(datum::pi * df / p), 2) / std::pow(l, 3);
  return res;
}

double d_k_per_s(const vec &x, const vec &y, const double &l, const double &s, const int &p) {
  double df = norm(x - y, 2);
  double res = 2 * s * exp(- 2 * std::pow(sin(datum::pi * df / p), 2) / std::pow(l, 2));
  return res;
}

double d_k_se_l(const vec &x, const vec &y, const double &l, const double &s) {
  double df = std::pow(norm(x - y, 2), 2);
  double res = std::pow(s, 2) * df / std::pow(l, 3) * exp(- df / (2 * std::pow(l, 2)));
  return res;
}

double d_k_se_s(const vec &x, const vec &y, const double &l, const double &s) {
  double df = std::pow(norm(x - y, 2), 2);
  double res = 2 * s * exp(- df / (2 * std::pow(l, 2)));
  return res;
}

double d_k_se_m(const vec &x, const vec &y, const double &m) {
  double df = pow(norm(x - y, 2), 2);
  double res = df / pow(m, 3) * exp(- df / (2 * pow(m, 2)));
  return res;
}

double d_k_IM_partial(const double &SUK_i, const double &dSUK_i, const double &SUK_j, const double &dSUK_j,
                      const double &ell, const double &bandwidth, const double &nb) {
  double res = (pow((SUK_i - SUK_j), 2) / bandwidth - (SUK_i - SUK_j) * (dSUK_i - dSUK_j)) / pow(sq2pi * nb * ell, 2);
  return res;
}

// matrix functions
mat d_K_per_l(const mat &M, const mat &N, const double &l, const double &s, const int &p, const bool &equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = d_k_per_l(M.row(r).t(), M.row(c).t(), l, s, p);
      }
    }
    K = K + K.t();

  } else {

    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_per_l(M.row(r).t(), N.row(c).t(), l, s, p);
      }
    }
  }
  return K;
}

mat d_K_per_s(const mat &M, const mat &N, const double &l, const double &s, const int &p, const bool &equal_matrices) {
  mat K;

  if (equal_matrices == 1) {

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

    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_per_s(M.row(r).t(), N.row(c).t(), l, s, p);
      }
    }
  }
  return K;
}

mat d_K_se_l(const mat &M, const mat &N, const double &l, const double &s, const bool &equal_matrices) {
  mat K;

  if (equal_matrices == 1) {

    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = d_k_se_l(M.row(r).t(), M.row(c).t(), l, s);
      }
    }
    K = K + K.t();

  } else {

    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_se_l(M.row(r).t(), N.row(c).t(), l, s);
      }
    }
  }
  return K;
}

mat d_K_se_s(const mat &M, const mat &N, const double &l, const double &s, const bool &equal_matrices) {
  mat K;

  if (equal_matrices == 1) {

    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = d_k_se_s(M.row(r).t(), M.row(c).t(), l, s);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = 2 * s;
    }

  } else {

    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_se_s(M.row(r).t(), N.row(c).t(), l, s);
      }
    }
  }
  return K;
}

mat d_K_se_m(const mat &M, const mat &N, const double &m, const bool &equal_matrices) {
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

mat d_K_IM_partial(const vec &SUK_vec_x, const vec &dSUK_vec_x, const vec &SUK_vec_y, const vec &dSUK_vec_y,
  const double &ell, const double &bandwidth, const double &n, const bool &equal_vec) {
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


// COMPLETE COVARIANCE MATRIX AND DERIVATIVES ----------------------------------
cube K_full(
  const mat &X_time, const mat &Y_time, const vec &par_time, const int &period,
  const mat &X_geo, const mat &Y_geo, const vec &par_geo,
  const mat &X, const mat &Y, const mat &SUK_X, const mat &dSUK_X, const mat &SUK_Y, const mat &dSUK_Y,
                   const vec &l_IM_vec, const vec &b_vec, const vec &n_vec, const vec &l_other_vec,
                   const bool &equal_mx, const bool &compute_d) {
  // create size variables
  int l_IM_size = l_IM_vec.n_elem;
  int b_size = b_vec.n_elem;
  int n_size = n_vec.n_elem;
  int l_other_size = l_other_vec.n_elem;
  int l_tot = l_IM_size + l_other_size; // amplitude is the last parameter
  // initial checks
  if (par_time.n_elem != 2) Rcout << "Wrong number of parameters for time.\n";
  if (par_geo.n_elem != 2) Rcout << "Wrong number of parameters for geography.\n";
  if (l_IM_size != b_size) Rcout << "Unequal number of lenghscales and bandwidth for IM.\n";
  if (l_IM_size != n_size) Rcout << "Unequal number of lenghscales and n for IM.\n";
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
  // initialize covariance derivative terms
  cube d_time_cube;
  cube d_geo_cube;
  cube d_cube;
  cube b_cube;
  cube out;
  // set sizes
  if (equal_mx == 1) {
    K_ARD_0.set_size(X.n_rows, X.n_rows);
    K_ARD_0.ones();
    K_ARD_cube.set_size(X.n_rows, X.n_rows, l_tot - 1);
    if (compute_d) {
      d_time_cube.set_size(X.n_rows, X.n_rows, 2);
      d_geo_cube.set_size(X.n_rows, X.n_rows, 2);
      d_cube.set_size(X.n_rows, X.n_rows, l_tot - 1); // l_tot - 1 because only needed for characteristic length-scales (no amplitude)
      b_cube.set_size(X.n_rows, X.n_rows, b_size); // b_size: one per bandwidth
      out.set_size(X.n_rows, X.n_rows, 1 + 2 + 2 + l_tot + b_size);
    } else {
      out.set_size(X.n_rows, X.n_rows, 1); // 1 for covariance only
    }
  } else {
    K_ARD_0.set_size(X.n_rows, Y.n_rows);
    K_ARD_0.ones();
    K_ARD_cube.set_size(X.n_rows, Y.n_rows, l_tot - 1);
    if (compute_d) {
      d_time_cube.set_size(X.n_rows, Y.n_rows, 2);
      d_geo_cube.set_size(X.n_rows, Y.n_rows, 2);
      d_cube.set_size(X.n_rows, Y.n_rows, l_tot - 1);
      b_cube.set_size(X.n_rows, Y.n_rows, b_size);
      out.set_size(X.n_rows, Y.n_rows, 1 + 2 + 2 + l_tot + b_size);
    } else {
      out.set_size(X.n_rows, Y.n_rows, 1);
    }
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
  out.slice(0) = K_per(X_time, Y_time, par_time(0), par_time(1), period, equal_mx) +
                 K_se(X_geo, Y_geo, par_geo(0), par_geo(1), equal_mx) +
                 pow(l_vec(l_tot - 1), 2) * K_ARD_0;

  // COMPUTE DERIVATIVE MATRICES
  // compute all individual derivative matrices, first for characteristic length-scales l
  if (compute_d) {
    out.slice(1) = d_K_per_l(X_time, Y_time, par_time(0), par_time(1), period, equal_mx);
    out.slice(2) = d_K_per_s(X_time, Y_time, par_time(0), par_time(1), period, equal_mx);
    out.slice(3) = d_K_se_l(X_geo, Y_geo, par_geo(0), par_geo(1), equal_mx);
    out.slice(4) = d_K_se_s(X_geo, Y_geo, par_geo(0), par_geo(1), equal_mx);
    // COMPUTE DERIVATIVES FOR LENGTH-SCALES OF ARD
    for (size_t i = 0; i < l_tot - 1; i++) {
      d_cube.slice(i) = d_K_se_m(X.col(i), Y.col(i), l_vec(i), equal_mx);
    }
    for (size_t i = 0; i < l_tot - 1; i++) {
      mat temp;
      temp.copy_size(K_ARD_0);
      temp.ones();
      for (size_t j = 0; j < l_tot - 1; j++) {
        if (i == j) {
          temp = temp % d_cube.slice(j);
        } else {
          temp = temp % K_ARD_cube.slice(j);
        }
      }
      out.slice(5 + i) = pow(l_vec(l_tot - 1), 2) * temp;
    }
    // Last entry in the outpust list is for the derivative of the amplitude s
    out.slice(4 + l_tot) = 2 * l_vec(l_tot - 1) * K_ARD_0;

    // COMPUTE DERIVATIVES FOR BANDWIDTHS
    for (size_t i = 0; i < b_size; i++) {
      b_cube.slice(i) = d_K_IM_partial(SUK_X.col(i), dSUK_X.col(i), SUK_Y.col(i), dSUK_Y.col(i), l_IM_vec(i), b_vec(i), n_vec(i), equal_mx);
    }
    for (size_t i = 0; i < b_size; i++) {
      out.slice(5 + l_tot + i) = pow(l_vec(l_tot - 1), 2) * K_ARD_0 % b_cube.slice(i);
    }
  }
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

mat d2_likelihood_m (const vec &y, const vec &f, const mat &M, const vec &lambda, const int &m_size) {
  mat mx_out = zeros<mat>(m_size, m_size);
  for (size_t j = 0; j < m_size; j++) {
    for (size_t k = 0; k < m_size; k++) {
      mx_out(j, k) = - dot(lambda, M.col(j) % M.col(k));
    }
  }
  return(mx_out);
}

vec d3_likelihood_m (const vec &y, const vec &f, const mat &M, const vec &lambda, const int &m_size, const int i) {
  vec v_out = zeros<vec>(m_size);
  for (size_t k = 0; k < m_size; k++) {
    v_out(k) = - dot(lambda, M.col(k) % M.col(k) % M.col(i));
  }
  return(v_out);
}

vec d_f_m (const vec &f, const vec &y,
  const mat &K_mm, const mat &W, const mat &M, const vec &lambda, const int &m_size, mat const &d_mx) {
    vec v_out = zeros<vec>(m_size);
    sp_mat I_m = zeros<sp_mat>(m_size, m_size);
    I_m.diag().ones();
    vec d_lik_vec = d_likelihood_m(y, f, M, lambda, m_size);
    v_out = solve(I_m + K_mm * W, d_mx * d_lik_vec);
    return(v_out);
  }


// [[Rcpp::export]]
List d_logZ_m_mixed (
  const vec &l_IM_vec, const vec &b_vec, const vec &n_vec, const vec &l_other_vec,
  const vec &y, const vec &f,
  const mat &X, const mat &SUK_X, const mat &dSUK_X, const mat &Y, const mat &SUK_Y, const mat &dSUK_Y,
  const double &jitter, const bool &compute_d) {
  // create size variables
  int m_size = f.n_elem;
  int l_IM_size = l_IM_vec.n_elem;
  int l_other_size = l_other_vec.n_elem;
  int l_tot = l_IM_size + l_other_size; // amplitude is the last parameter
  // initialize values
  List out(2);
  double logZ;
  // mm part
  cube all_K_mm = K_full(X, X, SUK_X, dSUK_X, SUK_X, dSUK_X, l_IM_vec, b_vec, n_vec, l_other_vec, 1, compute_d); // covariance and derivatives for mm
  sp_mat jitter_mx = zeros<sp_mat>(m_size, m_size);
  jitter_mx.diag().fill(jitter);
  mat K_mm = all_K_mm.slice(0) + jitter_mx;
  mat chol_K_mm = chol(K_mm);
  mat inv_chol_K_mm = inv(trimatu(chol_K_mm));
  mat inv_K_mm = inv_chol_K_mm * inv_chol_K_mm.t();
  // nm part
  cube all_K_nm = K_full(Y, X, SUK_Y, dSUK_Y, SUK_X, dSUK_X, l_IM_vec, b_vec, n_vec, l_other_vec, 0, compute_d); // covariance and derivatives for nm
  mat K_nm = all_K_nm.slice(0);
  // other matrices and input
  mat M = K_nm * inv_K_mm;
  vec lambda = exp(M * f);
  mat W = - d2_likelihood_m(y, f, M, lambda, m_size);
  mat A = inv_K_mm + W;
  mat chol_A = chol(A);
  mat inv_chol_A = inv(trimatu(chol_A));
  mat inv_A = inv_chol_A * inv_chol_A.t();
  vec inv_K_mm_f = inv_K_mm * f;
  // LOG-LIKELIHOOD
  logZ = log_d_pois(y, lambda) - .5 * dot(f, inv_K_mm_f) - sum(log(chol_K_mm.diag())) - sum(log(chol_A.diag()));
  out(0) = logZ;
  // GRADIENT OF LOG-LIKELIHOOD
  if (compute_d) {
    vec gradient = zeros<vec>(l_tot + l_IM_size); // l_IM_size to account for bandwidths
    // calculate last term of gradient : sapply(seq(m.size), function(h) {diag(inv.A) %*% d3.likelihood.m(y.sample, f.m.hat, M, h)})
    vec inv_A_d3_v = zeros<vec>(m_size);
    for (size_t h = 0; h < m_size; h++) {
      inv_A_d3_v(h) = dot(inv_A.diag(), d3_likelihood_m(y, f, M, lambda, m_size, h));
    }
    // calculate gradient entries
    for (size_t p = 1; p <= l_tot + l_IM_size; p++) { // l_IM_size to account for bandwidths
      // diagonal for first trace entry
      vec diag_1 = zeros<vec>(m_size);
      for (size_t i = 0; i < m_size; i++) {
        diag_1(i) = dot(inv_K_mm.row(i), all_K_mm.slice(p).col(i));
      }
      // gradient entry
      gradient(p - 1) = as_scalar(
        (y.t() * all_K_nm.slice(p)) * inv_K_mm_f -
        ((y.t() * K_nm) * inv_K_mm) * (all_K_mm.slice(p) * inv_K_mm_f) -
        (all_K_nm.slice(p) * inv_K_mm_f - K_nm * (inv_K_mm * (all_K_mm.slice(p) * inv_K_mm_f))).t() * exp(K_nm * inv_K_mm_f) +
        .5 * (
          inv_K_mm_f.t() * all_K_mm.slice(p) * inv_K_mm_f -
          sum(diag_1) +
          trace(inv_K_mm * inv_A * inv_K_mm * all_K_mm.slice(p)) -
          d_f_m(f, y, K_mm, W, M, lambda, m_size, all_K_mm.slice(p)).t() * inv_A_d3_v
        )
      );
    }
    out(1) = gradient;
  }
  return(out);
}
