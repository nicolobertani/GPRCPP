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

double d_k_se_m(const vec &x, const vec &y, const double &m) {
  double df = pow(norm(x - y, 2), 2);
  double res = df / pow(m, 3) * exp(- df / (2 * pow(m, 2)));
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

mat d_K_se(const mat &M, const mat &N, const double &m, const bool &equal_matrices) {
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
        K(r, c) = d_k_se_m(M.row(r).t(), M.row(c).t(), m);
      }
    }
    K = K + K.t();

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_se_m(M.row(r).t(), N.row(c).t(), m);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}

cube k_ARD(const mat &X, const mat &Y, const vec &p_vec, const bool &equal_mx, const bool &compute_d) {
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
    if (compute_d) {
      d_cube.set_size(X.n_rows, X.n_rows, p - 1); // p-1 because only needed for characteristic length-scales
      out.set_size(X.n_rows, X.n_rows, p + 1); // 1 for covariance plus p for each derivative
    } else {
      out.set_size(X.n_rows, X.n_rows, 1); // 1 for covariance plus p for each derivative
    }
  } else {
    K_0.set_size(X.n_rows, Y.n_rows);
    K_0.ones();
    K_cube.set_size(X.n_rows, Y.n_rows, p - 1);
    if (compute_d) {
      d_cube.set_size(X.n_rows, Y.n_rows, p - 1); // p-1 because only needed for characteristic length-scales
      out.set_size(X.n_rows, Y.n_rows, p + 1); // 1 for covariance plus p for each derivative
    } else {
      out.set_size(X.n_rows, Y.n_rows, 1); // 1 for covariance plus p for each derivative
    }
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
  if (compute_d) {
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
  }
  return out;
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

mat d2_likelihood_m (const vec &y, const vec &f, const mat &M, const vec &lambda, const int &m_size) {
  mat mx_out = zeros<mat>(m_size, m_size);
  #pragma omp parallel for collapse(2) schedule(static)
  for (size_t j = 0; j < m_size; j++) {
    for (size_t k = 0; k < m_size; k++) {
      mx_out(j, k) = - dot(lambda, M.col(j) % M.col(k));
    }
  }
  return(mx_out);
}

vec d3_likelihood_m (const vec &y, const vec &f, const mat &M, const vec &lambda, const int &m_size, const int i) {
  vec v_out = zeros<vec>(m_size);
  #pragma omp parallel for schedule(static)
  for (size_t k = 0; k < m_size; k++) {
    v_out(k) = - dot(lambda, M.col(k) % M.col(k) % M.col(i));
  }
  return(v_out);
}

vec d_f_m (const vec &par_hat, const vec &f, const vec &y,
  const mat &K_mm, const mat &W, const mat &M, const vec &lambda, const int &m_size, mat const &d_mx) {
    vec v_out = zeros<vec>(m_size);
    sp_mat I_m = zeros<sp_mat>(m_size, m_size);
    I_m.diag().ones();
    vec d_lik_vec = d_likelihood_m(y, f, M, lambda, m_size);
    v_out = solve(I_m + K_mm * W, d_mx * d_lik_vec);
    return(v_out);
  }


// [[Rcpp::export]]
List d_logZ_m (const vec &par, const vec &y, const vec &f, const mat &X, const mat &Y, const double &jitter, const bool &compute_d) {
  // initialize values
  int m_size = f.n_elem;
  int p_size = par.n_elem;
  List out(2);
  vec grad = zeros<vec>(p_size);
  double logZ;
  // mm part
  cube all_K_mm = k_ARD(X, X, par, 1, compute_d); // covariance and derivatives for mm
  sp_mat jitter_mx = zeros<sp_mat>(m_size, m_size);
  jitter_mx.diag().fill(jitter);
  mat K_mm = all_K_mm.slice(0) + jitter_mx;
  mat chol_K_mm = chol(K_mm);
  mat inv_chol_K_mm = inv(trimatu(chol_K_mm));
  mat inv_K_mm = inv_chol_K_mm * inv_chol_K_mm.t();
  // nm part
  cube all_K_nm = k_ARD(Y, X, par, 0, compute_d); // covariance and derivatives for nm
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
  // calculate last term of gradient : sapply(seq(m.size), function(h) {diag(inv.A) %*% d3.likelihood.m(y.sample, f.m.hat, M, h)})
  if (compute_d) {
    vec inv_A_d3_v = zeros<vec>(m_size);
    for (size_t h = 0; h < m_size; h++) {
      inv_A_d3_v(h) = dot(inv_A.diag(), d3_likelihood_m(y, f, M, lambda, m_size, h));
    }
    // calculate gradient entries
    for (size_t p = 1; p <= p_size; p++) {
      // diagonal for first trace entry
      vec diag_1 = zeros<vec>(m_size);
      for (size_t i = 0; i < m_size; i++) {
        diag_1(i) = dot(inv_K_mm.row(i), all_K_mm.slice(p).col(i));
      }
      // gradient entry
      grad(p - 1) = as_scalar(
        (y.t() * all_K_nm.slice(p)) * inv_K_mm_f -
        ((y.t() * K_nm) * inv_K_mm) * (all_K_mm.slice(p) * inv_K_mm_f) -
        (all_K_nm.slice(p) * inv_K_mm_f - K_nm * (inv_K_mm * (all_K_mm.slice(p) * inv_K_mm_f))).t() * exp(K_nm * inv_K_mm_f) +
        .5 * (
          inv_K_mm_f.t() * all_K_mm.slice(p) * inv_K_mm_f -
          sum(diag_1) +
          trace(inv_K_mm * inv_A * inv_K_mm * all_K_mm.slice(p)) -
          d_f_m(par, f, y, K_mm, W, M, lambda, m_size, all_K_mm.slice(p)).t() * inv_A_d3_v
        )
      );
    }
    out(1) = grad;
  }
  return(out);
}
