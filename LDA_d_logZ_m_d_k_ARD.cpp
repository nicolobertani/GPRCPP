#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

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


// [[Rcpp::export]]
