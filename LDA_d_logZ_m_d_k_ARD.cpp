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
/*
d.f.m <- function(par.hat, f.m = f.m.hat) {
  K.mm <- K.mm.f(par.hat)
  chol.K.mm <- chol(K.mm)
  inv.K.mm <- chol2inv(chol.K.mm)
  K.nm <- K.nm.f(par.hat)
  M <- K.nm %*% inv.K.mm
  W <- - d2.likelihood.m(y.sample, f.m, M)
  d.K.list <- self.derivative.list(par.hat)
  dfm <- lapply(d.K.list, function(d.mx) {
    solve(diag(m.size) + K.mm %*% W, d.mx %*% d.likelihood.m(y.sample, f.m, M))
  })
  return(dfm)
}
*/
