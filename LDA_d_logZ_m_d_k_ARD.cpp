#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
vec rcpp_d3_likelihood_m (const vec &y, const vec &f, const mat &M, const int i) {
  vec E = exp(M * f);
  int m_size = f.n_elem;
  vec v_out = zeros<vec>(m_size);
  for (size_t k = 0; k < m_size; k++) {
    v_out(k) = - dot(E, M.col(k) % M.col(k) % M.col(i));
  }
  return(v_out);
}
