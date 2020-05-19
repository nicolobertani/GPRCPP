#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double rcpp_step_sigma_sq_svi(const double sigma_sq,
  const vec &y, const vec &m, const mat &S, const mat &L,
  const mat &K_mn, const mat &inv_K_mm, const mat &K_tilde, const int batch_size) {
  double output;
  vec z;
  mat sl_term2;
  vec sl(m.n_elem);
  double trace_sl;

  sl_term2 = sigma_sq * (L - inv_K_mm);
  for (size_t i = 0; i < m.n_elem; i++) {
    sl(i) = dot(S.row(i), sl_term2.col(i));
  }
  trace_sl = sum(sl);
  z = y - K_mn.t() * (inv_K_mm * m);
  output = - batch_size / (2 * sigma_sq) + (dot(z, z) + trace(K_tilde) + trace_sl) / (2 * pow(sigma_sq, 2));

  return output;
}
