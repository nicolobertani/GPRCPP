#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double sq2pi = pow(2 * datum::pi, .5);

// [[Rcpp::export]]
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
