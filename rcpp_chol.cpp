#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
mat rcpp_chol(const mat M) {
  // Rcout << "I am running.\n"; // progress message
  mat cM = chol(M);
  return cM;
}
