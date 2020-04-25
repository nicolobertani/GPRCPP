#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*double k_se(const vec x, const vec y, const double l = 1, const double s = 1) {
  double res = std::pow(s, 2) * exp(pow(norm(x - y, 2), 2) / (2 * std::pow(l, 2)));
  return res;
}*/

double k_se_diff(double df, const double l, const double s) {
  double res = std::pow(s, 2) * exp(- df / (2 * std::pow(l, 2)));
  return res;
}

// [[Rcpp::export]]
mat rcpp_k_se(const mat &A, const mat &B,
  const double &l = 1, const double &s = 1,
  bool equal_matrices = 1) {

  mat K;
  if (equal_matrices) {

    Rcout << "I am executing TRUE." << "\n";

    K.set_size(A.n_rows, A.n_rows);
    K.fill(k_se_diff(1, 1, 1));

  } else {

    Rcout << "I am executing FALSE." << "\n";

    K.set_size(A.n_rows, B.n_rows);
    K.fill(k_se_diff(1, 1, 1));

  }

  Rcout << "I ACTUALLY GOT TO THE END!\n";
  return K;
}
