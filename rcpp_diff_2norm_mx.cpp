#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double entry_diff (const vec x, const vec y) {
  double res = norm(x - y, 2);
  return res;
}

// [[Rcpp::export]]
mat rcpp_diff_2norm_mx(const mat M, const mat N, const bool equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat D;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    D.set_size(M.n_rows, M.n_rows);
    D.fill(0);
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        D(r, c) = entry_diff(M.row(r).t(), M.row(c).t());
      }
    }
    D = D + D.t();

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    D.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        D(r, c) = entry_diff(M.row(r).t(), N.row(c).t());
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return D;
}
