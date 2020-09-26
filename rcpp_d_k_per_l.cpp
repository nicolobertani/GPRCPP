/*
GPRCPP: Gaussian Processes with R and C++ via Rcpp and RcppArmadillo
Copyright (C) 2020  Nicolò Bertani

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double d_k_per_l(const vec &x, const vec &y, const double l, const double s, const int p) {
  double df = norm(x - y, 2);
  double res = std::pow(s, 2) * exp(- 2 * std::pow(sin(datum::pi * df / p), 2) / std::pow(l, 2)) * 4 * std::pow(sin(datum::pi * df / p), 2) / std::pow(l, 3);
  return res;
}


// [[Rcpp::export]]
mat rcpp_d_k_per_l(const mat &M, const mat &N, const double l, const double s, const int p, bool equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = d_k_per_l(M.row(r).t(), M.row(c).t(), l, s, p);
      }
    }
    K = K + K.t();
    // fill diag
    /*for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = d_k_per_l(M.row(i).t(), M.row(i).t(), l, s, p);
    }*/

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_per_l(M.row(r).t(), N.row(c).t(), l, s, p);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}
