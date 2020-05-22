#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double k_lin(const vec &x, const vec &y, const double l, const double s) {
  double res = std::pow(s, 2) * dot(x - l, y - l);
  return res;
}

double k_se(const vec &x, const vec &y, const double m) {
  double df = pow(norm(x - y, 2), 2);
  double res = exp(- df / (2 * pow(m, 2)));
  return res;
}

double d_k_lin_l(const vec &x, const vec &y, const double l, const double s) {
  double res = - std::pow(s, 2) * (sum(x - l) + sum(y - l));
  return res;
}

double d_k_lin_s(const vec &x, const vec &y, const double l, const double s) {
  double res = 2 * s * dot(x - l, y - l);
  return res;
}

double d_k_se_m(const vec &x, const vec &y, const double m) {
  double df = pow(norm(x - y, 2), 2);
  double res = df / pow(m, 3) * exp(- df / (2 * pow(m, 2)));
  return res;
}

mat rcpp_k_lin(const mat &M, const mat &N, const double l, const double s, const bool equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = k_lin(M.row(r).t(), M.row(c).t(), l, s);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = k_lin(M.row(i).t(), M.row(i).t(), l, s);
    }

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = k_lin(M.row(r).t(), N.row(c).t(), l, s);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}

mat rcpp_k_se(const mat &M, const mat &N, const double m, const bool equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = k_se(M.row(r).t(), M.row(c).t(), m);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = k_se(M.row(i).t(), M.row(i).t(), m);
    }

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = k_se(M.row(r).t(), N.row(c).t(), m);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}

mat rcpp_d_k_lin_l(const mat &M, const mat &N, const double l, const double s, const bool equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = d_k_lin_l(M.row(r).t(), M.row(c).t(), l, s);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = d_k_lin_l(M.row(i).t(), M.row(i).t(), l, s);
    }

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_lin_l(M.row(r).t(), N.row(c).t(), l, s);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}

mat rcpp_d_k_lin_s(const mat &M, const mat &N, const double l, const double s, const bool equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = r + 1; c < M.n_rows; c++) {
        K(r, c) = d_k_lin_s(M.row(r).t(), M.row(c).t(), l, s);
      }
    }
    K = K + K.t();
    // fill diag
    for (int i = 0; i < M.n_rows; i++) {
      K(i,i) = d_k_lin_s(M.row(i).t(), M.row(i).t(), l, s);
    }

  } else {

    // Rcout << "Matrices are NOT equal.\n"; // progress message
    K.set_size(M.n_rows, N.n_rows);
    // fill everything
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_lin_s(M.row(r).t(), N.row(c).t(), l, s);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}

mat rcpp_d_k_se_m(const mat &M, const mat &N, const double m, const bool equal_matrices) {
  // Rcout << "I am running.\n"; // progress message
  mat K;

  if (equal_matrices == 1) {

    // Rcout << "Matrices are equal.\n"; // progress message
    K.set_size(M.n_rows, M.n_rows);
    K.zeros();
    // fill upper triangular wo diag
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
    for (int r = 0; r < M.n_rows; r++) {
      for (int c = 0; c < N.n_rows; c++) {
        K(r, c) = d_k_se_m(M.row(r).t(), N.row(c).t(), m);
      }
    }
  }
  // Rcout << "Ok, I got to the end.\n"; // progress message
  return K;
}

mat rcpp_d_k_linse_l(const mat &M, const mat &N, const mat &P, const mat &Q, const double l, const double s, const double m, const bool equal_matrices) {
  mat d_K_lin_l, K_se, res;
  d_K_lin_l = rcpp_d_k_lin_l(M, N, l, s, equal_matrices);
  K_se = rcpp_k_se(P, Q, m, equal_matrices);
  res = d_K_lin_l % K_se;
  return res;
}

mat rcpp_d_k_linse_s(const mat &M, const mat &N, const mat &P, const mat &Q, const double l, const double s, const double m, const bool equal_matrices) {
  mat d_K_lin_s, K_se, res;
  d_K_lin_s = rcpp_d_k_lin_s(M, N, l, s, equal_matrices);
  K_se = rcpp_k_se(P, Q, m, equal_matrices);
  res = d_K_lin_s % K_se;
  return res;
}

mat rcpp_d_k_linse_m(const mat &M, const mat &N, const mat &P, const mat &Q, const double l, const double s, const double m, const bool equal_matrices) {
  mat K_lin, d_K_se_m, res;
  K_lin = rcpp_k_lin(M, N, l, s, equal_matrices);
  d_K_se_m = rcpp_d_k_se_m(P, Q, m, equal_matrices);
  res = K_lin % d_K_se_m;
  return res;
}

mat d_sum_Lambda_fun(const mat &M1, const mat &M1D, const mat &M2, const mat &M2D, const double sigma_sq) {
  mat D;
  D = M1D * M2  * M2.t()  * M1 +
      M1  * M2D * M2.t()  * M1 + // problem here
      M1  * M2  * M2D.t() * M1 +
      M1  * M2  * M2.t()  * M1D;
  D = D / sigma_sq;
  return D;
}

// [[Rcpp::export]]
double rcpp_step_k_linse_svi(const mat &M, const mat &N, const mat &P, const mat &Q,
  const double l, const double s, const double m, double sigma_sq,
  const mat &inv_K_mm, const mat &K_mn, const vec &y, const vec &mu, const mat &S, const int par_id) {
    // Rcout << "I am running.\n"; // progress message
    mat d_K_mm(M.n_rows, M.n_rows, fill::zeros);
    mat d_inv_K_mm(M.n_rows, M.n_rows, fill::zeros);
    mat d_K_mn(M.n_rows, N.n_rows, fill::zeros);
    vec diag_d_K_nn(N.n_rows);
    double tr_d_K_nn;
    mat d_Q_nn(N.n_rows, N.n_rows, fill::zeros);
    mat d_sum_Lambda(M.n_rows, M.n_rows, fill::zeros);
    vec S_d_sum_Lambda(mu.n_elem);
    double output;

    if (par_id == 1) {

      d_K_mm = rcpp_d_k_linse_l(M, N, P, Q, l, s, m, 1);
      d_K_mn = rcpp_d_k_linse_l(M, N, P, Q, l, s, m, 0);
      for (size_t i = 0; i < N.n_rows; i++) {
        diag_d_K_nn(i) = d_k_lin_l(N.row(i).t(), N.row(i).t(), l, s) * k_se(Q.row(i).t(), Q.row(i).t(), m);
      }
    } else if (par_id == 2) {

      d_K_mm = rcpp_d_k_linse_s(M, N, P, Q, l, s, m, 1);
      d_K_mn = rcpp_d_k_linse_s(M, N, P, Q, l, s, m, 0);
      for (size_t i = 0; i < N.n_rows; i++) {
        diag_d_K_nn(i) = d_k_lin_s(N.row(i).t(), N.row(i).t(), l, s) * k_se(Q.row(i).t(), Q.row(i).t(), m);
      }
    } else {

      d_K_mm = rcpp_d_k_linse_m(M, N, P, Q, l, s, m, 1);
      d_K_mn = rcpp_d_k_linse_m(M, N, P, Q, l, s, m, 0);
      for (size_t i = 0; i < N.n_rows; i++) {
        diag_d_K_nn(i) = k_lin(N.row(i).t(), N.row(i).t(), l, s) * d_k_se_m(Q.row(i).t(), Q.row(i).t(), m);
      }
    }

    d_inv_K_mm = - inv_K_mm * d_K_mm * inv_K_mm;
    tr_d_K_nn = sum(diag_d_K_nn);
    d_Q_nn = (d_K_mn.t() * inv_K_mm * K_mn) + (K_mn.t() * d_inv_K_mm * K_mn) + (K_mn.t() * inv_K_mm * d_K_mn);
    d_sum_Lambda = d_sum_Lambda_fun(inv_K_mm, d_inv_K_mm, K_mn, d_K_mn, sigma_sq);
    for (size_t i = 0; i < mu.n_elem; i++) {
      S_d_sum_Lambda(i) = dot(S.row(i), d_sum_Lambda.col(i));
    }

    output = as_scalar(
      (y - K_mn.t() * (inv_K_mm * mu)).t() * (d_K_mn.t() * (inv_K_mm * mu) + K_mn.t() * (d_inv_K_mm * mu)) / sigma_sq -
      (tr_d_K_nn - trace(d_Q_nn)) / (2 * sigma_sq) -
      sum(S_d_sum_Lambda) / 2 +
      (trace((inv_K_mm * S * inv_K_mm - inv_K_mm) * d_K_mm) - mu.t() * (d_inv_K_mm * mu)) / 2
    );
    // Rcout << "Ok, I got to the end.\n"; // progress message
    return output;
}
