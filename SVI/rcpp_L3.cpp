#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double rcpp_L3(const double sigma_sq,
  const vec &y, const vec &m, const mat &S, const mat &L,
  const mat &K_mn, const mat &inv_K_mm, const mat &K_tilde, const int batch_size) {
  double output;

  Function dnorm("dnorm")

  sum(sapply(seq(batch_size), function(i) {
      dnorm(y.svi[i], t(K.mn[, i]) %*% (inv.K.mm %*% m.svi), sqrt(sigma.sq.svi), log = T)
    })) -
      .5 / sigma.sq.svi * sum(diag(K.tilde)) -
      .5 * sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% (Lambda - inv.K.mm)[, j]})) -
      .5 * (sum(sapply(seq(m.size), function(j) {inv.K.mm[j, ] %*% S.svi[, j]})) +
              t(m.svi) %*% (inv.K.mm %*% m.svi) - m.size) -
      sum(log(diag(U.mm))) + sum(log(diag(U.S.svi)))

  return output;
}
