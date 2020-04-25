library(RcppArmadillo)
Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_diff_2norm_mx.cpp")
Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_k_per_diff.cpp")
Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_k_per.cpp")
setwd("~/Dropbox (INSEAD)/Crime Modeling/")
source("GP functions - 20200417.R")


X <- as.matrix(c(0, 1))
all.equal(
  k.per(X, X, l = .4, s = 4, p = 2),
  rcpp_k_per_diff(rcpp_diff_2norm_mx(X, as.matrix(0), T), .4, 4, 2, T)
)


X <- matrix(rnorm(1000), ncol = 1)
all.equal(
  k.per(X, X, l = .4, s = 4, p = 2),
  rcpp_k_per_diff(rcpp_diff_2norm_mx(X, as.matrix(0), T), .4, 4, 2, T)
)

system.time(
  for (i in 1:10) {
    rcpp_k_per_diff(rcpp_diff_2norm_mx(X, as.matrix(0), T), .4, 4, 2, T)
  }
)
system.time(
  for (i in 1:10) {
    k.per(X, X, l = .4, s = 4, p = 2)
  }
)


X1 <- matrix(rnorm(400), ncol = 2)
X2 <- matrix(rnorm(600), ncol = 2)

all.equal(
  k.per(X1, X2, .8, 1.2, 2),
  rcpp_k_per_diff(rcpp_diff_2norm_mx(X1, X2, F), .8, 1.2, 2, F)
)

system.time(
  for (i in 1:10) {
    k.per(X1, X2, .8, 1.2, 2)
  }
)
system.time(
  for (i in 1:10) {
    rcpp_k_per_diff(rcpp_diff_2norm_mx(X1, X2, F), .8, 1.2, 2, F)
  }
)



all.equal(
  rcpp_k_per(as.matrix(c(1, 2, 4)), as.matrix(0), .4, 4, 2, T),
  rcpp_k_per_diff(rcpp_diff_2norm_mx(as.matrix(c(1, 2, 4)), as.matrix(0), T), .4, 4, 2, T)
)
all.equal(
  rcpp_k_per(X, as.matrix(0), .4, 4, 2, T),
  rcpp_k_per_diff(rcpp_diff_2norm_mx(X, as.matrix(0), T), .4, 4, 2, T)
)
all.equal(
  rcpp_k_per(X1, X2, .4, 4, 2, F),
  rcpp_k_per_diff(rcpp_diff_2norm_mx(X1, X2, F), .4, 4, 2, F)
)
all.equal(
  rcpp_k_per(X, X, 1, 1, 1, T),
  rcpp_k_per(X, X, 1, 1, 1, F)
)
