library(RcppArmadillo)
Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_k_lin.cpp")
setwd("~/Dropbox (INSEAD)/Crime Modeling/")
source("GP functions - 20200417.R")

X <- matrix(rnorm(10), ncol = 2)
rcpp_k_lin(X, X, .5, 1, T)
k.lin(X, X, .5, 1)

X <- matrix(rnorm(2000), ncol = 2)
all.equal(
  k.lin(X, l = .5, s = .2),
  rcpp_k_lin(X, X, .5, .2, T)
)

system.time(
  for (i in 1:1000) {
    rcpp_k_lin(X, X, .5, .2, T)
  }
)
system.time(
  for (i in 1:10) {
    k.lin(X, l = .5, s = .2)
  }
)

X1 <- matrix(rnorm(1400), ncol = 2)
X2 <- matrix(rnorm(600), ncol = 2)

all.equal(
  k.lin(X1, X2, .5, .2),
  rcpp_k_lin(X1, X2, .5, .2, F)
)
all.equal(
  k.lin(X1, X2, .8, 1.2),
  rcpp_k_lin(X1, X2, .8, 1.2, F)
)

all.equal(
  rcpp_k_lin(X, X, 1, 1, T),
  rcpp_k_lin(X, X, 1, 1, F)
)
