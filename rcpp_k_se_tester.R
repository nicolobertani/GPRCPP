library(RcppArmadillo)
library(microbenchmark)
library(purrr)
Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_diff_sq_mx.cpp")
Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_k_se_diff.cpp")
setwd("~/Dropbox (INSEAD)/Crime Modeling/")
source("GP functions - 20200417.R")

# rcpp_k_se_diff(cbind(c(1.2, 1.2, 4), c(1.2, 1.2, 0), c(2, 2, 2)), 1, 1, 1)
# rcpp_diff_sq_mx(as.matrix(c(1, 2, 3, 10, 0, 1.44)), as.matrix(c(2, 0)), 1)
# rcpp_diff_sq_mx(cbind(c(1, 2, 3, 10, 0, 1.44), c(1, 2, 3, 10, 0, 1.44)), cbind(c(2, 0), c(2, 0)), 1)
# rcpp_diff_sq_mx(cbind(c(1, 2, 3, 10, 0, 1.44), c(1, 2, 3, 10, 0, 1.44)), cbind(c(2, 0), c(2, 0)), F)
# rcpp_k_se_diff(rcpp_diff_sq_mx(cbind(c(1, 2, 3, 10, 0, 1.44), c(1, 2, 3, 10, 0, 1.44)), cbind(c(2, 0), c(2, 0)), F), 1, 1, F)


X <- matrix(rnorm(2000), ncol = 2)

all.equal(
  k.se(X),
  rcpp_k_se_diff(rcpp_diff_sq_mx(X, as.matrix(0), T), 1, 1, T)
)

system.time(
  for (i in 1:10) {
    rcpp_k_se_diff(rcpp_diff_sq_mx(X, as.matrix(0), T), 1, 1, T)
  } 
)
system.time(
  for (i in 1:10) {
    k.se(X)
  } 
)


X1 <- matrix(rnorm(1400), ncol = 2)
X2 <- matrix(rnorm(600), ncol = 2)

all.equal(
  k.se(X1, X2, .8, 1.2),
  rcpp_k_se_diff(rcpp_diff_sq_mx(X1, X2, F), .8, 1.2, F)
)

system.time(
  for (i in 1:10) {
  k.se(X1, X2, .8, 1.2)
  } 
)
system.time(
  for (i in 1:10) {
  rcpp_k_se_diff(rcpp_diff_sq_mx(X1, X2, F), .8, 1.2, F)
  } 
)

Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_k_se.cpp")

X <- matrix(rnorm(20), ncol = 2)
all.equal(
  k.se(X),
  rcpp_k_se(X, as.matrix(1), 1, 1, T)
)


all.equal(
  rcpp_k_se(as.matrix(c(1,2, 4)), as.matrix(0), 1, 1, T),
  rcpp_k_se_diff(rcpp_diff_sq_mx(as.matrix(c(1,2, 4)), as.matrix(0), T), 1, 1, T)
)
all.equal(
  rcpp_k_se(X, as.matrix(0), 1, 1, T),
  rcpp_k_se_diff(rcpp_diff_sq_mx(X, as.matrix(0), T), 1, 1, T)
)
all.equal(
  rcpp_k_se(X1, X2, 1, 1, F),
  rcpp_k_se_diff(rcpp_diff_sq_mx(X1, X2, F), 1, 1, F)
)
all.equal(
  rcpp_k_se(X, X, 1, 1, T),
  rcpp_k_se(X, X, 1, 1, F)
)
