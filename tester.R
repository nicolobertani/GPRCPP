library(RcppArmadillo)
# Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_k_se_diff.cpp")
# rcpp_k_se_diff(cbind(c(1.2, 1.2, 4), c(1.2, 1.2, 0), c(2, 2, 2)), 1, 1, 1)

Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_diff_mx.cpp")
rcpp_diff_mx(as.matrix(c(1, 2, 3, 10, 0, 1.44)), as.matrix(c(2, 0)), 1)
rcpp_diff_mx(cbind(c(1, 2, 3, 10, 0, 1.44), c(1, 2, 3, 10, 0, 1.44)),
cbind(c(2, 0), c(2, 0)), 1)
