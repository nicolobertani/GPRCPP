library(RcppArmadillo)
Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_k_se_diff.cpp")
rcpp_k_se_diff(cbind(c(1.2, 1.2, 4), c(1.2, 1.2, 0), c(2, 2, 2)), 1, 1, 1)

# k.se(sqrt(2), 0)