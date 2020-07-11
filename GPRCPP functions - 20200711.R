if (!requireNamespace('RcppArmadillo')) install.packages('RcppArmadillo')
library(RcppArmadillo)
setwd("~/Dropbox (INSEAD)/Crime Modeling/CPP code/")
# setwd("~/crime.modeling/Rscripts/helper_functions/")


# SE - Function and derivatives ----------------------------------------------------

# Rcpp::sourceCpp("rcpp_diff_sq_mx.cpp")
# Rcpp::sourceCpp("rcpp_k_se_diff.cpp")
Rcpp::sourceCpp("rcpp_k_se.cpp")
Rcpp::sourceCpp("rcpp_d_k_se_l.cpp")
Rcpp::sourceCpp("rcpp_d_k_se_s.cpp")


# Periodic - Function and derivatives -------------------------------------

# Rcpp::sourceCpp("rcpp_diff_2norm_mx.cpp")
# Rcpp::sourceCpp("rcpp_k_per_diff.cpp")
Rcpp::sourceCpp("rcpp_k_per.cpp")
Rcpp::sourceCpp("rcpp_d_k_per_l.cpp")
Rcpp::sourceCpp("rcpp_d_k_per_s.cpp")


# # Linear - Function and derivatives -----------------------------------
# 
# Rcpp::sourceCpp("rcpp_k_lin.cpp")
# # Rcpp::sourceCpp("rcpp_d_k_lin_l.cpp")
# # Rcpp::sourceCpp("rcpp_d_k_lin_s.cpp")
# 
# # SExSE - Function and derivatives ----------------------------------------------------
# 
# Rcpp::sourceCpp("rcpp_k_sese.cpp")
# 
# 
# # PerxSE - Function and derivatives ----------------------------------------------------
# 
# Rcpp::sourceCpp("rcpp_k_perse.cpp")
# 
# 
# # LinxSE - Function and derivatives ----------------------------------------------------
# 
# Rcpp::sourceCpp("rcpp_k_linse.cpp")
# 
# 

# ARD ---------------------------------------------------------------------

Rcpp::sourceCpp("rcpp_d_k_ARD.cpp")

