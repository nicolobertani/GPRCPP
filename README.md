# GPRCPP: Gaussian Processes with R and C++

In this directory you can find C++ implementations of the most common covariance functions for modelling with Gaussian Process priors and Gaussian or non-Gaussian likelihoods.
These functions include the squared exponential, the periodic, the linear, and the ARD (automatic relevance determination) kernels.
The implementation uses Rcpp and RcppArmadillo and can be used in a general R script.
This is particularly handy since it allows to use the optimization routines that come with R, as well as the more convienient visualization, etc. 
The speed improvement offered by these functions in C++ is of roughly two orders of magnitude compared to the R equivalents.

An example of how to use this code for GP regression with hyperparameter estimation is exemplified in EXAMPLE.R.
The algorithms are taken from *Rasmussen and Williams, 2006, Gaussian Processes for Machine Learning* (freely available [here](http://www.gaussianprocess.org/gpml/)).

This code is shared freely and is provided without any guarantee whatsoever. 
Contributions and suggestions are welcome.
