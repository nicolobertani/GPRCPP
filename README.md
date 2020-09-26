# GPRCPP: Gaussian Processes with R and C++ via Rcpp and RcppArmadillo

In this directory you can find C++ implementations from scratch of the most common covariance functions for modeling with Gaussian Process priors and Gaussian or non-Gaussian likelihoods.
These functions include the Squared-Exponential, the Periodic, and the Linear kernels.

The implementation uses Rcpp and RcppArmadillo and can be used in a general R script.
This is particularly handy since it allows to use the optimization routines that come with R, as well as the more convenient visualization, etc.
The speed improvement offered by these functions in C++ is of roughly two orders of magnitude compared to the R equivalents.
More complex combinations of these simple kernels can be structured in R directly (e.g. the Automatic Relevance Kernel via multiplication of SEs).

For an introduction to these kernels, see:

*David Duvenaud, James Lloyd, Roger Grosse, Joshua Tenenbaum, Ghahramani Zoubin*,
**Structure Discovery in Nonparametric Regression through Compositional Kernel Search**,
Proceedings of the 30th International Conference on Machine Learning, PMLR 28(3):1166-1174, 2013.

For an exhaustive coverage of the topic of Gaussian Processes, the reference is:

*Rasmussen and Williams*,
**Gaussian Processes for Machine Learning**,
Volume 2, MIT press Cambridge, MA, 2006.

This is freely available [here](http://www.gaussianprocess.org/gpml/).

This project is ongoing. Contributions and suggestions are welcome.

### License

GPRCPP: Gaussian Processes with R and C++ via Rcpp and RcppArmadillo
Copyright (C) 2020  Nicolò Bertani

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>. You should also find a copy of this license in this repository.

Please use Github Issues to get in touch.
