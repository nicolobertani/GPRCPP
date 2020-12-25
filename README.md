Overview
--------

In this directory you can find C++ implementations from scratch of the
most common covariance functions for modeling with Gaussian Process
priors and Gaussian or non-Gaussian likelihoods. These functions include
the Squared-Exponential, the Periodic, and the Linear kernels.

The implementation uses Rcpp and RcppArmadillo and can be used in a
general R script. This is particularly handy since it allows to use the
optimization routines that come with R, as well as the more convenient
visualization, etc. The speed improvement offered by these functions in
C++ is of roughly two orders of magnitude compared to the R equivalents.
More complex combinations of these simple kernels can be structured in R
directly (e.g. the Automatic Relevance Kernel via multiplication of
SEs).

For an introduction to these kernels, see:

*David Duvenaud, James Lloyd, Roger Grosse, Joshua Tenenbaum, Ghahramani
Zoubin*, **Structure Discovery in Nonparametric Regression through
Compositional Kernel Search**, Proceedings of the 30th International
Conference on Machine Learning, PMLR 28(3):1166-1174, 2013.

For an exhaustive coverage of the topic of Gaussian Processes, the
reference is:

*Rasmussen and Williams*, **Gaussian Processes for Machine Learning**,
Volume 2, MIT press Cambridge, MA, 2006.

This is freely available [here](http://www.gaussianprocess.org/gpml/).

This project is ongoing. Contributions and suggestions are welcome.

License
-------

GPRCPP: Gaussian Processes with R and C++ via Rcpp and RcppArmadillo

Copyright (C) 2020 Nicolò Bertani

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or any later
version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see
<a href="https://www.gnu.org/licenses/" class="uri">https://www.gnu.org/licenses/</a>.
You should also find a copy of this license in this repository.

Please open an [*issue*](https://github.com/nicolobertani/GPRCPP/issues)
if you have any trouble with the code.

Details
-------

To my knowledge, existing R packages for Gaussian Processes only cover a
limited number of models (GP regression and a few more things). To work
with less vanilla models (e.g. count data, large-scale approximations),
you probably will have to code things from scratch. In my opinion, using
packages is not the best approach to fully understand a model, anyhow.
So this is probably a good thing.

Regardless of the GP model you are using, you will need functions to
compute the covariance matrix and its derivatives. These functions must
*for*-loop through every entry of the matrix and calculate the
appropriate value. Since these matrices tend to be large, this operation
generates non-negligible computational burden. Additionally, this burden
grows with the number of covariates and of kernels in the covariance
function.

In R, *for* loops are very slow. Even using looping functions and
exploiting matrix symmetry, the creation of these matrices ends up
taking the lion share of the computing time, in my experience. Because
of this, these functions are natural candidates to be implemented in
C++. They can then be easily embedded in an R estimation procedure,
replacing the equivalent R functions, with significant speed gains.

I am going to do three things here:

1.  Micro-bench the R and C++ implementation of the Squared-Exponential
    covariance function.

2.  Showcase how they can be embedded in and speed up two R routine to
    perform GP regression and estimation of a GP model with count data.
    The implementations follow chapters 2 and 3 of the aforementioned
    Rasmussen and Williams’ book.

<!-- ### Microbenching R and C++ implementations -->
<!-- ### Example 1: GP regression - GP prior and Normal likelihood -->
<!-- ### Example 2: GP prior with Poisson likelihood -->
