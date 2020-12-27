# GPRCPP: Gaussian Processes with R and C++ via Rcpp and RcppArmadillo

## Overview

In this directory you can find C++ implementations from scratch of the
most common covariance functions and their derivatives for modeling with
Gaussian Process priors. These functions include the
squared-exponential, the periodic, and the linear kernels. Derivatives
are provided with respect to the characteristic length-scale and
amplitude.

The implementation uses Rcpp and RcppArmadillo and can be used in a
general R script. This is particularly handy since it allows to use the
efficient functions in any GP model and to reply on R for optimization
routines, convenient data handling and visualization, etc. The speed
improvement offered by these functions in C++ is of roughly two orders
of magnitude compared to the R equivalents. More complex combinations of
these simple kernels can be structured in R directly (e.g. building the
Automatic Relevance Determination kernel by multiplicatying SEs).

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

## License

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
with this program. If not, see <https://www.gnu.org/licenses/>. You
should also find a copy of this license in this repository.

Please open an [*issue*](https://github.com/nicolobertani/GPRCPP/issues)
if you have any trouble with the code.

## Details

To my knowledge, existing R packages for Gaussian Processes only cover a
limited number of models (GP regression and a few more things). To work
with less vanilla models (e.g. count data, large-scale approximations),
you probably will have to code things from scratch. In my opinion, this
is probably a good thing and the recommended approach if you actually
want to understand the model.

Regardless of the GP model you are using, you will need functions to
compute the covariance matrix and its derivatives. These functions must
*for*-loop through every entry of these matrices and calculate the
appropriate value. Since these matrices tend to be large, this operation
generates non-negligible computational burden. Additionally, the burden
grows with the number of covariates and of kernels in the covariance
function.

In R, *for* loops are known to be slow. Using looping functions and
exploiting matrix symmetry does not seem to solve the problem. In my
experience, the creation of the covariance matrix and its derivatives
ends up taking the lion share of the computing time (yes, more that
matrix inversion). Because of this, relying of efficient C++ functions
to create these matrices seems a natural way to proceed. Using Rcpp and
RcppArmadillo, these functions can then be easily embedded in any GP
estimation procedure in R, with significant speed gains. As mentioned,
here you will find C++ implementations the following covariance kernels
and their derivatives: squared-exponential, periodic, and linear. I am
going to:

1.  Explain the naming conventions of the scripts in this repository and
    of the functions they import in R.

2.  Micro-bench the R and C++ implementations of the Squared-Exponential
    covariance function.

3.  Showcase how they can be embedded in and speed up an R routine to
    perform GP regression, following Chapter 2 of the aforementioned
    Rasmussen and Williams’ book.

### Naming convention and function arguments

The names of the scripts follow this structure:

1.  everything is prefixed by *rcpp\_*.

2.  *k\_* for kernels, *d\_k\_* for kernel derivatives.

3.  *se* for squared-exponential, *per* for periodic, *lin* for linear
    kernel.

4.  kernel derivatives are suffixed by *\_l* for the derivative with
    respect to the characteristic length-scale, and *\_s* for amplitude.

For instance, the script containing the derivative with respect to the
length-scale of a squared-exponential kernel is *rcpp\_d\_se\_l*.

The scripts can be loaded using the function `Rcpp::sourceCpp()` and
will import homonymous functions. Every function will require two
matrices (vectors need to me input using the R function `as.matrix()`),
*l* and *s* parameters, and a boolen value to indicate whether the 2
matrices are the same and symmetry can be exploited. For the periodic
kernel, after *l* and *s* you need to input a value *p* for the period.
All functions return a matrix.

### Microbenching R and C++ implementations

To give you a sense of how much faster the C++ implementation can be, I
am going to compare it to two R alternatives. All functions will
calculate the covariance function between a matrix and itself, where the
covariance function is a simple squared-exponential kernel. I am using a
2-column matrix of 100 rows.

``` r
X <- matrix(rnorm(200), ncol = 2)
```

The first R implementation simply *for*-loops through the entire matrix,
without exploiting symmetry (called “*for-loop*”).

``` r
# R for-loop IMPLEMENTATION
k.se.forloop <- function(X1, X2 = NULL, l = 1, s = 1) {
  k <- function(x, y) {
    s ^ 2 * exp(- sum((x - y) ^ 2) / (2 * l ^ 2)) 
  }
  x1 <- as.matrix(X1)
  x2 <- as.matrix(X2)
  output <- matrix(NA, nrow = nrow(x1), ncol = nrow(x2))
  for (i in seq(nrow(x1))) {
    for (j in seq(nrow(x2))) {
      output[i, j] <- k(x1[i, ], x2[j, ])
    }
  }
  return(output)
}
```

The second R implementation uses looping function `sapply` and exploits
symmetry when possible (called “*looping*”).

``` r
# R looping IMPLEMENTATION
k.se.looping <- function(X1, X2 = NULL, l = 1, s = 1) {
  k <- function(x, y) {
    s ^ 2 * exp(- sum((x - y) ^ 2) / (2 * l ^ 2)) 
  }
  x1 <- as.matrix(X1)
  if (!length(X2)) { # exploit symmetry
    output <- matrix(0, nrow = nrow(x1), ncol = nrow(x1))
    sapply(seq(nrow(x1)), function(i) {
      sapply(seq(nrow(x1))[seq(nrow(x1)) > i], function(j) {
        output[i, j] <<- k(x1[i, ], x1[j, ]) 
      })})
    output <- output + t(output)
    diag(output) <- s ^ 2
  } else {
    x2 <- as.matrix(X2)
    output <- matrix(0, nrow = nrow(x1), ncol = nrow(x2))
    sapply(seq(nrow(x1)), function(i) {
      sapply(seq(nrow(x2)), function(j) {
        output[i, j] <<- k(x1[i, ], x2[j, ]) 
      })})
  }
  return(output)
}
```

The C++ alternative is simply loaded as follows.

``` r
# LOAD CPP IMPLEMENTATION CALLED rcpp_k_se
Rcpp::sourceCpp("rcpp_k_se.cpp")
```

Without exploiting symmetry, the runtimes are the following.

``` r
microbenchmark::microbenchmark(  
  "R for-loop" = (k.se.forloop(X, X, l = .8, s = 1.2)),
  "R looping" = (k.se.looping(X, X, l = .8, s = 1.2)),
  "C++" = (rcpp_k_se(X, X, l = .8, s = 1.2, equal_matrices = F)),
  check = "equal"
)
```

    Unit: microseconds
           expr   min    lq  mean median    uq   max neval
     R for-loop 17317 19739 21597  20691 22102 44029   100
      R looping 24738 28261 29918  29255 31102 51326   100
            C++   518   586   638    606   648  1880   100

The best R alternative is about 36 times slower (median runtime). Using
looping functions does not seem to help.

If, in addition, we exploit symmetry, we get:

``` r
microbenchmark::microbenchmark(
  "R for-loop" = (k.se.forloop(X, X, l = .8, s = 1.2)),
  "R looping" = (k.se.looping(X, l = .8, s = 1.2)),
  "C++" = (rcpp_k_se(X, X, l = .8, s = 1.2, equal_matrices = T)),
  check = "equal"
)
```

    Unit: microseconds
           expr   min    lq  mean median    uq   max neval
     R for-loop 16809 19254 21641  20764 22444 41212   100
      R looping 12652 15055 16459  15893 16810 30052   100
            C++   249   305   342    322   353   581   100

where, in the C++ implementation, we reduce runtime by a bit less than
half.

### Example: GP regression - GP prior and Normal likelihood

    Loading required package: MASS

``` r
# SIMULATE DATA
population.size <- 5e2
x.limits <- c(-5, 5)
X <- as.matrix(seq(x.limits[1], x.limits[2], length.out = population.size))
K.true <- rcpp_k_se(X, X, 2, .4, T)
y <- mvrnorm(1, rep(0, population.size), K.true)

# SIMULATE SAMPLE
sample.size <- 1e2
positions <- sort(sample(seq(population.size), sample.size))
X.sample <- as.matrix(X[positions])
y.sample <- y[positions] + rnorm(sample.size)
```
