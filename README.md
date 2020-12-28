# GPRCPP: Gaussian Processes with R and C++ via Rcpp and RcppArmadillo

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
Automatic Relevance Determination kernel by multiplying SEs).

Using the efficient C++ functions made available in this repository
instead of the R equivalents, runtime for model estimation in GP
regression is reduced by more than 90%. All details are discussed below.

``` r
microbenchmark::microbenchmark(
  fit.model(use.cpp = T),
  fit.model(use.cpp = F),
  times = 10
)
```

    Unit: seconds
                       expr   min      lq   mean median     uq    max neval
     fit.model(use.cpp = T)  7.24   9.581  17.02  15.42  20.21  32.46    10
     fit.model(use.cpp = F) 45.31 169.495 224.90 208.25 316.16 432.82    10

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
    perform GP regression, following Chapters 2 and 5 of the
    aforementioned Rasmussen and Williams’ book.

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
*l* and *s* parameters, and a Boolean value to indicate whether the 2
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
           expr     min      lq    mean  median      uq     max neval
     R for-loop 17591.9 19077.7 19864.7 20158.7 20501.6 23877.4   100
      R looping 26631.4 27513.4 28554.7 28469.4 28957.8 56102.9   100
            C++   555.4   578.3   590.2   583.4   598.5   661.2   100

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
           expr     min      lq    mean  median      uq     max neval
     R for-loop 18116.5 19018.9 20111.2 19904.4 20396.2 49684.6   100
      R looping 14314.0 14730.2 15337.7 15035.1 16058.3 17000.3   100
            C++   294.5   303.4   317.4   308.1   323.4   386.1   100

where, in the C++ implementation, we reduce runtime by a bit less than
half.

### Example: GP regression - GP prior and Normal likelihood

To perform GP regression with hyperparameter estimation, the model
evidence or marginal likelihood is iteratively maximized. Computing the
posterior of the GP functionals, the model evidence and its derivatives
requires to update the covariance function and its derivatives at every
iteration. The code below implements GP regression. Details are omitted
and can be found, as said, in Rasmussen and Williams (2006), chapters 2
and 5. The estimation algorithm is written in R. We compare runtime of
the same algorithm when using R or C++ functions to create the
covariance matrix and derivatives. The covariance function is set to a
simple squared-exponential, as above in the microbenching section. The
derivatives can be coded in R or loaded from the C++ implementations as
follows.

``` r
# R DERIVATIVES OF SE
d.k.se.l <- function(X1, X2 = NULL, l = 1, s = 1) {
  k <- function(x, y) {
    s ^ 2 * (sum((x - y) ^ 2) / (l ^ 3)) * exp(- sum((x - y) ^ 2) / (2 * l ^ 2)) # changed here
  }
  x1 <- as.matrix(X1)
  if (!length(X2)) {
    output <- matrix(0, nrow = nrow(x1), ncol = nrow(x1))
    sapply(seq(nrow(x1)), function(i) {
      sapply(seq(nrow(x1))[seq(nrow(x1)) > i], function(j) {
        output[i, j] <<- k(x1[i, ], x1[j, ]) 
      })})
    output <- output + t(output)
    diag(output) <- 0 # changed here
  } else {
    x2 <- as.matrix(X2)
    combs <- expand.grid(seq(nrow(x1)), seq(nrow(x2)))
    output <- mapply(function(i, j) {
      k(x1[i, ], x2[j, ]) 
    }, i = combs[, 1], j = combs[, 2])
    output <- matrix(output, nrow = nrow(x1))
  }
  return(output)
}

d.k.se.s <- function(X1, X2 = NULL, l = 1, s = 1) {
  k <- function(x, y) {
    2 * s * exp(- sum((x - y) ^ 2) / (2 * l ^ 2)) # changed here
  }
  x1 <- as.matrix(X1)
  if (!length(X2)) {
    output <- matrix(0, nrow = nrow(x1), ncol = nrow(x1))
    sapply(seq(nrow(x1)), function(i) {
      sapply(seq(nrow(x1))[seq(nrow(x1)) > i], function(j) {
        output[i, j] <<- k(x1[i, ], x1[j, ]) 
      })})
    output <- output + t(output)
    diag(output) <- 2 * s # changed here
  } else {
    x2 <- as.matrix(X2)
    combs <- expand.grid(seq(nrow(x1)), seq(nrow(x2)))
    output <- mapply(function(i, j) {
      k(x1[i, ], x2[j, ]) 
    }, i = combs[, 1], j = combs[, 2])
    output <- matrix(output, nrow = nrow(x1))
  }
  return(output)
}

# LOAD CPP DERIVATIVES OF SE
Rcpp::sourceCpp("rcpp_d_k_se_l.cpp")
Rcpp::sourceCpp("rcpp_d_k_se_s.cpp")
```

For this comparison, I consider a random sample of 500 observations
drawn from a GP prior with standard Normal noise. True characteristic
length-scale *l* is set to \(1.2\), and amplitude *s* to \(.8\).

To fit the observations and recover the hyperparameters, the following
algorithm can be used.

``` r
# FUNCTION TO SIMULATE SAMPLE AND FIT MODEL
fit.model <- function(use.cpp) {
  # SIMULATE SAMPLE
  sample.size <<- 5e2
  x.limits <- c(-5, 5)
  X <- as.matrix(seq(x.limits[1], x.limits[2], length.out = sample.size))
  if (use.cpp) {
    K.true <- rcpp_k_se(X, X, 1.2, .8, T)
  } else {
    K.true <- k.se.looping(X, l = 1.2, s = .8)
  }
  y <- mvrnorm(1, rep(0, sample.size), K.true)
  
  X.sample <<- as.matrix(X)
  y.sample <<- y + rnorm(sample.size)
  
  if (use.cpp) {
    k.hat.call <- parse(text = "rcpp_k_se(X.sample, X.sample, par.hat[1], par.hat[2], T)")
    step.calls <- lapply(list(
      'sum(sapply(seq(sample.size), function(i, D) {M[i, ] %*% D[, i]}, D = rcpp_d_k_se_l(X.sample, X.sample, par.hat[1], par.hat[2], T))) / 2',
      'sum(sapply(seq(sample.size), function(i, D) {M[i, ] %*% D[, i]}, D = rcpp_d_k_se_s(X.sample, X.sample, par.hat[1], par.hat[2], T))) / 2',
      'sum(diag(M)) / 2'
    ), function(i) parse(text = i))
  } else {
    k.hat.call <- parse(text = "k.se.looping(X.sample, l = par.hat[1], s = par.hat[2])")
    step.calls <- lapply(list(
      'sum(sapply(seq(sample.size), function(i, D) {M[i, ] %*% D[, i]}, D = d.k.se.l(X.sample, l = par.hat[1], s = par.hat[2]))) / 2',
      'sum(sapply(seq(sample.size), function(i, D) {M[i, ] %*% D[, i]}, D = d.k.se.s(X.sample, l = par.hat[1], s = par.hat[2]))) / 2',
      'sum(diag(M)) / 2'
    ), function(i) parse(text = i))
  }
  
  par.hat <<- rep(1, length(step.calls))
  threshold <- 1e-4
  iteration <- 1
  log.lik.sequence <- -Inf

  repeat({
    # update iteration
    iteration <- iteration + 1
    # print(c(iteration, log.lik.sequence[iteration - 1], par.hat))
    
    K.hat <<- eval(k.hat.call)
    U <<- chol(K.hat + diag(tail(par.hat, 1), sample.size)) 
    inv.K.hat <<- chol2inv(U)
    a <<- inv.K.hat %*% y.sample
    M <<- a %*% t(a) - inv.K.hat
    step <- unlist(lapply(step.calls, eval))
    par.hat <<- par.hat + 1e-3 * step
    
    # compute batch likelihood
    log.lik.sequence[iteration] <- 
      - sample.size / 2 * log(2 * pi) - sum(log(diag(U))) - t(y.sample) %*% (inv.K.hat %*% y.sample) / 2 
    if (abs(diff(tail(log.lik.sequence, 2))) < threshold) break
  })
  
  K.hat <- eval(k.hat.call)
  U <- chol(K.hat + diag(par.hat[3], sample.size)) 
  z <- forwardsolve(t(U), y.sample)
  a <- backsolve(U, z)
  mu.n <- c(K.hat %*% a) # POSTERIOR PREDICTIVE MEAN
  v <- forwardsolve(t(U), K.hat)
  Sigma <- zapsmall(K.hat - t(v) %*% v) # POSTERIOR PREDICTIVE VARIANCE
}
```

Using 10 different random samples, we can compare runtime when using C++
or R functions.

``` r
microbenchmark::microbenchmark(
  fit.model(use.cpp = T),
  fit.model(use.cpp = F),
  times = 10
)
```

    Unit: seconds
                       expr   min      lq   mean median     uq    max neval
     fit.model(use.cpp = T)  7.24   9.581  17.02  15.42  20.21  32.46    10
     fit.model(use.cpp = F) 45.31 169.495 224.90 208.25 316.16 432.82    10

Runtime is reduced by more than 90%. This is surprising, since the
theoretical bottleneck of GP modeling is matrix inversion and
determinant calculation. These results, however, clearly show the
importance of also considering the practical aspects of model
estimation. In fact, matrix inversion (or Choleski decomposition) in R
relies on efficient LAPACK routines. User-defined functions do not
benefit from this level of efficiency, and can impair computational
performance more than one would think (at least, more that I thought).

In summary, the functions you find in this repo replace the
computationally expensive, non-efficient parts of GP modeling. I hope
you find them helpful.
