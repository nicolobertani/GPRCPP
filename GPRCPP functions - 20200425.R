library(RcppArmadillo)

# SE - Function and derivatives ----------------------------------------------------

Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_k_se.cpp")
Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_d_k_se_l.cpp")
Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_d_k_se_s.cpp")


# Periodic - Function and derivatives -------------------------------------

Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_k_per.cpp")
Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_d_k_per_l.cpp")
Rcpp::sourceCpp("~/Dropbox (INSEAD)/Crime Modeling/CPP code/rcpp_d_k_per_s.cpp")


d.k.per.l <- function(X1, X2 = NULL, l = 1, s = 1, p = 1) {
  k <- function(x, y) {
    s ^ 2 * exp(- 2 * sin(pi * norm(x - y, "2") / p) ^ 2 / (l ^ 2)) * (4 * sin(pi * norm(x - y, "2") / p) ^ 2) / (l ^ 3)   # changed here
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
# sum(k.per(X.m, l = 2 + 1e-6, s = 4) - k.per(X.m, l = 2 - 1e-6, s = 4)) / 2e-6
# sum(d.k.per.l(X.m, l = 2, s = 4))
# sum(k.per(X.m, X.n, l = 2 + 1e-6) - k.per(X.m, X.n, l = 2 - 1e-6)) / 2e-6
# sum(d.k.per.l(X.m, X.n, l = 2))

d.k.per.s <- function(X1, X2 = NULL, l = 1, s = 1, p = 1) {
  k <- function(x, y) {
    2 * s * exp(- 2 * sin(pi * norm(x - y, "2") / p) ^ 2 / (l ^ 2))  # changed here
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
# sum(k.per(X.m, s = 8 + 1e-6) - k.per(X.m, s = 8 - 1e-6)) / 2e-6
# sum(d.k.per.s(X.m, s = 8))
# sum(k.per(X.m, X.n, s = 8 + 1e-6) - k.per(X.m, X.n, s = 8 - 1e-6)) / 2e-6
# sum(d.k.per.s(X.m, X.n, s = 8))


# Linear - Function and derivatives -----------------------------------

k.lin <- function(X1, X2 = NULL, l = 1, s = 1) {
  k <- function(x, y) {
    s ^ 2 * t(x - l) %*% (y - l)
  }
  x1 <- as.matrix(X1)
  if (!length(X2)) {
    output <- matrix(0, nrow = nrow(x1), ncol = nrow(x1))
    sapply(seq(nrow(x1)), function(i) {
      sapply(seq(nrow(x1))[seq(nrow(x1)) > i], function(j) {
        output[i, j] <<- k(x1[i, ], x1[j, ]) 
      })})
    output <- output + t(output)
    diag(output) <- sapply(seq(nrow(x1)), function(j) {k(x1[j, ], x1[j, ])})
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

d.k.lin.l <- function(X1, X2 = NULL, l = 1, s = 1) {
  k <- function(x, y) {
    - s ^ 2 * (sum(x - l) + sum(y - l))   # changed here
  }
  x1 <- as.matrix(X1)
  if (!length(X2)) {
    output <- matrix(0, nrow = nrow(x1), ncol = nrow(x1))
    sapply(seq(nrow(x1)), function(i) {
      sapply(seq(nrow(x1))[seq(nrow(x1)) > i], function(j) {
        output[i, j] <<- k(x1[i, ], x1[j, ]) 
      })})
    output <- output + t(output)
    diag(output) <- sapply(seq(nrow(x1)), function(j) {k(x1[j, ], x1[j, ])}) # changed here
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
# sum(k.lin(X.m, l = -8 + 1e-6) - k.lin(X.m, l = -8 - 1e-6)) / 2e-6
# sum(d.k.lin.l(X.m, l = -8))
# sum(k.lin(X.m, X.n, l = -8 + 1e-6) - k.lin(X.m, X.n, l = -8 - 1e-6)) / 2e-6
# sum(d.k.lin.l(X.m, X.n, l = -8))

d.k.lin.s <- function(X1, X2 = NULL, l = 1, s = 1) {
  k <- function(x, y) {
    2 * s * t(x - l) %*% (y - l)  # changed here
  }
  x1 <- as.matrix(X1)
  if (!length(X2)) {
    output <- matrix(0, nrow = nrow(x1), ncol = nrow(x1))
    sapply(seq(nrow(x1)), function(i) {
      sapply(seq(nrow(x1))[seq(nrow(x1)) > i], function(j) {
        output[i, j] <<- k(x1[i, ], x1[j, ]) 
      })})
    output <- output + t(output)
    diag(output) <- sapply(seq(nrow(x1)), function(j) {k(x1[j, ], x1[j, ])}) # changed here
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
# sum(k.lin(X.m, s = 8 + 1e-6) - k.lin(X.m, s = 8 - 1e-6)) / 2e-6
# sum(d.k.lin.s(X.m, s = 8))
# sum(k.lin(X.m, X.n, s = 8 + 1e-6) - k.lin(X.m, X.n, s = 8 - 1e-6)) / 2e-6
# sum(d.k.lin.s(X.m, X.n, s = 8))


# Titsias 2009 - functions ------------------------------------------------

VLB.tits <- function(X.m, l.tits, s.tits, sigma.sq.tits, jitter.v = 0) {
  K.mm <<- k.se(X.m, l = l.tits, s = s.tits) + diag(jitter.v, m.size)
  U.mm <<- chol(K.mm)
  K.mn <<- k.se(X.m, X.sample, l = l.tits, s = s.tits)
  UK.mn <- forwardsolve(t(U.mm), K.mn)
  Q.nn <- t(UK.mn) %*% UK.mn
  Sigma.U <- chol(diag(sigma.sq.tits, sample.size) + Q.nn)
  inv.Sigma.U.y <- forwardsolve(t(Sigma.U), y.sample)
  output <- - sample.size / 2 * log(2 * pi) - 
    sum(log(diag(t(Sigma.U)))) -
    t(inv.Sigma.U.y) %*% inv.Sigma.U.y / 2 -
    sum(sapply(X.sample, function(x.i) {k.se(x.i, x.i, l = l.tits, s = s.tits)})) / (2 * sigma.sq.tits) +
    sum(diag(Q.nn)) / (2 * sigma.sq.tits)
  return(output)
}
step.l.tits <- function() {
  d.K.mm <- d.k.se.l(X.m, l = par.tits[1], s = par.tits[2])
  d.K.mn <- d.k.se.l(X.m, X.sample, l = par.tits[1], s = par.tits[2])
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn - 
    t(K.mn) %*% inv.K.mm %*% d.K.mm %*% inv.K.mm %*% K.mn + 
    t(K.mn) %*% inv.K.mm %*% d.K.mn
  output <- sum(sapply(seq(sample.size), function(i) {M[i, ] %*% d.Q.nn[, i]})) / 2 + 
    sum(diag(d.Q.nn)) / (2 * par.tits[3])
  return(output)
}
step.s.tits <- function() {
  d.K.mm <- d.k.se.s(X.m, l = par.tits[1], s = par.tits[2])
  d.K.mn <- d.k.se.s(X.m, X.sample, l = par.tits[1], s = par.tits[2])
  tr.d.K.nn <- sum(sapply(X.sample, d.k.se.s, l = par.tits[1], s = par.tits[2]))
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn - 
    t(K.mn) %*% inv.K.mm %*% d.K.mm %*% inv.K.mm %*% K.mn + 
    t(K.mn) %*% inv.K.mm %*% d.K.mn
  output <- sum(sapply(seq(sample.size), function(i) {M[i, ] %*% d.Q.nn[, i]})) / 2 -
    tr.d.K.nn / (2 * par.tits[3]) +
    sum(diag(d.Q.nn)) / (2 * par.tits[3])
  return(output)
}
step.sigma.sq.tits <- function() {
  output <- - sum(diag(inv.K)) / 2 + 
    t(a) %*% a / 2 + 
    sum(sapply(X.sample, function(x.i) {k.se(x.i, x.i, l = par.tits[1], s = par.tits[2])})) / (2 * par.tits[3] ^ 2) - 
    sum(diag(t(K.mn) %*% inv.K.mm %*% K.mn)) / (2 * par.tits[3] ^ 2)
  return(output)
}
step.x.tits <- function(r, c = 1) {
  d.K.mm <- d.k.se.x(X.m, l = par.tits[1], s = par.tits[2], r = r, c = c)
  d.K.mn <- d.k.se.x(X.m, X.sample, l = par.tits[1], s = par.tits[2], r = r, c = c)
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn - 
    t(K.mn) %*% inv.K.mm %*% d.K.mm %*% inv.K.mm %*% K.mn + 
    t(K.mn) %*% inv.K.mm %*% d.K.mn
  output <- sum(sapply(seq(sample.size), function(i) {M[i, ] %*% d.Q.nn[, i]})) / 2 + 
    sum(diag(d.Q.nn)) / (2 * par.tits[3])
  return(output)
}


# Hensman et alii 2013 - functions ----------------------------------------

# # CHECKS
# Lambda_i <- function(i) {(inv.K.mm %*% K.mn[, i]) %*% (t(K.mn[, i]) %*% inv.K.mm) / sigma.sq.svi}
# all.equal(Lambda, Reduce('+', lapply(seq(batch.size), Lambda_i)) + inv.K.mm)
# L3 <- function() {
#   sum(sapply(seq(batch.size), function(i) {
#     dnorm(y.svi[i], t(K.mn[, i]) %*% (inv.K.mm %*% m.svi), sqrt(sigma.sq.svi), log = T)
#   })) -
#     sum(diag(K.tilde)) / (2 * sigma.sq.svi) -
#     sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% (Lambda - inv.K.mm)[, j]})) / 2 -
#     (sum(diag(inv.K.mm %*% S.svi)) + t(m.svi) %*% (inv.K.mm %*% m.svi) - m.size + log(det(K.mm)) - log(det(S.svi))) / 2
# }
# # END OF CHECKS

L3 <- function(sigma.sq.svi) {
  sum(sapply(seq(batch.size), function(i) {
    dnorm(y.svi[i], t(K.mn[, i]) %*% (inv.K.mm %*% m.svi), sqrt(sigma.sq.svi), log = T)
  })) -
    .5 / sigma.sq.svi * sum(diag(K.tilde)) -
    .5 * sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% (Lambda - inv.K.mm)[, j]})) -
    .5 * (sum(sapply(seq(m.size), function(j) {inv.K.mm[j, ] %*% S.svi[, j]})) + 
            t(m.svi) %*% (inv.K.mm %*% m.svi) - m.size) - 
    sum(log(diag(U.mm))) + sum(log(diag(U.S.svi)))
}

covariance.eval <- function(which.K) {
  Reduce('+', mclapply(covariance.calls[[which.K]], eval))
}

update.L3 <- function(covariance.calls, sigma.sq.svi, jitter.v = 1e-2) {
  K.mm <<- covariance.eval('K.mm')
  U.mm <<- chol(K.mm)
  inv.K.mm <<- chol2inv(U.mm)
  K.mn <<- covariance.eval('K.mn')
  K.nn <<- covariance.eval('K.nn')
  K.tilde <<- K.nn - t(K.mn) %*% inv.K.mm %*% K.mn
  Lambda <<- inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm / sigma.sq.svi + inv.K.mm
  U.S.svi <<- chol(S.svi)
  return(L3(sigma.sq.svi))
}

# without momentum
update.variational.parameters <- function(sigma.sq.svi, step.size = 1e-2) {
  if (length(step.size) == 1) step.size <- rep(step.size, 2)
  inv.S <- chol2inv(U.S.svi)
  updates <- list(
    inv.K.mm %*% (K.mn %*% y.svi) / sigma.sq.svi - inv.S %*% m.svi,
    inv.S / 2 - Lambda / 2
  )
  theta.new <- list(
    inv.S %*% m.svi + step.size[1] * updates[[1]],
    - inv.S / 2 + step.size[2] * updates[[2]]
  )
  S.new <- chol2inv(chol(- 2 * theta.new[[2]]))
  S.new <- chol2inv(chol(- 2 * theta.new[[2]]))
  m.new <- S.new %*% theta.new[[1]]
  output <- list(c(m.new), S.new)
  return(output)
}
# # with momentum
# update.variational.parameters <- function(m, S, step.size = 1, momentum = 0, sigma.sq.svi) {
#   if (length(step.size) == 1) step.size <- rep(step.size, 2)
#   if (length(momentum) == 1) momentum <- rep(momentum, 2)
#   inv.S <- chol2inv(chol(S))
#   if ('old.updates' %in% ls()) {
#     updates <- list(
#       (1 - momentum[1]) * (inv.K.mm %*% (K.mn %*% y.svi) / sigma.sq.svi - inv.S %*% m) + momentum * old.updates[[1]],
#       (1 - momentum[2]) * (inv.S / 2 - Lambda / 2) + momentum * old.updates[[2]]
#     )} else {
#       updates <- list(
#         inv.K.mm %*% (K.mn %*% y.svi) / sigma.sq.svi - inv.S %*% m,
#         inv.S / 2 - Lambda / 2
#       )}
#   old.updates <<- updates
#   theta.new <- list(
#     inv.S %*% m + step.size[1] * updates[[1]],
#     - inv.S / 2 + step.size[2] * updates[[2]]
#   )
#   S.new <- chol2inv(chol(- 2 * theta.new[[2]]))
#   m.new <- S.new %*% theta.new[[1]]
#   output <- list(c(m.new), S.new)
#   return(output)
# }

step.k.sigma.sq.svi <- function(sigma.sq.svi = 1) {
  z <- y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)
  output <- - batch.size / (2 * sigma.sq.svi) + 
    (t(z) %*% z + sum(diag(K.tilde)) + 
       sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% (sigma.sq.svi * (Lambda - inv.K.mm))[, j]}))) / (2 * sigma.sq.svi ^ 2)
  return(output)
}

# Hensman - SE+ ------------------------------------------------------------

covariance.calls <- lapply(list(
  K.mm = "k.se(X.m, l = l1, s = s1) + k.se(X.m, l = l2, s = s2) + diag(jitter.v, m.size)",
  K.mn = "k.se(X.m, X.n, l = l1, s = s1) + k.se(X.m, X.n, l = l2, s = s2)",
  K.nn = "k.se(X.n, l = l1, s = s1) + k.se(X.n, l = l2, s = s2)"
), function(i) parse(text = i))

step.k.se.l.svi <- function(X.m, X.n, relevant.l, relevant.s, sigma.sq.svi) {
  d.K.mm <- d.k.se.l(X.m, l = relevant.l, s = relevant.s)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.se.l(X.m, X.n, l = relevant.l, s = relevant.s)
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    sum(diag(-d.Q.nn)) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}

step.k.se.s.svi <- function(X.m, X.n, relevant.l, relevant.s, sigma.sq.svi) {
  d.K.mm <- d.k.se.s(X.m, l = relevant.l, s = relevant.s)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.se.s(X.m, X.n, l = relevant.l, s = relevant.s)
  tr.d.K.nn <- 2 * relevant.s * nrow(as.matrix(X.n)) # diag is constant
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    (tr.d.K.nn - sum(diag(d.Q.nn))) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}


# Hensman - Per+ -----------------------------------------------------------

step.k.per.l.svi <- function(X.m, X.n, relevant.l, relevant.s, sigma.sq.svi, period) {
  d.K.mm <- d.k.per.l(X.m, l = relevant.l, s = relevant.s, p = period)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.per.l(X.m, X.n, l = relevant.l, s = relevant.s, p = period)
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    sum(diag(-d.Q.nn)) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}

step.k.per.s.svi <- function(X.m, X.n, relevant.l, relevant.s, sigma.sq.svi, period) {
  d.K.mm <- d.k.per.s(X.m, l = relevant.l, s = relevant.s, p = period)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.per.s(X.m, X.n, l = relevant.l, s = relevant.s, p = period)
  tr.d.K.nn <- 2 * relevant.s * nrow(as.matrix(X.n))
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    (tr.d.K.nn - sum(diag(d.Q.nn))) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}


# Hensman - Lin+ -----------------------------------------------------------

step.k.lin.l.svi <- function(X.m, X.n, relevant.l, relevant.s, sigma.sq.svi) {
  d.K.mm <- d.k.lin.l(X.m, l = relevant.l, s = relevant.s)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.lin.l(X.m, X.n, l = relevant.l, s = relevant.s)
  # tr.d.K.nn <- sum(diag(d.k.lin.l(X.n, l = relevant.l, s = relevant.s)))
  tr.d.K.nn <- sum(sapply(seq(nrow(as.matrix(X.n))), function(r) {- 2 * relevant.s ^ 2 * sum(X.n[r, ] - relevant.l)}))
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    (tr.d.K.nn - sum(diag(d.Q.nn))) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}

step.k.lin.s.svi <- function(X.m, X.n, relevant.l, relevant.s, sigma.sq.svi) {
  d.K.mm <- d.k.lin.s(X.m, l = relevant.l, s = relevant.s)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.lin.s(X.m, X.n, l = relevant.l, s = relevant.s)
  # tr.d.K.nn <- sum(diag(d.k.lin.s(X.n, l = relevant.l, s = relevant.s)))
  tr.d.K.nn <- sum(sapply(seq(nrow(as.matrix(X.n))), function(r) {2 * relevant.s * sum((X.n[r, ] - relevant.l) ^ 2)}))
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    (tr.d.K.nn - sum(diag(d.Q.nn))) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}


# Hensman - SExSE ------------------------------------------------------------

step.k.sese.l.svi <- function(X.m.1, X.n.1, X.m.2 = X.m.1, X.n.2 = X.n.1, 
                              relevant.l, relevant.s, relevant.m, sigma.sq.svi) {
  d.K.mm <- d.k.se.l(X.m.1, l = relevant.l, s = relevant.s) * k.se(X.m.2, l = relevant.m)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.se.l(X.m.1, X.n.1, l = relevant.l, s = relevant.s) * k.se(X.m.2, X.n.2, l = relevant.m)
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    sum(diag(-d.Q.nn)) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}

step.k.sese.s.svi <- function(X.m.1, X.n.1, X.m.2 = X.m.1, X.n.2 = X.n.1, 
                              relevant.l, relevant.s, relevant.m, sigma.sq.svi) {
  d.K.mm <- d.k.se.s(X.m.1, l = relevant.l, s = relevant.s) * k.se(X.m.2, l = relevant.m)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.se.s(X.m.1, X.n.1, l = relevant.l, s = relevant.s) * k.se(X.m.2, X.n.2, l = relevant.m)
  # tr.d.K.nn <- sum(diag(d.k.se.s(X.n.1, l = relevant.l, s = relevant.s) * k.se(X.n.2, l = relevant.m)))
  tr.d.K.nn <- 2 * relevant.s * nrow(as.matrix(X.n.1)) # diag is constant
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    (tr.d.K.nn - sum(diag(d.Q.nn))) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}

step.k.sese.m.svi <- function(X.m.1, X.n.1, X.m.2 = X.m.1, X.n.2 = X.n.1, 
                              relevant.l, relevant.s, relevant.m, sigma.sq.svi) {
  d.K.mm <- k.se(X.m.1, l = relevant.l, s = relevant.s) * d.k.se.l(X.m.2, l = relevant.m)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- k.se(X.m.1, X.n.1, l = relevant.l, s = relevant.s) * d.k.se.l(X.m.2, X.n.2, l = relevant.m)
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    sum(diag(-d.Q.nn)) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}


# Hensman - PerxSE --------------------------------------------------------

step.k.PerSE.l.svi <- function(X.m.1, X.n.1, X.m.2 = X.m.1, X.n.2 = X.n.1, 
                               relevant.l, relevant.s, relevant.m, sigma.sq.svi, period) {
  d.K.mm <- d.k.per.l(X.m.1, l = relevant.l, s = relevant.s, p = period) * k.se(X.m.2, l = relevant.m)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.per.l(X.m.1, X.n.1, l = relevant.l, s = relevant.s, p = period) * k.se(X.m.2, X.n.2, l = relevant.m)
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    sum(diag(-d.Q.nn)) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}

step.k.PerSE.s.svi <- function(X.m.1, X.n.1, X.m.2 = X.m.1, X.n.2 = X.n.1, 
                               relevant.l, relevant.s, relevant.m, sigma.sq.svi, period) {
  d.K.mm <- d.k.per.s(X.m.1, l = relevant.l, s = relevant.s, p = period) * k.se(X.m.2, l = relevant.m)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.per.s(X.m.1, X.n.1, l = relevant.l, s = relevant.s, p = period) * k.se(X.m.2, X.n.2, l = relevant.m)
  # tr.d.K.nn <- sum(diag(d.k.per.s(X.n.1, l = relevant.l, s = relevant.s, p = period) * k.se(X.n.2, l = relevant.m)))
  tr.d.K.nn <- 2 * relevant.s * nrow(as.matrix(X.n.1)) 
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    (tr.d.K.nn - sum(diag(d.Q.nn))) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}

step.k.PerSE.m.svi <- function(X.m.1, X.n.1, X.m.2 = X.m.1, X.n.2 = X.n.1, 
                               relevant.l, relevant.s, relevant.m, sigma.sq.svi, period) {
  d.K.mm <- k.per(X.m.1, l = relevant.l, s = relevant.s, p = period) * d.k.se.l(X.m.2, l = relevant.m)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- k.per(X.m.1, X.n.1, l = relevant.l, s = relevant.s, p = period) * d.k.se.l(X.m.2, X.n.2, l = relevant.m)
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    sum(diag(-d.Q.nn)) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}


# Hensman - LinxSE -----------------------------------------------------------

step.k.LinSE.l.svi <- function(X.m.1, X.n.1, X.m.2 = X.m.1, X.n.2 = X.n.1, 
                               relevant.l, relevant.s, relevant.m, sigma.sq.svi) {
  d.K.mm <- d.k.lin.l(X.m.1, l = relevant.l, s = relevant.s) * k.se(X.m.2, l = relevant.m)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.lin.l(X.m.1, X.n.1, l = relevant.l, s = relevant.s) * k.se(X.m.2, X.n.2, l = relevant.m)
  # tr.d.K.nn <- sum(diag(d.k.lin.l(X.n.1, l = relevant.l, s = relevant.s) * k.se(X.n.2, l = relevant.m)))
  tr.d.K.nn <- sum(sapply(seq(nrow(as.matrix(X.n.1))), function(r) {- 2 * relevant.s ^ 2 * sum(X.n.1[r, ] - relevant.l)}))
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    (tr.d.K.nn - sum(diag(d.Q.nn))) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}

step.k.LinSE.s.svi <- function(X.m.1, X.n.1, X.m.2 = X.m.1, X.n.2 = X.n.1,
                               relevant.l, relevant.s, relevant.m, sigma.sq.svi) {
  d.K.mm <- d.k.lin.s(X.m.1, l = relevant.l, s = relevant.s) * k.se(X.m.2, l = relevant.m)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- d.k.lin.s(X.m.1, X.n.1, l = relevant.l, s = relevant.s) * k.se(X.m.2, X.n.2, l = relevant.m)
  # tr.d.K.nn <- sum(diag(d.k.lin.s(X.n.1, l = relevant.l, s = relevant.s) * k.se(X.n.2, l = relevant.m)))
  tr.d.K.nn <- sum(sapply(seq(nrow(as.matrix(X.n.1))), function(r) {2 * relevant.s * sum((X.n.1[r, ] - relevant.l) ^ 2)}))
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    (tr.d.K.nn - sum(diag(d.Q.nn))) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2
  return(output)
}

step.k.LinSE.m.svi <- function(X.m.1, X.n.1, X.m.2 = X.m.1, X.n.2 = X.n.1,
                               relevant.l, relevant.s, relevant.m, sigma.sq.svi) {
  d.K.mm <- k.lin(X.m.1, l = relevant.l, s = relevant.s) * d.k.se.l(X.m.2, l = relevant.m)
  d.inv.K.mm <- - inv.K.mm %*% d.K.mm %*% inv.K.mm
  d.K.mn <- k.lin(X.m.1, X.n.1, l = relevant.l, s = relevant.s) * d.k.se.l(X.m.2, X.n.2, l = relevant.m)
  d.Q.nn <-  t(d.K.mn) %*% inv.K.mm %*% K.mn + t(K.mn) %*% d.inv.K.mm %*% K.mn + t(K.mn) %*% inv.K.mm %*% d.K.mn
  d.sum.Lambda_i <- (d.inv.K.mm %*% K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% d.K.mn %*% t(K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(d.K.mn) %*% inv.K.mm + 
                       inv.K.mm %*% K.mn %*% t(K.mn) %*% d.inv.K.mm) / sigma.sq.svi 
  
  output <- t(y.svi - t(K.mn) %*% (inv.K.mm %*% m.svi)) %*% (t(d.K.mn) %*% (inv.K.mm %*% m.svi) + t(K.mn) %*% (d.inv.K.mm %*% m.svi)) / sigma.sq.svi - 
    sum(diag(-d.Q.nn)) / (2 * sigma.sq.svi) - 
    sum(sapply(seq(m.size), function(j) {S.svi[j, ] %*% d.sum.Lambda_i[, j]})) / 2 + 
    (sum(diag((inv.K.mm %*% S.svi %*% inv.K.mm - inv.K.mm) %*% d.K.mm)) - t(m.svi) %*% d.inv.K.mm %*% m.svi) / 2 
  return(output)
}








