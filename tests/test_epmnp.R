library(fsprobit)
library(epmgpr)

set.seed(1)

check_epmnp_case = function(m, choice_index, tol = 1e-6)
{
  M = matrix(rnorm(m * m), nrow = m)
  Sigma = crossprod(M) / m + diag(m) * 0.5
  mu = rnorm(m)
  constraints = utility_shift_constraints(m + 1)

  if (choice_index == 0) {
    A = constraints[[1]]
    lb = rep(-Inf, m)
    ub = rep(0, m)
  } else {
    A = constraints[[choice_index + 1]]
    lb = rep(0, m)
    ub = rep(Inf, m)
  }

  dense = epmgpr::moments2(mu, Sigma, lb, ub, A)
  sparse = epmnp(mu, Sigma, choice_index)

  stopifnot(length(sparse$mu) == m)
  stopifnot(all(dim(sparse$Sigma) == c(m, m)))
  stopifnot(max(abs(sparse$Sigma - t(sparse$Sigma))) < tol)
  stopifnot(max(abs(as.vector(dense$mu) - as.vector(sparse$mu))) < tol)
  stopifnot(max(abs(dense$Sigma - sparse$Sigma)) < tol)
}

for (choice_index in 0:2) {
  check_epmnp_case(2, choice_index)
}

for (choice_index in 0:5) {
  check_epmnp_case(5, choice_index)
}

for (choice_index in c(0, 1, 4, 10)) {
  check_epmnp_case(10, choice_index)
}

mu = c(0, 0)
Sigma = diag(2)

stopifnot(inherits(try(epmnp(mu, Sigma, -1), silent = TRUE), "try-error"))
stopifnot(inherits(try(epmnp(mu, Sigma, 3), silent = TRUE), "try-error"))
stopifnot(inherits(try(epmnp(mu, diag(3), 0), silent = TRUE), "try-error"))
