library(epmgpr)

#' utility constraints under reference identification for shift invariance
#' constraints are on n_choices-1 non-base case choices
#' base case is assumed to be first factor level
#'
#' @param n_choices number of choices
#' @return list of constraint matrices implied by each possible choice
utility_shift_constraints = function(n_choices)
{
  constraints = vector(mode = "list", length = n_choices)

  # base case
  constraints[[1]] = diag(n_choices-1)

  for (choice in 2:n_choices) {
    idx = choice - 1
    A = -1*diag(n_choices-1)
    A[, idx] = 1
    constraints[[choice]] = A
  }

  # other cases
  # D = -1*diag(n_choices-2)
  #
  # # # choice 2
  # A = cbind(1, D)
  # A = rbind(A, c(1, rep(0, n_choices-2)))
  # constraints[[2]] = A
  #
  # # choices 3...n_choices-1
  # if (n_choices > 3)
  # {
  #   for (idx in 3:(n_choices-1))
  #   {
  #     A = cbind(D[, 1:(idx-2)], rep(1, n_choices-2), D[, (idx-1):(n_choices-2)])
  #     A = rbind(A, c( rep(0, idx-2), 1, rep(0, n_choices-idx) ))
  #     constraints[[idx]] = A
  #   }
  # }
  #
  # # # choice n_choices
  # A = cbind(D, 1)
  # A = rbind(A, c(rep(0, n_choices-2), 1))
  # constraints[[n_choices]] = A
  #
  return(constraints)
}

#' @param constraints list of constraint matrices, with first element corresponding to base choice
#' @importFrom epmgpr moments2
#' @param Sigma \eqn{m-1 \times m-1} covariance matrix
#' @return list of approximate means and covariances
utility_ep_approximations = function(constraints, Sigma)
{
  n_choices = length(constraints)

  approximations = vector(mode = "list", length = n_choices)
  zeros = rep(0, n_choices-1)
  infs = rep(Inf, n_choices-1)

  mu = rep(0, n_choices-1)

  # base choice approximation
  lb = -infs
  ub = zeros
  mom = moments2(mu, Sigma, lb, ub, constraints[[1]])
  approximations[[1]] = list(mu = mom$mu, Sigma = mom$Sigma)

  # non base choice approximations
  lb = zeros
  ub = infs
  for (i in 2:n_choices)
  {
    mom = moments2(mu, Sigma, lb, ub, constraints[[i]])
    approximations[[i]] = list(mu = mom$mu, Sigma = mom$Sigma)
  }

  return(approximations)
}

#' @param y factor, the individual's observed choice
initial_mc_point = function(y) {
  m = nlevels(y)
  idx = as.numeric(y)

  # base class case
  if (idx == 1) {
    util = rep(-1, m-1)
    # util = rtnorm(n = m-1, mu = 0, sd = 1, lb = -Inf, ub = 0)
  }
  else if (idx == m) {
    util = c(rep(-1, m-2), 1)

  }
  else {
    util = c(rep(-1, idx-2), 1, rep(-1, m-idx))
  }
  return(util)
}
