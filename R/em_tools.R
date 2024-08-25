library("CVXR")

# E step methods
library(tmvtnorm)
library(tmg)

#' Compute per observation 1st and 2nd moments for multinomial probit
#' Input can be subset of observations. Elements of X, y, and constraints must
#' correspond to same observations.
#' @param beta coefficient
#' @param Sigma utility covariance
#' @param x m x p matrix, covariate
#' @param y factor response
#' @param A constraint matrix
mnp_ep_moments = function(Xbeta, Sigma, y, A, transform = FALSE)
{
  m = nrow(Xbeta)
  base_choice = levels(y)[1]

  if (transform)
  {
    if (y == base_choice)
    {
      lb = rep(-Inf, m)
      ub = rep(0, m)

      utility_moments = epmgpr::moments(
        lb = lb,
        ub = ub,
        mu = Xbeta,
        Sigma = Sigma
      )
    }
    else
    {
      lb = rep(0, m)
      ub = rep(Inf, m)

      # transform to axis aligned
      AXbeta = A %*% Xbeta
      utility_moments = epmgpr::moments(
        lb = -AXbeta,
        ub = ub,
        mu = rep(0, m),
        Sigma = A %*% tcrossprod(Sigma, A)
      )

      # transform back
      utility_moments$mu = A %*% utility_moments$mu + Xbeta
      utility_moments$Sigma = A %*% tcrossprod(utility_moments$Sigma, A)
    }
  }
  else
  {
    if (y == base_choice)  # base case
    {
      lb = rep(-Inf, m)
      ub = rep(0, m)

      # in the base case we have axis aligned problem
      # utility_moments = epmgpr::moments(lb, ub, Xbeta, Sigma)
      utility_moments = moments2(Xbeta, Sigma, lb, ub, A)
    }
    else  # all other choices
    {
      lb = rep(0, m)
      ub = rep(Inf, m)
      utility_moments = moments2(Xbeta, Sigma, lb, ub, A)
    }
  }

  return(utility_moments)
}

mnp_hmc_moments = function(Xbeta, Precision, y, A, n_mc)
{
  m = nrow(Xbeta)

  initial_point = initial_mc_point(y)

  base_choice = levels(y)[1]
  if (y == base_choice)  # base case
  {
    lb = rep(-Inf, m)
    ub = rep(0, m)

    samples =
      rtmg(
        n=n_mc,
        M=Precision,
        r=as.vector(Precision %*% Xbeta),
        initial = initial_point,
        f = -A,
        g = ub
      )
  }
  else  # all other choices
  {
    lb = rep(0, m)
    ub = rep(Inf, m)

    samples =
      rtmg(
        n=n_mc,
        M=Precision,
        r=as.vector(Precision %*% Xbeta),
        initial = initial_point,
        f = A,
        g = lb
      )
  }
  moments = list(mu = colMeans(samples), Sigma = cov(samples))

  return(moments)
}


Prec_newton_estimation = function(
    E_sample_cov,
    true_trace,
    max_newton_iter,
    newton_tol
)
{
  eigen_decomp = eigen(E_sample_cov)

  eigen_min = eigen_decomp$values[length(eigen_decomp$values)]

  y = eigen_min - .1

  # Newton iterations to optimize M step bound
  k = 1
  dy = Inf
  while (k < max_newton_iter && abs(dy) > newton_tol) {
    s_vec = eigen_decomp$values - y
    # print(paste0("svec:", s_vec))
    fv = sum(1 / s_vec) - true_trace
    df = sum(1 / s_vec^2)
    y_new = y - (fv / df)
    dy = y_new - y
    y = y_new
    k = k+1

    # newton_history[k, iter] = y
  }
  s_vec = (eigen_decomp$values - y)

  Precision_new =
    eigen_decomp$vectors %*% diag(1 / s_vec) %*% t(eigen_decomp$vectors)
  # Sigma_new =
  # eigen_decomp$vectors %*% diag(s_vec) %*% t(eigen_decomp$vectors)

  return(Precision_new)
  # return(Sigma_new)
}

Prec_cvx_estimation = function(
    E_sample_cov,
    true_trace,
    scale_iden_method,
    topleft_value = 1
)
{
  m = dim(E_sample_cov)[1]
  Prec = Variable(m, m, PSD=TRUE)
  objective = log_det(Prec) - matrix_trace(Prec %*% E_sample_cov)
  if (scale_iden_method == "topleft")
  {
    constr = list(
      Prec[1,1] == topleft_value
    )
  }
  else if (scale_iden_method == "trace")
  {
    constr = list(
      matrix_trace(Prec) == true_trace
    )
  }
  prob = Problem(Maximize(objective), constr)
  cvx_solution = CVXR::psolve(prob)
  Precision_new = cvx_solution$getValue(Prec)

  return(Precision_new)
}
