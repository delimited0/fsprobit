library("CVXR")

# E step methods
library(tmvtnorm)
library(tmg)
library(MomTrunc)


#' Compute per observation 1st and 2nd moments for multinomial probit
#' Input can be subset of observations. Elements of X, y, and constraints must
#' correspond to same observations.
#' @importFrom epmgpr moments
#' @param beta coefficient
#' @param Sigma utility covariance
#' @param x m x p matrix, covariate
#' @param y factor response
#' @param A constraint matrix
#' @return list of approximate means and covariances
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
      # AXbeta = A %*% Xbeta
      AXbeta = vec_relative_choice(y, AXbeta)
      # browser()
      utility_moments = epmgpr::moments(
        lb = -AXbeta,
        ub = ub,
        mu = rep(0, m),
        # Sigma = A %*% tcrossprod(Sigma, A)
        Sigma = sandwich_choice(y, Sigma)
      )

      # transform back
      # utility_moments$mu = A %*% utility_moments$mu + Xbeta
      # utility_moments$Sigma = A %*% tcrossprod(utility_moments$Sigma, A)
      utility_moments$mu = vec_relative_choice(y, utility_moments$mu) + Xbeta
      utility_moments$Sigma = sandwich_choice(y, utility_moments$Sigma)
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

mnp_epmnp_moments = function(Xbeta, Sigma, y)
{
  choice_index = as.integer(y) - 1
  epmnp(as.vector(Xbeta), Sigma, choice_index)
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
      tmg::rtmg(
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
      tmg::rtmg(
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

mnp_momtrunc_moments = function(Xbeta, Sigma, y, A)
{
  m = nrow(Xbeta)
  base_choice = levels(y)[1]

  kappa = rep(1, nrow(Xbeta))

  if (y == base_choice)
  {
    lb = rep(-Inf, m)
    ub = rep(0, m)


    utility_moments = MomTrunc::meanvarTMD(
      lower = lb,
      upper = ub,
      mu = Xbeta,
      Sigma = Sigma,
      lambda=0,
      tau=0,
      dist = "normal"
    )

    utility_moments$mu = utility_moments$mean
    utility_moments$Sigma = utility_moments$varcov
  }
  else
  {
    ub = rep(Inf, m)
    # transform to axis aligned
    AXbeta = A %*% Xbeta
    lb = -AXbeta

    utility_moments = MomTrunc::meanvarTMD(
      lower = lb,
      upper = ub,
      mu = rep(0, m),
      Sigma = A %*% tcrossprod(Sigma, A),
      lambda=0,
      tau=0,
      dist = "normal"
    )

    # transform back
    utility_moments$mu = A %*% utility_moments$mean + Xbeta
    utility_moments$Sigma = A %*% tcrossprod(utility_moments$varcov, A)
  }

  utility_moments$mean = NULL
  utility_moments$EYY = NULL
  utility_moments$varcov = NULL

  return(utility_moments)
}

mnp_met_moments = function(Xbeta, Sigma, y, A, n_mc)
{
  m = nrow(Xbeta)

  # initial_point = initial_mc_point(y)

  base_choice = levels(y)[1]
  if (y == base_choice)  # base case
  {
    lb = rep(-Inf, m)
    ub = rep(0, m)

    samples = TruncatedNormal::rtmvnorm(
        n=n_mc,
        mu=as.vector(Xbeta),
        sigma=Sigma,
        lb=lb,
        ub=ub,
    )
  }
  else  # all other choices
  {
    ub = rep(Inf, m)

    # transform to axis aligned
    AXbeta = A %*% Xbeta

    samples = TruncatedNormal::rtmvnorm(
      n=n_mc,
      mu=rep(0, m),
      sigma=A %*% tcrossprod(Sigma, A),
      lb=-AXbeta,
      ub=ub
    )

    # transform back
    # recycle Xbeta to add mean per sample
    samples = t( A %*% t(samples) + as.vector(Xbeta) )
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

#' Estimate a covariance matrix subject to a trace constraint
#'
#' Performs the covariance M-step by using the eigenvectors of the expected
#' sample covariance and solving for the Lagrange multiplier associated with
#' the trace constraint.
#'
#' @param E_sample_cov Expected sample covariance matrix from the E-step.
#' @param true_trace Required trace of the updated covariance matrix.
#' @return The updated covariance matrix.
Cov_newton_estimation = function(E_sample_cov, true_trace)
{
  if (!is.matrix(E_sample_cov) ||
      nrow(E_sample_cov) != ncol(E_sample_cov) ||
      nrow(E_sample_cov) == 0) {
    stop("E_sample_cov must be a non-empty square matrix")
  }
  if (any(!is.finite(E_sample_cov))) {
    stop("E_sample_cov must contain only finite values")
  }
  if (length(true_trace) != 1 ||
      !is.finite(true_trace) ||
      true_trace <= 0) {
    stop("true_trace must be a finite positive scalar")
  }

  symmetry_tol = 100 * .Machine$double.eps *
    max(1, max(abs(E_sample_cov)))
  if (max(abs(E_sample_cov - t(E_sample_cov))) > symmetry_tol) {
    stop("E_sample_cov must be symmetric")
  }

  eigen_decomp = eigen(E_sample_cov, symmetric = TRUE)
  s_vec = eigen_decomp$values
  eigen_tol = 100 * .Machine$double.eps * max(1, max(abs(s_vec)))
  if (any(s_vec <= eigen_tol)) {
    stop("E_sample_cov must be positive definite")
  }

  sigma_at = function(lambda) {
    2 * s_vec / (1 + sqrt(1 + 4 * lambda * s_vec))
  }
  trace_difference = function(lambda) {
    sum(sigma_at(lambda)) - true_trace
  }

  trace_at_zero = sum(s_vec)
  trace_tol = sqrt(.Machine$double.eps) * max(1, true_trace)
  if (abs(trace_at_zero - true_trace) <= trace_tol) {
    sigma_vec = s_vec
  } else {
    if (trace_at_zero > true_trace) {
      lower = 0
      upper = 1
      while (trace_difference(upper) > 0) {
        upper = upper * 2
        if (!is.finite(upper)) {
          stop("Could not bracket the covariance trace multiplier")
        }
      }
    } else {
      upper = 0
      lower_limit = -1 / (4 * max(s_vec))
      lower = lower_limit * (1 - sqrt(.Machine$double.eps))
      if (trace_difference(lower) < 0) {
        stop(paste0(
          "No covariance update on the specified solution branch has trace ",
          true_trace
        ))
      }
    }

    lambda = uniroot(
      trace_difference,
      interval = c(lower, upper),
      tol = .Machine$double.eps^0.75
    )$root
    sigma_vec = sigma_at(lambda)
  }

  Sigma_new = eigen_decomp$vectors %*%
    (sigma_vec * t(eigen_decomp$vectors))
  Sigma_new = (Sigma_new + t(Sigma_new)) / 2

  return(Sigma_new)
}

Prec_cvx_estimation = function(
    E_sample_cov,
    true_trace,
    scale_iden_method,
    topleft_value = 1
)
{
  m = dim(E_sample_cov)[1]
  Prec = Variable(c(m, m), PSD=TRUE)
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
  cvx_solution = psolve(prob)
  Precision_new = cvx_solution$getValue(Prec)

  return(Precision_new)
}
