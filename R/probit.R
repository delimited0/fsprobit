# source('constraint_tools.R')
# source('em_tools.R')
library(doRNG)
library(doParallel)
library(tictoc)

#' precision parameterized, direct opt
#' n is number of observations
#' m is number of choices
#' p is covariate dimension
#' @param X \eqn{n \times m \times p} array, covariates
#' @param Y \eqn{n \times 1} factor vector of observed choices
mnp_probit = function(
    X, Y,
    true_trace, n_choices = nlevels(Y),
    Sigma_init = NULL, beta_init = NULL,
    tol = .1,
    newton_tol = 1e-2, max_newton_iter = 15,
    max_iter = 15,
    verbose = 0,
    shift_iden_method = "ref",
    scale_iden_method = "topleft",
    topleft_value = 1,
    E_method = "EP",
    n_mc = 25,
    E_sample_rate = 1,
    true_beta = NULL, true_Sigma = NULL,
    M_method = "Newton",
    update_beta = TRUE,
    conv_metric = "precision",
    penalty = NULL,
    record_history = FALSE,
    nugget = 0,
    transform=FALSE
)
{
  # boilerplate set up ----
  # if (verbose)
  #   options(progressr.enable = TRUE)

  tic()
  n_obs = length(Y)
  p = dim(X)[3]
  # p = 1

  Ts = NULL

  # identified constraints
  if (shift_iden_method == "ref")
  {
    # constraints = lapply(Y, function(y) { utility_shift_constraints(y) })
    constraints = utility_shift_constraints(n_choices)
    m = n_choices-1
    base_choice = levels(Y)[1]
  }
  else
  {
    constraints = lapply(Y, function(y) { utility_constraints_zero(y) })
    m = n_choices
    Ts = matrix(-1/(m-1), nrow = m, ncol=m)
    diag(Ts) = 1
  }

  # initial parameter estimates
  if (is.null(Sigma_init))
    Sigma = diag(m)
  else
    Sigma = Sigma_init
  Precision = solve(Sigma)

  if (is.null(beta_init))
    beta = as.matrix(rep(0, p))
  else
    beta = beta_init

  if (is.null(penalty))
    penalty = 2*m

  m_step_bound = -1e6

  iter = 1
  # dSigma = Inf

  dmetric = Inf

  if (record_history)
  {
    llik = rep(NA, max_iter)
    Sigma_history = array(NA, dim = c(max_iter+1, m, m))
    Prec_history = array(NA, dim = c(max_iter+1, m, m))
    beta_history = matrix(NA, nrow = max_iter+1, ncol = p)

    Sigma_history[1, , ] = Sigma
    Prec_history[1, , ] = Precision
    beta_history[1, ] = beta
  }
  else
  {
    llik = NULL
    Sigma_history = NULL
    Prec_history = NULL
    beta_history = NULL
  }

  # gls_a_history = rep(NA, max_iter)
  # gls_b_history = rep(NA, max_iter)ß

  # newton_history = matrix(NA, nrow = max_newton_iter, ncol = max_iter)

  infs = rep(Inf, m)
  zeros = rep(0, m)

  # Loop body ----
  while (dmetric > tol & iter <= max_iter)
  {
    # obs_moments = vector(mode = "list", length = n_obs)
    E_sample_cov = matrix(0, nrow = m, ncol = m)
    gls_a = 0
    gls_b = 0

    if (E_sample_rate == 1)
      obs_set = 1:n_obs
    else
      obs_set = sample(1:n_obs, size = floor(E_sample_rate * n_obs), replace = FALSE)

    # per observation computations ----
    obs_moments =
      foreach(
        i = obs_set
        # .errorhandling = 'pass'
        # .inorder = TRUE
        # .combine = "list"
      ) %dorng%
      {
        # print(paste0("doing obs ", i))
        beta_e = beta
        Sigma_e = Sigma

        if (!is.null(true_beta))
          beta_e = true_beta
        if (!is.null(true_Sigma))
          Sigma_e = true_Sigma

        Xbeta = X[i, , ] %*% beta_e
        y = Y[i]

        A = constraints[[y]]

        # browser()
        if (E_method == 'EP') {
          # browser(expr = {i == 463})
          utility_moments = mnp_ep_moments(Xbeta, Sigma_e, y, A, transform)
        }
        else if (E_method == "HMC") {
          # browser(expr = {i == 45})
          utility_moments = mnp_hmc_moments(Xbeta, Precision, y, A, n_mc)
        }
        else if (E_method == "LINESS")
          utility_moments = mnp_ess_moments(Xbeta, Sigma_e, y, A, n_mc)
        else if (E_method == "Gibbs")
          utility_moments = mnp_gibbs_moments(Xbeta, Precision, y, A, n_mc)
        else
          stop('Moments must be one of EP, HMC, Gibbs, or LINESS')

        # accumulate m step quantities
        # tXPrecision = crossprod(X[i, , ], Precision)
        # gls_a = tXPrecision %*% X[i, , ]
        # gls_b = tXPrecision %*% moments$mu
        #
        # # E[Z | -]
        # E_sample_cov = moments$Sigma + tcrossprod(moments$mu - Xbeta)
        #
        # return(
        #   list(
        #     gls_a = gls_a,
        #     gls_b = gls_b,
        #     E_sample_cov = E_sample_cov
        #   )
        # )

        # print(typeof(utility_moments))

        # browser()

        # browser(
        #   expr = {typeof(utility_moments$mu) == "double"}
        # )

        # result_list = list(
        #   mu = utility_moments$mu,
        #   Sigma = utility_moments$Sigma
        # )

        # print(typeof(result_list))

        # return(result_list)
        return(utility_moments)
      }

    # propagate updated values to outer environment
    # gls_a = E_obs_quantities$gls_a
    # gls_b = E_obs_quantities$gls_b
    # E_sample_cov = E_obs_quantities$E_sample_cov
    # obs_moments = E_obs_quantities

    # browser()

    # record old parameters
    beta_old = beta
    Sigma_old = Sigma
    Precision_old = Precision
    m_step_bound_old = m_step_bound

    # M step ----
    # beta estimation
    for (i in 1:length(obs_moments))
    {
      obs_idx = obs_set[i]

      tXPrecision = crossprod(X[obs_idx, ,], Precision)
      gls_a = gls_a + tXPrecision %*% X[obs_idx, ,]
      gls_b = gls_b + tXPrecision %*% obs_moments[[i]]$mu
    }
    beta_new = solve(gls_a, gls_b)

    # E[Q | -]
    for (i in 1:length(obs_moments))
    {
      obs_idx = obs_set[i]
      Xbeta = X[obs_idx, ,] %*% beta_new
      E_sample_cov =
        E_sample_cov + obs_moments[[i]]$Sigma +
        tcrossprod(obs_moments[[i]]$mu - Xbeta)
    }
    E_sample_cov = E_sample_cov / length(obs_moments)
    # browser()

    #### precision estimation ####
    # trace constraint
    if (M_method == "Newton")
    {
      Precision_new = Prec_newton_estimation(
        # Sigma_new = Prec_newton_estimation(
        E_sample_cov = E_sample_cov,
        true_trace = true_trace,
        max_newton_iter = max_newton_iter,
        newton_tol = newton_tol
      )
      Precision_new = Precision_new + nugget*diag(m)
      Sigma_new = solve(Precision_new)
    }
    else if (M_method == "CVX")
    {
      Precision_new = Prec_cvx_estimation(
        E_sample_cov = E_sample_cov,
        true_trace = true_trace,
        scale_iden_method = scale_iden_method,
        topleft_value = topleft_value
      )
      Precision_new = Precision_new + nugget*diag(m)
      Sigma_new = solve(Precision_new)
    }
    else
    {
      stop("Invalid m step method")
    }

    m_step_bound_new = - (sum(diag(E_sample_cov %*% Precision_new)) + determinant(Precision_new)$modulus )

    # update parameters
    beta = beta_new
    Precision = Precision_new
    Sigma = Sigma_new
    m_step_bound = m_step_bound_new

    #### record history ####
    if (record_history)
    {
      llik[iter] = m_step_bound_new
      beta_history[iter+1, ] = beta_new
      Sigma_history[iter+1, , ] = Sigma_new
      Prec_history[iter+1, , ] = Precision_new
    }

    #### update convergence ####
    if (conv_metric == "precision")
    {
      param_new = c(Precision_new, beta_new)
      param_old = c(Precision_old, beta_old)
      # dmetric = max(abs(param_new - param_old) / abs(param_old))
      dmetric = max(abs(param_new - param_old))
    }
    else if (conv_metric == "covariance")
    {
      param_new = c(Sigma_new, beta_new)
      param_old = c(Sigma_old, beta_old)
      # dmetric = max(abs(param_new  - param_old) / abs(param_old))
      dmetric = max(abs(param_new - param_old))
    }
    else if (conv_metric == "mbound")
    {
      dmetric = abs(m_step_bound_new - m_step_bound_old)
      m_step_bound = m_step_bound_new
    }
    else
      stop("invalid convergence metric")

    if (iter %% verbose == 0)
    {
      # print(paste0("EM iteration: ", iter))
      print(paste0("EM iteration: ", iter, ", M step bound: ", m_step_bound_new))
    }

    #### update parameters ####
    iter = iter + 1
    beta = beta_new
    Sigma = Sigma_new
    # Precision = Precision_new
  }
  elapsed = toc(quiet = (!(verbose != 0)))

  if (verbose != 0)
    print(paste0("Total iterations: ", iter))

  if (record_history) {
    Sigma_history = Sigma_history[1:(iter), , ]
    Prec_history = Prec_history[1:(iter), , ]
    beta_history = beta_history[1:(iter), ]
    llik = llik[1:(iter-1)]
  }

  return(list(
    Sigma = Sigma,
    Precision = Precision,
    beta = beta,
    m_step_bound = m_step_bound,
    iters = iter,
    Sigma_history = Sigma_history,
    Prec_history = Prec_history,
    beta_history = beta_history,
    llik = llik,
    # newton_history = newton_history[, 1:(iter-1)],
    # gls_a_history = gls_a_history[1:iter],
    # gls_b_history = gls_b_history[1:iter],
    E_method = E_method,
    M_method = M_method,
    E_sample_rate = E_sample_rate,
    shift_iden_method = shift_iden_method,
    scale_iden_method = scale_iden_method,
    time = elapsed$toc - elapsed$tic
  ))
}
