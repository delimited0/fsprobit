#' Accelerated multinomial probit EM
#' n is number of observations
#' m is number of choices
#' p is covariate dimension
#' @param X \eqn{n \times m \times p} array, covariates
#' @param Y \eqn{n \times 1} factor vector of observed choices
mnp_probit_accelerated = function(
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
    E_method = "EPMNP",
    n_mc = 25,
    E_sample_rate = 1,
    true_beta = NULL, true_Sigma = NULL,
    M_method = "Newton",
    M_damping = 0,
    update_beta = TRUE,
    conv_metric = "precision",
    penalty = NULL,
    record_history = FALSE,
    nugget = 0,
    transform=FALSE
)
{
  # boilerplate set up ----
  tictoc::tic()
  n_obs = length(Y)
  p = dim(X)[3]

  Ts = NULL

  if (!(E_method %in% c("EPMNP", "HMC", "MET")))
    stop("mnp_probit_accelerated supports E_method = 'EPMNP', 'HMC', or 'MET'")
  if (E_method %in% c("HMC", "MET") && shift_iden_method != "ref")
    stop(paste0(E_method, " requires shift_iden_method = 'ref'"))
  if (E_method %in% c("HMC", "MET") && (!is.numeric(n_mc) || length(n_mc) != 1 ||
      is.na(n_mc) || n_mc < 2 || n_mc != as.integer(n_mc)))
    stop("n_mc must be a single integer of at least 2 for HMC or MET")

  # identified dimensions
  if (shift_iden_method == "ref")
  {
    m = n_choices-1
    base_choice = levels(Y)[1]
  }
  else
  {
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

  Y_int = as.integer(Y)

  # Loop body ----
  while (dmetric > tol & iter <= max_iter)
  {
    if (E_sample_rate == 1)
      obs_set = 1:n_obs
    else
      obs_set = sample(1:n_obs, size = floor(E_sample_rate * n_obs), replace = FALSE)

    # record old parameters
    beta_old = beta
    Sigma_old = Sigma
    Precision_old = Precision
    m_step_bound_old = m_step_bound

    beta_e = beta
    Sigma_e = Sigma

    if (!is.null(true_beta))
      beta_e = true_beta
    if (!is.null(true_Sigma))
      Sigma_e = true_Sigma

    E_quantities = mnp_probit_accumulate_cpp(
      X = X,
      Y = Y_int,
      obs_set = obs_set,
      beta = beta_e,
      Sigma = Sigma_e,
      Precision = Precision,
      E_method = E_method,
      n_mc = as.integer(n_mc)
    )

    # M step ----
    # beta estimation
    beta_new = solve(E_quantities$gls_a, E_quantities$gls_b)

    # E[Q | -]
    E_sample_cov = E_quantities$E_second_moment
    mu_X_beta = matrix(0, nrow = m, ncol = m)
    X_beta_cross = matrix(0, nrow = m, ncol = m)

    for (k in 1:p)
      mu_X_beta = mu_X_beta + beta_new[k, 1] * E_quantities$mu_X_sum[, , k]

    for (k in 1:p)
    {
      for (l in 1:p)
        X_beta_cross =
          X_beta_cross + beta_new[k, 1] * beta_new[l, 1] * E_quantities$X_cross_sum[, , k, l]
    }

    E_sample_cov = E_sample_cov - mu_X_beta - t(mu_X_beta) + X_beta_cross
    E_sample_cov = E_sample_cov / E_quantities$n_used

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

    # m_step_bound_new = - (sum(diag(E_sample_cov %*% Precision_new)) + determinant(Precision_new)$modulus )

    # update parameters
    beta = beta_new * (1 - M_damping) + M_damping * beta_old
    Precision = Precision_new * (1 - M_damping) + M_damping * Precision_old
    Sigma = Sigma_new * (1 - M_damping) + M_damping * Sigma_old
    m_step_bound = - (sum(diag(E_sample_cov %*% Precision)) + determinant(Precision)$modulus )

    #### record history ####
    if (record_history)
    {
      llik[iter] = m_step_bound
      beta_history[iter+1, ] = beta
      Sigma_history[iter+1, , ] = Sigma
      Prec_history[iter+1, , ] = Precision
    }

    #### update convergence ####
    if (conv_metric == "precision")
    {
      param_new = c(Precision, beta)
      param_old = c(Precision_old, beta_old)
      # dmetric = max(abs(param_new - param_old) / abs(param_old))
      dmetric = max(abs(param_new - param_old))
    }
    else if (conv_metric == "precision_relative")
    {
      param_new = c(Precision, beta)
      param_old = c(Precision_old, beta_old)
      # dmetric = max(abs(param_new - param_old) / abs(param_old))
      dmetric = max(abs(param_new - param_old) / abs(param_old))
    }
    else if (conv_metric == "covariance")
    {
      param_new = c(Sigma, beta)
      param_old = c(Sigma_old, beta_old)
      # dmetric = max(abs(param_new  - param_old) / abs(param_old))
      dmetric = max(abs(param_new - param_old))
    }
    else if (conv_metric == "mbound")
    {
      dmetric = abs(m_step_bound - m_step_bound_old)
    }
    else if (conv_metric == "beta")
    {
      dmetric = max(abs(beta - beta_old))
    }
    else if (conv_metric == "beta_relative")
    {
      dmetric = max(abs(beta - beta_old) / abs(beta_old))
    }
    else
      stop("invalid convergence metric")

    if (verbose != 0 && iter %% verbose == 0)
    {
      # print(paste0("EM iteration: ", iter))
      print(paste0("EM iteration: ", iter, ", M step bound: ", m_step_bound))
      print(paste0("coefficient estimate: ", beta))
    }

    #### update parameters ####
    iter = iter + 1
    # beta = beta_new
    # Sigma = Sigma_new
    # Precision = Precision_new

    if (!isSymmetric.matrix(Sigma))
    {
      print("Non symmetric covariance, stopping")
      break
    }
  }
  elapsed = tictoc::toc(quiet = (!(verbose != 0)))

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
    time = elapsed$toc - elapsed$tic,
    conv_metric = conv_metric,
    tol = tol
  ))
}
