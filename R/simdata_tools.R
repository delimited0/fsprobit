# dataset simulation
library(mvtnorm)
library(CVXR)

#' generate latent utilties for mnp choice model
#' @param means n observations x m choices matrix of individual utility means
#' @param precision m choices x m choices precision matrix of utility
rchoicemvn = function(means, precision)
{
  m = ncol(means)  # number of choices
  n = nrow(means)  # number of observations
  Z = matrix(rnorm(m*n), nrow = m, ncol = n)
  U = chol(precision)
  UinvZ = backsolve(U, Z)
  t(UinvZ + t(means))
}

#' Fast vectorized multivariate normal generation with observation-specific means
#' @param means matrix of means (n_obs x n_vars)
#' @param sigma covariance matrix (n_vars x n_vars)
#' @param n_obs number of observations
#' @param n_vars number of variables
fast_mvrnorm = function(means, sigma, n_obs = NULL, n_vars = NULL) {
  if (is.null(n_obs)) n_obs = nrow(means)
  if (is.null(n_vars)) n_vars = ncol(means)
  
  # Generate standard normal random variables
  Z = matrix(rnorm(n_obs * n_vars), nrow = n_vars, ncol = n_obs)
  
  # Cholesky decomposition (computed once)
  L = chol(sigma)
  
  # Transform: L %*% Z + means^T
  # L %*% Z gives us the correlated random variables
  # Adding t(means) broadcasts the observation-specific means
  result = t(L) %*% Z + t(means)
  
  return(t(result))  # Return as n_obs x n_vars
}

# simulate latent utility with given covariance
# generate covariates from uniform distribution
#' @param n_obs number of observations
#' @param n_choices number of choices
#' @param x_range domain of uniform covariates
#' @param coef_true true coefficients
#' @param Sigma_true true covariance matrix
#' @param seed random seed
#' @returns list with covariates X, latent utilities Z, observed choices Y_cat, and identification transformed covariates X_iden
generate_choice_data = function(
    n_obs,
    n_choices,
    x_range,
    coef_true,
    Sigma_true,
    seed = 1
)
{
  p = nrow(coef_true)

  # data generation
  set.seed(seed)
  X = array(
    data = runif(n_obs * n_choices * p,
                 min = x_range[1], max = x_range[2]),
    dim = c(n_obs, n_choices, p),
    dimnames = list("obs" = paste0("obs_", 1:n_obs),
                    "choice" = paste0("choice_", 1:n_choices),
                    "covariate" = paste0("covariate_", 1:p))
  )
  
  # Vectorized mean calculation: compute X[i,,] %*% coef_true for all i
  # Reshape X to (n_obs * n_choices, p) for matrix multiplication
  X_reshaped = matrix(X, nrow = n_obs * n_choices, ncol = p)
  means_flat = X_reshaped %*% coef_true
  # Reshape back to (n_obs, n_choices)
  means = matrix(means_flat, nrow = n_obs, ncol = n_choices)
  
  # Use fast vectorized multivariate normal generation
  Z = fast_mvrnorm(means, Sigma_true)
  
  Y = apply(Z, 1, function(z) {
    return(which.max(z))
  })
  Y_cat = factor(Y, levels = 1:n_choices)

  # Vectorized identification transforms
  X_iden = array(NA, dim = c(n_obs, n_choices-1, p))
  iden_mat = cbind(-1, diag(n_choices-1))
  for (i in 1:n_obs) {
    X_iden[i, , ] = iden_mat %*% X[i, ,]
  }

  result = list(
    X = X,
    Z = Z,
    Y_cat = Y_cat,
    X_iden = X_iden
  )
}

# simulate relative utiltiies and choices given identified covariance and coefficients
# generate covariates from uniform distribution
#' @param n_obs number of observations
#' @param x_range domain of uniform covariates
#' @param coef_true true coefficients
#' @param Sigma_iden identified covariance matrix
#' @param seed random seed
#' @returns list with covariates X, relative utilities Z_rel, observed choices Y_cat, and identification transformed covariates X_iden
generate_identified_choice_data = function(
    n_obs,
    x_range,
    coef_true,
    Sigma_iden,
    seed = 1
)
{
  p = nrow(coef_true)

  n_choices = nrow(Sigma_iden)+1

  # simulate covariates from uniform in original dimension
  set.seed(seed)
  covariates = array(
    data = runif(n_obs * (n_choices-1) * p,
                 min = x_range[1], max = x_range[2]),
    dim = c(n_obs, n_choices-1, p),
    dimnames = list("obs" = paste0("obs_", 1:n_obs),
                    "choice" = paste0("choice_", 1:(n_choices-1)),
                    "covariate" = paste0("covariate_", 1:p))
  )

  # Vectorized mean calculation for relative utilities
  # Reshape covariates to (n_obs * (n_choices-1), p) for matrix multiplication
  covariates_reshaped = matrix(covariates, nrow = n_obs * (n_choices-1), ncol = p)
  means_flat = covariates_reshaped %*% coef_true
  # Reshape back to (n_obs, n_choices-1)
  means = matrix(means_flat, nrow = n_obs, ncol = n_choices-1)
  
  # Use fast vectorized multivariate normal generation
  relative_utilities = fast_mvrnorm(means, Sigma_iden)
  
  Y = apply(relative_utilities, 1, function(z) {
    if (all(z < 0))
    {
      return(1)
    }
    else
    {
      return(which.max(z)+1)
    }
  })
  Y_cat = factor(Y, levels = 1:n_choices)

  result = list(
    X = covariates,
    Z_rel = relative_utilities,
    Y_cat = Y_cat
  )
}

#' generate covariates from custom probability distribution. Provide identified
#' covariance as argument.
#' @param n_obs number of observations
#' @param sampler function to generate samples from some distribution
#' @export
generate_custom_identified_choice_data = function(
    n_obs,
    sampler,
    coef_true,
    Sigma_iden,
    seed = 1
)
{
  p = nrow(coef_true)

  n_choices = nrow(Sigma_iden)+1

  # simulate covariates from custom sampler
  set.seed(seed)
  covariates = array(
    data = sampler(n = n_obs * (n_choices-1) * p),
    dim = c(n_obs, n_choices-1, p),
    dimnames = list("obs" = paste0("obs_", 1:n_obs),
                    "choice" = paste0("choice_", 1:(n_choices-1)),
                    "covariate" = paste0("covariate_", 1:p))
  )

  # Vectorized mean calculation for relative utilities
  # Reshape covariates to (n_obs * (n_choices-1), p) for matrix multiplication
  covariates_reshaped = matrix(covariates, nrow = n_obs * (n_choices-1), ncol = p)
  means_flat = covariates_reshaped %*% coef_true
  # Reshape back to (n_obs, n_choices-1)
  means = matrix(means_flat, nrow = n_obs, ncol = n_choices-1)
  
  # Use fast vectorized multivariate normal generation
  relative_utilities = fast_mvrnorm(means, Sigma_iden)
  
  Y = apply(relative_utilities, 1, function(z) {
    if (all(z < 0))
    {
      return(1)
    }
    else
    {
      return(which.max(z)+1)
    }
  })
  Y_cat = factor(Y, levels = 1:n_choices)

  result = list(
    X = covariates,
    Z_rel = relative_utilities,
    Y_cat = Y_cat
  )
}

# given identified precision, generate data according to corresponding not identified precision
generate_target_identified_choice_data = function(
    n_obs,
    x_range,
    coef_true,
    Prec_iden,
    seed = 1
)
{
  p = nrow(coef_true)
  n_choices = nrow(Prec_iden)+1  # original number of choices

  # given the identified precision, recover the not identified precision
  # we assume the reference choice is the first one
  iden_mat = cbind(-1, diag(n_choices-1))
  m = n_choices
  Prec = Variable(m, m, PSD=TRUE)

  objective = 1
  constr = list(
    iden_mat %*% Prec %*% t(iden_mat) == Prec_iden
  )

  prob = Problem(Maximize(objective), constr)
  cvx_solution = CVXR::psolve(prob, abstol = 1e-10, num_iter=100)
  Prec_not_iden = cvx_solution$getValue(Prec)

  # simulate covariates from uniform in original dimension
  set.seed(seed)
  covariates = array(
    data = runif(n_obs * (n_choices) * p,
                 min = x_range[1], max = x_range[2]),
    dim = c(n_obs, n_choices, p),
    dimnames = list("obs" = paste0("obs_", 1:n_obs),
                    "choice" = paste0("choice_", 1:(n_choices)),
                    "covariate" = paste0("covariate_", 1:p))
  )

  # Vectorized mean calculation: compute covariates[i,,] %*% coef_true for all i
  # Reshape covariates to (n_obs * n_choices, p) for matrix multiplication
  covariates_reshaped = matrix(covariates, nrow = n_obs * n_choices, ncol = p)
  means_flat = covariates_reshaped %*% coef_true
  # Reshape back to (n_obs, n_choices)
  means = matrix(means_flat, nrow = n_obs, ncol = n_choices)

  # Use optimized rchoicemvn function for precision-based sampling
  utilities = rchoicemvn(means, Prec_not_iden)

  Y = apply(utilities, 1, function(z) {
    return(which.max(z))
  })
  Y_cat = factor(Y, levels = 1:n_choices)

  # identification transforms for covariates
  covariates_iden = array(NA, dim = c(n_obs, n_choices-1, p))
  for (i in 1:n_obs) {
    covariates_iden[i, , ] = iden_mat %*% covariates[i, ,]
  }


  result = list(
    X = covariates,
    Z = utilities,
    Y_cat = Y_cat,
    X_iden = covariates_iden,
    Prec_not_iden = Prec_not_iden
  )
}
