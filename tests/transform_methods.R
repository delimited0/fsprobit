library(fsprobit)
library(mvtnorm)
library(doParallel)
library(MomTrunc)

p = 1
n_obs = 2000
n_choices = 25

# tol = 0.0005
tol = 1e-2
conv_metric = "beta_relative"
max_iter = 300

# true parameters ----
# Prec_iden = .2 * diag(n_choices-1) + .8 * rep(1, n_choices-1) %*% t(rep(1, n_choices-1))
# Sigma_iden = solve(Prec_iden)
Sigma_iden <- matrix(0.5, nrow = n_choices - 1, ncol = n_choices - 1)
diag(Sigma_iden) <- 1
Prec_iden <- solve(Sigma_iden)

# coef_true = as.matrix(c(2))
coef_true = as.matrix(c(1,-1))

# initial parameters
Sigma_init = diag(n_choices-1)
coef_init = as.matrix(c(.2, .2))

# covariate mean and sd
n_mean = 0
n_sd = 1

simdata = generate_custom_identified_choice_data(
  n_obs = n_obs,
  sampler = function(n) rnorm(n, mean=n_mean, sd=n_sd),
  coef_true = coef_true,
  Sigma_iden = solve(Prec_iden),
  seed = 1
)

constraints = utility_shift_constraints(n_choices)

# fit model ----
cl <- makeCluster(8) # Example: use all but one core
registerDoParallel(cl)
# registerDoSEQ()
probit_ep =  mnp_probit(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = Sigma_init,
  E_method = "EP",
  E_sample_rate = 1,
  M_method = "Newton",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  newton_tol = 1e-3,
  max_newton_iter = 50,
  max_iter = max_iter,
  conv_metric = conv_metric,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  verbose=5,
  record_history=TRUE,
  transform=TRUE
)
stopCluster(cl)

cl <- makeCluster(8) # Example: use all but one core
registerDoParallel(cl)
probit_ep_untransform = mnp_probit(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = Sigma_init,
  E_method = "EP",
  E_sample_rate = 1,
  M_method = "Newton",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  newton_tol = 1e-3,
  max_newton_iter = 50,
  max_iter = max_iter,
  conv_metric = "beta",
  n_mc = 1000,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  verbose=5,
  record_history=TRUE,
  transform=FALSE
)
stopCluster(cl)

cl <- makeCluster(8) # Example: use all but one core
registerDoParallel(cl)
probit_hmc = mnp_probit(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = Sigma_init,
  E_method = "HMC",
  E_sample_rate = 1,
  M_method = "Newton",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  newton_tol = 1e-3,
  max_newton_iter = 50,
  max_iter = max_iter,
  conv_metric = conv_metric,
  n_mc = 1000,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  verbose=5,
  record_history=TRUE,
  transform=FALSE
)
stopCluster(cl)




# registerDoParallel(8)
# registerDoSEQ()
probit_momtrunc = mnp_probit(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = Sigma_init,
  E_method = "MomTrunc",
  E_sample_rate = 1,
  M_method = "Newton",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  newton_tol = 1e-3,
  max_newton_iter = 50,
  max_iter = max_iter,
  conv_metric = conv_metric,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  verbose=5,
  record_history=TRUE
)

probit_met = mnp_probit(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = Sigma_init,
  E_method = "MET",
  E_sample_rate = 1,
  M_method = "Newton",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  newton_tol = 1e-3,
  max_newton_iter = 50,
  max_iter = max_iter,
  conv_metric = conv_metric,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  verbose=5,
  record_history=TRUE
)


Xbeta = matrix(c(0.1, 0.2, 0.3))
# Sigma = diag(n_choices-1)
Sigma = .3 * diag(n_choices-1) + .7 * rep(1, n_choices-1) %*% t(rep(1, n_choices-1))
y = factor(1)
A = constraints[[y]]

m = nrow(A)
lb = rep(0, m)
ub = rep(Inf, m)

MomTrunc::meanvarTMD(
  lower = lb,
  upper = ub,
  mu = Xbeta,
  lambda=0,
  tau=0,
  Sigma=Sigma,
  dist='normal'
)

epmgpr::moments(
  lb = lb,
  ub = ub,
  mu = Xbeta,
  Sigma = Sigma
)

samples = TruncatedNormal::rtmvnorm(
  n=1000,
  mu=as.vector(Xbeta),
  sigma=Sigma,
  lb=lb,
  ub=ub
)
colMeans(samples)
cov(samples)

# transform back
# recycle Xbeta to add mean per sample
samples = t( A %*% t(samples) + as.vector(Xbeta) )





mnp_momtrunc_moments(Xbeta, Sigma, y, A)

moments = mnp_met_moments(Xbeta, Sigma, y, A, 2000)
ep_moments = mnp_ep_moments(Xbeta, Sigma, y, A)

n_mc =100
m = nrow(Xbeta)
samples = TruncatedNormal::rtmvnorm(
  n=n_mc,
  mu=as.vector(Xbeta),
  sigma=Sigma,
  lb=rep(-Inf, m),
  ub=rep(0, m),
)


# block -------------------------------------------------------------------
library(Matrix)
p = 1
n_obs = 2000
n_choices = 10

tol = 1e-5
conv_metric = "precision"
relerr_tol = .1
max_iter = 500

set.seed(12)

# Function to generate a random positive definite matrix
generate_cov_matrix <- function(size) {
  M <- matrix(rnorm(size^2), nrow = size)
  cov_matrix <- crossprod(M)  # Make it symmetric positive definite
  cov_matrix <- cov_matrix / max(abs(cov_matrix))  # Normalize for stability
  return(cov_matrix)
}

generate_block_diag_cov <- function(n) {
  # Generate random block sizes summing to n
  sizes <- c()
  remaining <- n
  while (remaining > 0) {
    size <- sample(1:min(remaining, 5), 1)  # Limit block size to 5 for variability
    sizes <- c(sizes, size)
    remaining <- remaining - size
  }

  # Construct the block diagonal matrix
  blocks <- lapply(sizes, generate_cov_matrix)
  cov_matrix <- bdiag(blocks)

  return(as.matrix(cov_matrix))
}

Prec_iden = generate_block_diag_cov(n_choices-1)
Sigma_iden = solve(Prec_iden)

coef_true = as.matrix(c(2))

# initial parameters
Sigma_init = diag(n_choices-1)
coef_init = as.matrix(c(.2))

# covariate mean and sd
n_mean = 0
n_sd = 1

simdata = generate_custom_identified_choice_data(
  n_obs = n_obs,
  sampler = function(n) rnorm(n, mean=n_mean, sd=n_sd),
  coef_true = coef_true,
  Sigma_iden = solve(Prec_iden),
  seed = 1
)

probit_ep_iden =  mnp_probit(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = Sigma_init,
  E_method = "EP",
  E_sample_rate = 1,
  M_method = "Newton",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  newton_tol = 1e-3,
  max_newton_iter = 50,
  max_iter = max_iter,
  conv_metric = conv_metric,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  verbose=5,
  record_history=TRUE,
  transform=TRUE
)




