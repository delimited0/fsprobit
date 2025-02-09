library(fsprobit)
library(mvtnorm)
library(doParallel)

p = 1
n_obs = 2000
n_choices = 4

tol = .001
conv_metric = "precision"
relerr_tol = .1
max_iter = 500

# true parameters ----
Prec_iden = .2 * diag(n_choices-1) + .8 * rep(1, n_choices-1) %*% t(rep(1, n_choices-1))
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

constraints = utility_shift_constraints(n_choices)

# fit model ----
# registerDoParallel(4)
registerDoSEQ()
probit_trace_iden = mnp_probit(
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
  record_history=TRUE
)


Xbeta = matrix(c(0.1, 0.2, 0.3))
# Sigma = diag(n_choices-1)
Sigma = .3 * diag(n_choices-1) + .7 * rep(1, n_choices-1) %*% t(rep(1, n_choices-1))
y = factor(1)
A = constraints[[y]]

# mnp_momtrunc_moments(Xbeta, Sigma, y, A)

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



