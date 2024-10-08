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

# fit model ----
registerDoParallel(4)
probit_trace_iden = mnp_probit(
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
  verbose=5
)

# check result ---
probit_trace_iden$Precision
