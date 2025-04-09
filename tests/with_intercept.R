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


simdata = generate_identified_choice_data(
  n_obs = n_obs,
  x_range = c(-.5,.5),
  coef_true = coef_true,
  Sigma_iden = Sigma_init,
  seed = 1
)
