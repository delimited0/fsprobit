library(fsprobit)
library(mvtnorm)
library(doParallel)

p = 1
n_obs = 10000
n_choices = 6

tol = .001
conv_metric = "covariance"
max_iter = 500

random_covariance_matrix = function(size) {
  Z = matrix(rnorm(size * size), nrow = size)
  crossprod(Z) + diag(size) * 1e-6
}

# true parameters ----
# Prec_iden = diag(1, n_choices-1)
# Prec_iden[abs(row(Prec_iden) - col(Prec_iden)) == 1] = .5
# Sigma_iden = solve(Prec_iden)
Sigma_iden = diag(1, n_choices-1)
Sigma_iden[abs(row(Sigma_iden) - col(Sigma_iden)) == 1] = -.5
Prec_iden = solve(Sigma_iden)

coef_true = as.matrix(c(2))

# initial parameters
Sigma_init = random_covariance_matrix(n_choices-1)
coef_init = as.matrix(c(.2))

# covariate mean and sd
n_mean = 0
n_sd = 1

simdata = generate_custom_identified_choice_data(
  n_obs = n_obs,
  sampler = function(n) rnorm(n, mean=n_mean, sd=n_sd),
  coef_true = coef_true,
  Sigma_iden = Sigma_iden,
  seed = 1
)


# compare EP and HMC ----
ep_fit = mnp_probit_accelerated(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = Sigma_init,
  E_method = "EPMNP",
  M_method = "Covariance",
  n_choices = n_choices,
  true_trace = n_choices-1,
  tol = tol,
  max_iter = max_iter,
  conv_metric = conv_metric,
  shift_iden_method = "ref",
  scale_iden_method = "topleft",
  record_history = TRUE,
  verbose = 5
)

ep_fit$llik |> plot(type="l")

ep_fit$beta_history |> plot(type = "l")

sigma_history = ep_fit$Sigma_history
off_diag_idx = which(upper.tri(sigma_history[1, , ]), arr.ind = TRUE)

# Classify elements using the true banded covariance, rather than whether an
# estimate happens to be non-zero during optimization.
true_off_diag = Sigma_iden[off_diag_idx]
nonzero_off_diag_idx = off_diag_idx[true_off_diag != 0, , drop = FALSE]
zero_off_diag_idx = off_diag_idx[true_off_diag == 0, , drop = FALSE]

extract_covariance_history = function(index) {
  history = sapply(
    seq_len(nrow(index)),
    function(i) sigma_history[, index[i, 1], index[i, 2]]
  )
  if (is.null(dim(history)))
    history = matrix(history, ncol = 1)
  history
}

nonzero_off_diag_history = extract_covariance_history(nonzero_off_diag_idx)
zero_off_diag_history = extract_covariance_history(zero_off_diag_idx)
diagonal_history = sapply(
  seq_len(nrow(Sigma_iden)),
  function(i) sigma_history[, i, i]
)

old_mfrow = par(mfrow = c(3, 1))

diagonal_cols = seq_len(ncol(diagonal_history))
matplot(
  diagonal_history,
  type = "l",
  lty = 1,
  col = diagonal_cols,
  xlab = "Iteration",
  ylab = "Estimated covariance",
  main = "Diagonal covariance elements"
)
abline(h = 1, lty = 2, lwd = 2)
legend(
  "topright",
  legend = paste0("(", seq_len(ncol(diagonal_history)), ", ",
                  seq_len(ncol(diagonal_history)), ")"),
  col = diagonal_cols,
  lty = 1,
  cex = 0.8
)

nonzero_cols = seq_len(ncol(nonzero_off_diag_history))
matplot(
  nonzero_off_diag_history,
  type = "l",
  lty = 1,
  col = nonzero_cols,
  xlab = "Iteration",
  ylab = "Estimated covariance",
  main = "True nonzero off-diagonal covariance elements"
)
abline(h = unique(Sigma_iden[nonzero_off_diag_idx]), lty = 2, lwd = 2)
legend(
  "topright",
  legend = paste0(
    "(", nonzero_off_diag_idx[, 1], ", ",
    nonzero_off_diag_idx[, 2], ")"
  ),
  col = nonzero_cols,
  lty = 1,
  cex = 0.8
)

zero_cols = seq_len(ncol(zero_off_diag_history))
matplot(
  zero_off_diag_history,
  type = "l",
  lty = 1,
  col = zero_cols,
  xlab = "Iteration",
  ylab = "Estimated covariance",
  main = "True zero off-diagonal covariance elements"
)
abline(h = 0, lty = 2, lwd = 2)
legend(
  "topright",
  legend = paste0(
    "(", zero_off_diag_idx[, 1], ", ", zero_off_diag_idx[, 2], ")"
  ),
  col = zero_cols,
  lty = 1,
  cex = 0.8
)

par(old_mfrow)


hmc_fit = mnp_probit_accelerated(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = Sigma_init,
  E_method = "HMC",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  max_iter = max_iter,
  conv_metric = conv_metric,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  record_history = TRUE,
  verbose = 5
)

hmc_fit$beta_history |> plot(type = "l")
hmc_fit$Sigma_history[,1,1] |> plot(type = "l")
hmc_fit$Sigma_history[,1,2] |> plot(type = "l")

met_fit = mnp_probit_accelerated(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = Sigma_init,
  E_method = "MET",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  max_iter = max_iter,
  conv_metric = conv_metric,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  record_history = TRUE,
  verbose = 5
)

met_fit$beta_history |> plot(type = "l")
met_fit$Sigma_history[,1,1] |> plot(type = "l")
met_fit$Sigma_history[,1,2] |> plot(type = "l")

# fit model ----
registerDoParallel(4)
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
  n_mc = 200,
  max_iter = max_iter,
  conv_metric = conv_metric,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  verbose=5
)

# check result ---
probit_trace_iden$Precision

probit_trace_iden$Sigma

# another example ----
Prec_iden = diag(1, n_choices-1)
Prec_iden[abs(row(Prec_iden) - col(Prec_iden)) == 1] = -.1
Sigma_iden = solve(Prec_iden)

coef_true = as.matrix(c(2))

# initial parameters
Sigma_init = random_covariance_matrix(n_choices-1)
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
