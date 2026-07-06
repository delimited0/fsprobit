# Timing comparison for dense EP vs sparse MNP EP
#
# This is an interactive script. Run it line by line to inspect the simulated
# data, moment agreement, and timing results.

library(fsprobit)

cat("\n=== Setup ===\n")

set.seed(202407)

n_obs = 150
n_choices = 8
m = n_choices - 1
p = 2

cat("Observations:", n_obs, "\n")
cat("Choices:", n_choices, "\n")
cat("Identified latent utility dimension:", m, "\n")
cat("Covariates per utility:", p, "\n")

# Build a moderate positive definite covariance for identified utilities.
Sigma_seed = matrix(rnorm(m * m), nrow = m)
Sigma_iden = crossprod(Sigma_seed) / m + diag(m)
Sigma_iden = Sigma_iden * (m / sum(diag(Sigma_iden)))

coef_true = matrix(c(0.8, -0.4), nrow = p)

cat("\nTrue coefficients:\n")
print(coef_true)

cat("\nIdentified covariance matrix:\n")
print(round(Sigma_iden, 3))

cat("\n=== Simulate MNP Data ===\n")

simdata = generate_identified_choice_data(
  n_obs = n_obs,
  x_range = c(-1, 1),
  coef_true = coef_true,
  Sigma_iden = Sigma_iden,
  seed = 11
)

X = simdata$X
Y = simdata$Y_cat

cat("X dimensions:\n")
print(dim(X))

cat("\nChoice counts:\n")
print(table(Y))

cat("\nFirst five observed choices:\n")
print(head(Y, 5))

cat("\n=== Sanity Check One Observation ===\n")

constraints = utility_shift_constraints(n_choices)

compute_sparse_ep_moment = function(Xbeta, Sigma, y)
{
  choice_index = as.integer(y) - 1
  epmnp(as.vector(Xbeta), Sigma, choice_index)
}

obs_idx = 1
Xbeta = X[obs_idx, , ] %*% coef_true
y = Y[obs_idx]
A = constraints[[y]]

cat("Observation:", obs_idx, "\n")
cat("Choice:", as.character(y), "\n")
cat("Sparse EP choice_index:", as.integer(y) - 1, "\n")
cat("Xbeta:\n")
print(round(Xbeta, 4))

old_one = mnp_ep_moments(Xbeta, Sigma_iden, y, A)
new_one = compute_sparse_ep_moment(Xbeta, Sigma_iden, y)

cat("\nMax absolute mean difference, old vs sparse:\n")
print(max(abs(as.vector(old_one$mu) - as.vector(new_one$mu))))

cat("\nMax absolute covariance difference, old vs sparse:\n")
print(max(abs(old_one$Sigma - new_one$Sigma)))

cat("\n=== E-Step Moment Timing ===\n")

compute_old_ep_moments = function()
{
  lapply(seq_len(n_obs), function(i) {
    y_i = Y[i]
    Xbeta_i = X[i, , ] %*% coef_true
    A_i = constraints[[y_i]]
    mnp_ep_moments(Xbeta_i, Sigma_iden, y_i, A_i)
  })
}

compute_sparse_ep_moments = function()
{
  lapply(seq_len(n_obs), function(i) {
    y_i = Y[i]
    Xbeta_i = X[i, , ] %*% coef_true
    compute_sparse_ep_moment(Xbeta_i, Sigma_iden, y_i)
  })
}

# Warm up both paths so one-time loading/compilation effects are not counted.
invisible(compute_old_ep_moments())
invisible(compute_sparse_ep_moments())

cat("Timing old dense EP path...\n")
old_time = system.time(old_moments <- compute_old_ep_moments())
print(old_time)

cat("\nTiming new sparse MNP EP path...\n")
sparse_time = system.time(sparse_moments <- compute_sparse_ep_moments())
print(sparse_time)

elapsed = c(
  old_dense_ep = unname(old_time[["elapsed"]]),
  sparse_mnp_ep = unname(sparse_time[["elapsed"]])
)

cat("\nElapsed seconds:\n")
print(elapsed)

cat("\nSpeedup, old elapsed / sparse elapsed:\n")
print(elapsed[["old_dense_ep"]] / elapsed[["sparse_mnp_ep"]])

max_mu_diff = max(vapply(seq_len(n_obs), function(i) {
  max(abs(as.vector(old_moments[[i]]$mu) - as.vector(sparse_moments[[i]]$mu)))
}, numeric(1)))

max_sigma_diff = max(vapply(seq_len(n_obs), function(i) {
  max(abs(old_moments[[i]]$Sigma - sparse_moments[[i]]$Sigma))
}, numeric(1)))

cat("\nAgreement over all observations:\n")
cat("Max mean difference:", max_mu_diff, "\n")
cat("Max covariance difference:", max_sigma_diff, "\n")

cat("\n=== Optional Full mnp_probit Timing ===\n")

# Set this to TRUE and run the block below when you want to compare the full
# model-fitting call. This includes both E-step moments and the M-step, so it is
# less isolated than the timing above.
RUN_FULL_MODEL_TIMING = TRUE

if (RUN_FULL_MODEL_TIMING) {
  common_args = list(
    X = X,
    Y = Y,
    true_trace = m,
    Sigma_init = Sigma_iden,
    beta_init = coef_true,
    max_iter = 1,
    update_beta = TRUE,
    record_history = FALSE,
    verbose = 0
  )

  cat("Timing mnp_probit with E_method = 'EP'...\n")
  old_model_time = system.time({
    old_fit = do.call(mnp_probit, c(common_args, list(E_method = "EP")))
  })
  print(old_model_time)

  cat("\nTiming mnp_probit with E_method = 'EPMNP'...\n")
  sparse_model_time = system.time({
    sparse_fit = do.call(mnp_probit, c(common_args, list(E_method = "EPMNP")))
  })
  print(sparse_model_time)

  model_elapsed = c(
    old_dense_ep = unname(old_model_time[["elapsed"]]),
    sparse_mnp_ep = unname(sparse_model_time[["elapsed"]])
  )

  cat("\nFull model elapsed seconds:\n")
  print(model_elapsed)

  cat("\nFull model speedup, old elapsed / sparse elapsed:\n")
  print(model_elapsed[["old_dense_ep"]] / model_elapsed[["sparse_mnp_ep"]])
}
