library(fsprobit)

legacy_covariance = function(X, beta_new, quantities)
{
  m = dim(X)[2]
  p = dim(X)[3]
  mu_X_beta = matrix(0, m, m)
  X_beta_cross = matrix(0, m, m)

  for (k in seq_len(p)) {
    mu_X_beta = mu_X_beta + beta_new[k, 1] * quantities$mu_X_sum[, , k]
    for (l in seq_len(p)) {
      X_beta_cross = X_beta_cross +
        beta_new[k, 1] * beta_new[l, 1] *
        quantities$X_cross_sum[, , k, l]
    }
  }

  (quantities$E_second_moment - mu_X_beta - t(mu_X_beta) +
     X_beta_cross) / quantities$n_used
}

check_accumulators = function(
    n, m, p, method = "EPMNP", n_mc = 25L, subset = FALSE,
    tolerance = 1e-10, seed = 1L)
{
  set.seed(seed)
  X = array(rnorm(n * m * p), dim = c(n, m, p))
  Y = sample.int(m + 1L, n, replace = TRUE)
  Y[seq_len(m + 1L)] = seq_len(m + 1L)
  beta = matrix(rnorm(p), ncol = 1)
  sigma_seed = matrix(rnorm(m * m), m, m)
  Sigma = crossprod(sigma_seed) / m + diag(m)
  Precision = solve(Sigma)
  obs_set = if (subset) sample.int(n, max(1L, floor(0.7 * n))) else seq_len(n)

  set.seed(seed + 1000L)
  legacy = mnp_probit_accumulate_cpp(
    X, Y, obs_set, beta, Sigma, Precision, method, as.integer(n_mc), TRUE
  )
  legacy_rng = .Random.seed

  set.seed(seed + 1000L)
  optimized = mnp_probit_accumulate_cpp(
    X, Y, obs_set, beta, Sigma, Precision, method, as.integer(n_mc), FALSE
  )
  optimized_rng = .Random.seed

  beta_legacy = solve(legacy$gls_a, legacy$gls_b)
  beta_optimized = solve(optimized$gls_a, optimized$gls_b)
  covariance_legacy = legacy_covariance(X, beta_legacy, legacy)
  covariance_optimized = mnp_probit_residual_cov_cpp(
    X, obs_set, optimized$moment_mu, beta_optimized,
    optimized$E_second_moment
  )

  stopifnot(identical(legacy_rng, optimized_rng))
  stopifnot(length(legacy$moment_mu) == 0L)
  stopifnot(identical(dim(legacy$mu_X_sum), c(m, m, p)))
  stopifnot(identical(dim(legacy$X_cross_sum), c(m, m, p, p)))
  stopifnot(length(optimized$mu_X_sum) == 0L)
  stopifnot(length(optimized$X_cross_sum) == 0L)
  stopifnot(identical(dim(optimized$moment_mu), c(m, length(obs_set))))
  stopifnot(isTRUE(all.equal(legacy$gls_a, optimized$gls_a, tolerance = 0)))
  stopifnot(isTRUE(all.equal(legacy$gls_b, optimized$gls_b, tolerance = 0)))
  stopifnot(isTRUE(all.equal(
    legacy$E_second_moment, optimized$E_second_moment, tolerance = 0
  )))
  stopifnot(isTRUE(all.equal(
    beta_legacy, beta_optimized, tolerance = tolerance
  )))
  stopifnot(isTRUE(all.equal(
    covariance_legacy, covariance_optimized, tolerance = tolerance
  )))
  stopifnot(max(abs(covariance_optimized - t(covariance_optimized))) < tolerance)
}

# Deterministic EP cases cover base/non-base choices through balanced response
# draws, several matrix shapes, and both complete and sampled E-steps.
check_accumulators(30L, 2L, 1L, seed = 11L)
check_accumulators(80L, 5L, 3L, subset = TRUE, seed = 12L)
check_accumulators(60L, 9L, 6L, seed = 13L)

# Stochastic paths must consume exactly the same random stream because their
# moment calculations are unchanged by sufficient-statistic representation.
check_accumulators(
  12L, 2L, 2L, method = "HMC", n_mc = 20L,
  tolerance = 1e-9, seed = 21L
)
check_accumulators(
  10L, 2L, 2L, method = "MET", n_mc = 20L,
  tolerance = 1e-9, seed = 22L
)

# The legacy representation remains the default for direct callers.
X = array(rnorm(12), dim = c(3, 2, 2))
default_result = mnp_probit_accumulate_cpp(
  X, c(1L, 2L, 3L), 1:3, matrix(0, 2, 1), diag(2), diag(2),
  "EPMNP", 25L
)
stopifnot(identical(dim(default_result$X_cross_sum), c(2L, 2L, 2L, 2L)))

# Validation failures in the new path.
stopifnot(inherits(try(
  mnp_probit_accumulate_cpp(
    X, c(1L, 2L, 3L), integer(), matrix(0, 2, 1), diag(2), diag(2),
    "EPMNP", 25L, FALSE
  ), silent = TRUE
), "try-error"))
stopifnot(inherits(try(
  mnp_probit_residual_cov_cpp(
    X, 1:3, matrix(0, 3, 3), matrix(0, 2, 1), diag(2)
  ), silent = TRUE
), "try-error"))
stopifnot(inherits(try(
  mnp_probit_residual_cov_cpp(
    X, c(1L, 4L), matrix(0, 2, 2), matrix(0, 2, 1), diag(2)
  ), silent = TRUE
), "try-error"))

# Exercise the production integration for more than one EM iteration.
set.seed(31)
fit_X = array(rnorm(120 * 4 * 3), dim = c(120, 4, 3))
fit_Y = factor(sample.int(5, 120, replace = TRUE), levels = 1:5)
fit = mnp_probit_accelerated(
  fit_X, fit_Y, true_trace = 4, beta_init = matrix(0, 3, 1),
  Sigma_init = diag(4), E_method = "EPMNP", max_iter = 3,
  tol = 0, record_history = TRUE
)
stopifnot(all(is.finite(fit$beta)))
stopifnot(all(is.finite(fit$Sigma)))
stopifnot(all(is.finite(fit$Precision)))
stopifnot(fit$iters == 4)
