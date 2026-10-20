# Reproducible timing comparison for the legacy contraction tensors and the
# optimized stored-means second pass. Run against an installed development
# build with: Rscript tests/timing_accumulator_optimization.R

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

run_legacy = function(X, Y, beta, Sigma, Precision)
{
  quantities = mnp_probit_accumulate_cpp(
    X, Y, seq_len(length(Y)), beta, Sigma, Precision,
    "EPMNP", 25L, TRUE
  )
  beta_new = solve(quantities$gls_a, quantities$gls_b)
  list(
    beta = beta_new,
    covariance = legacy_covariance(X, beta_new, quantities),
    result_bytes = as.numeric(object.size(quantities))
  )
}

run_optimized = function(X, Y, beta, Sigma, Precision)
{
  obs_set = seq_len(length(Y))
  quantities = mnp_probit_accumulate_cpp(
    X, Y, obs_set, beta, Sigma, Precision,
    "EPMNP", 25L, FALSE
  )
  beta_new = solve(quantities$gls_a, quantities$gls_b)
  list(
    beta = beta_new,
    covariance = mnp_probit_residual_cov_cpp(
      X, obs_set, quantities$moment_mu, beta_new,
      quantities$E_second_moment
    ),
    result_bytes = as.numeric(object.size(quantities))
  )
}

elapsed = function(expression)
{
  gc()
  unname(system.time(force(expression))[["elapsed"]])
}

benchmark_case = function(label, n, m, p, repetitions)
{
  set.seed(100 + n + m + p)
  X = array(rnorm(n * m * p), dim = c(n, m, p))
  Y = sample.int(m + 1L, n, replace = TRUE)
  beta = matrix(rnorm(p), p, 1)
  sigma_seed = matrix(rnorm(m * m), m, m)
  Sigma = crossprod(sigma_seed) / m + diag(m)
  Precision = solve(Sigma)

  invisible(run_legacy(X, Y, beta, Sigma, Precision))
  invisible(run_optimized(X, Y, beta, Sigma, Precision))

  legacy_result = run_legacy(X, Y, beta, Sigma, Precision)
  optimized_result = run_optimized(X, Y, beta, Sigma, Precision)
  max_beta_difference = max(abs(legacy_result$beta - optimized_result$beta))
  max_covariance_difference = max(abs(
    legacy_result$covariance - optimized_result$covariance
  ))

  legacy_times = replicate(
    repetitions, elapsed(run_legacy(X, Y, beta, Sigma, Precision))
  )
  optimized_times = replicate(
    repetitions, elapsed(run_optimized(X, Y, beta, Sigma, Precision))
  )

  data.frame(
    case = label,
    n = n,
    m = m,
    p = p,
    repetitions = repetitions,
    legacy_seconds = median(legacy_times),
    optimized_seconds = median(optimized_times),
    speedup = median(legacy_times) / median(optimized_times),
    legacy_result_mb = legacy_result$result_bytes / 1024^2,
    optimized_result_mb = optimized_result$result_bytes / 1024^2,
    max_beta_difference = max_beta_difference,
    max_covariance_difference = max_covariance_difference
  )
}

repetitions = as.integer(Sys.getenv("FSPROBIT_BENCH_REPS", "5"))
if (is.na(repetitions) || repetitions < 1L) {
  stop("FSPROBIT_BENCH_REPS must be a positive integer")
}

cat("R:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("fsprobit:", as.character(packageVersion("fsprobit")), "\n")
commit = tryCatch(
  system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE),
  error = function(error) "unavailable"
)
cat("Commit:", commit[[1]], "\n")
cat("Repetitions:", repetitions, "\n\n")

results = rbind(
  benchmark_case("small_common", 10000L, 9L, 2L, repetitions),
  benchmark_case("large_m", 2000L, 49L, 5L, repetitions),
  benchmark_case("large_p", 2000L, 19L, 20L, repetitions)
)

print(results, row.names = FALSE, digits = 4)
