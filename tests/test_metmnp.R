library(fsprobit)

# Univariate base-choice moments have a closed form.
set.seed(11)
mu = 0.4
variance = 1.7
sd = sqrt(variance)
alpha = -mu / sd
mills = dnorm(alpha) / pnorm(alpha)
analytic_mean = mu - sd * mills
analytic_variance = variance * (1 - alpha * mills - mills^2)
univariate = metmnp(mu, matrix(variance), 0, 30000)
stopifnot(abs(univariate$mu - analytic_mean) < 0.025)
stopifnot(abs(univariate$Sigma - analytic_variance) < 0.035)

# Compare every reference-identified choice with an independent implementation.
if (requireNamespace("TruncatedNormal", quietly = TRUE)) {
  mu = c(0.3, -0.2, 0.1)
  Sigma = matrix(c(
    1.2, 0.25, -0.1,
    0.25, 0.9, 0.15,
    -0.1, 0.15, 1.1
  ), 3, 3)

  for (choice_index in 0:3) {
    if (choice_index == 0) {
      A = -diag(3)
    } else {
      A = -diag(3)
      A[, choice_index] = 1
      A[choice_index, choice_index] = 1
      stopifnot(max(abs(A %*% A - diag(3))) < 1e-12)
    }

    set.seed(100 + choice_index)
    result = metmnp(mu, Sigma, choice_index, 20000)
    set.seed(200 + choice_index)
    transformed = TruncatedNormal::rtmvnorm(
      n = 20000,
      mu = as.vector(A %*% mu),
      sigma = A %*% Sigma %*% t(A),
      lb = rep(0, 3),
      ub = rep(Inf, 3)
    )
    reference = t(solve(A, t(transformed)))
    stopifnot(max(abs(result$mu - colMeans(reference))) < 0.05)
    stopifnot(max(abs(result$Sigma - cov(reference))) < 0.06)
  }
}

# R's RNG controls the complete sampler.
set.seed(77)
first = metmnp(c(0.2, -0.1), matrix(c(1.3, 0.2, 0.2, 0.9), 2), 2, 100)
set.seed(77)
second = metmnp(c(0.2, -0.1), matrix(c(1.3, 0.2, 0.2, 0.9), 2), 2, 100)
stopifnot(identical(first, second))

# Input validation.
stopifnot(inherits(try(metmnp(c(0, 0), diag(3), 0, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(metmnp(c(0, 0), matrix(c(1, 1, 0, 1), 2), 0, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(metmnp(c(0, 0), diag(c(1, -1)), 0, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(metmnp(c(0, NA), diag(2), 0, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(metmnp(c(0, 0), diag(2), -1, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(metmnp(c(0, 0), diag(2), 3, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(metmnp(c(0, 0), diag(2), 0, 1), silent = TRUE), "try-error"))

# Accelerated EM dispatches MET through the C++ sufficient-statistic accumulator.
simdata = generate_identified_choice_data(
  n_obs = 80,
  x_range = c(-1, 1),
  coef_true = matrix(0.4),
  Sigma_iden = diag(2),
  seed = 9
)
set.seed(10)
fit = mnp_probit_accelerated(
  X = simdata$X,
  Y = simdata$Y_cat,
  true_trace = 2,
  n_choices = 3,
  E_method = "MET",
  n_mc = 20,
  max_iter = 1,
  tol = 0,
  verbose = 0
)
stopifnot(fit$E_method == "MET")
stopifnot(all(is.finite(fit$beta)))
stopifnot(all(is.finite(fit$Sigma)))

stopifnot(inherits(try(
  mnp_probit_accelerated(
    X = simdata$X,
    Y = simdata$Y_cat,
    true_trace = 2,
    n_choices = 3,
    E_method = "MET",
    shift_iden_method = "sum",
    max_iter = 1
  ),
  silent = TRUE
), "try-error"))
