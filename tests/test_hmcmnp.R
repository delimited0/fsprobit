library(fsprobit)

exact_univariate_moments = function(mu, variance, lower)
{
  sd = sqrt(variance)
  alpha = -mu / sd
  if (lower) {
    lambda = dnorm(alpha) / pnorm(alpha, lower.tail = FALSE)
    mean = mu + sd * lambda
    var = variance * (1 + alpha * lambda - lambda^2)
  } else {
    lambda = dnorm(alpha) / pnorm(alpha)
    mean = mu - sd * lambda
    var = variance * (1 - alpha * lambda - lambda^2)
  }
  c(mean = mean, var = var)
}

for (choice_index in 0:1) {
  set.seed(100 + choice_index)
  result = hmcmnp(
    mu = 0.35,
    Precision = matrix(1 / 1.7),
    choice_index = choice_index,
    n_mc = 20000
  )
  expected = exact_univariate_moments(0.35, 1.7, choice_index == 1)
  stopifnot(abs(result$mu - expected[["mean"]]) < 0.05)
  stopifnot(abs(result$Sigma[1, 1] - expected[["var"]]) < 0.07)
}

mu = c(0.3, -0.2, 0.1)
Precision = matrix(
  c(1.4, 0.2, 0.1, 0.2, 1.1, -0.15, 0.1, -0.15, 0.9),
  nrow = 3
)
constraints = utility_shift_constraints(4)

for (choice_index in 0:3) {
  A = constraints[[choice_index + 1]]
  f = if (choice_index == 0) -A else A
  initial = rep(-1, 3)
  if (choice_index > 0)
    initial[choice_index] = 1

  set.seed(700 + choice_index)
  result = hmcmnp(mu, Precision, choice_index, 10000)
  set.seed(900 + choice_index)
  reference_samples = tmg::rtmg(
    n = 10000,
    M = Precision,
    r = as.vector(Precision %*% mu),
    initial = initial,
    f = f,
    g = rep(0, 3),
    burn.in = 30
  )

  stopifnot(max(abs(result$mu - colMeans(reference_samples))) < 0.08)
  stopifnot(max(abs(result$Sigma - cov(reference_samples))) < 0.1)
  stopifnot(max(abs(result$Sigma - t(result$Sigma))) < 1e-12)
  stopifnot(min(eigen(result$Sigma, symmetric = TRUE)$values) > -1e-10)
}

set.seed(42)
first = hmcmnp(c(0.2, -0.1), matrix(c(1.3, 0.2, 0.2, 0.9), 2), 2, 100)
set.seed(42)
second = hmcmnp(c(0.2, -0.1), matrix(c(1.3, 0.2, 0.2, 0.9), 2), 2, 100)
stopifnot(identical(first, second))

stopifnot(inherits(try(hmcmnp(c(0, 0), diag(3), 0, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(hmcmnp(c(0, 0), matrix(c(1, 1, 0, 1), 2), 0, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(hmcmnp(c(0, 0), diag(c(1, -1)), 0, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(hmcmnp(c(0, NA), diag(2), 0, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(hmcmnp(c(0, 0), diag(2), -1, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(hmcmnp(c(0, 0), diag(2), 3, 10), silent = TRUE), "try-error"))
stopifnot(inherits(try(hmcmnp(c(0, 0), diag(2), 0, 1), silent = TRUE), "try-error"))

set.seed(11)
simdata = generate_identified_choice_data(
  n_obs = 30,
  x_range = c(-1, 1),
  coef_true = matrix(c(0.4, -0.2), nrow = 2),
  Sigma_iden = diag(2),
  seed = 12
)
fit = mnp_probit_accelerated(
  X = simdata$X,
  Y = simdata$Y_cat,
  true_trace = 2,
  Sigma_init = diag(2),
  beta_init = matrix(0, nrow = 2),
  E_method = "HMC",
  n_mc = 40,
  max_iter = 2,
  tol = 0,
  record_history = TRUE
)
stopifnot(all(is.finite(fit$beta)))
stopifnot(all(is.finite(fit$Sigma)))
stopifnot(fit$E_method == "HMC")
stopifnot(fit$iters == 3)

stopifnot(inherits(
  try(mnp_probit_accelerated(
    X = simdata$X,
    Y = simdata$Y_cat,
    true_trace = 2,
    E_method = "HMC",
    shift_iden_method = "zero"
  ), silent = TRUE),
  "try-error"
))
