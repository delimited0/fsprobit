# Performance Test and Validation Script for Optimized simdata_tools.R
# This script compares original vs optimized functions and validates correctness

library(mvtnorm)
library(CVXR)
source("R/simdata_tools.R")

# Function to create backup of original functions for comparison
create_original_functions <- function() {
  # Original generate_choice_data function (slow version)
  generate_choice_data_original <<- function(
      n_obs,
      n_choices,
      x_range,
      coef_true,
      Sigma_true,
      seed = 1
  ) {
    p = nrow(coef_true)
    
    set.seed(seed)
    X = array(
      data = runif(n_obs * n_choices * p,
                   min = x_range[1], max = x_range[2]),
      dim = c(n_obs, n_choices, p),
      dimnames = list("obs" = paste0("obs_", 1:n_obs),
                      "choice" = paste0("choice_", 1:n_choices),
                      "covariate" = paste0("covariate_", 1:p))
    )
    Z = t(apply(X, 1,
                function(x) rmvnorm(1, as.matrix(x) %*% coef_true, Sigma_true))
    )
    Y = apply(Z, 1, function(z) {
      return(which.max(z))
    })
    Y_cat = factor(Y, levels = 1:n_choices)
    
    X_iden = array(NA, dim = c(n_obs, n_choices-1, p))
    iden_mat = cbind(-1, diag(n_choices-1))
    for (i in 1:n_obs) {
      X_iden[i, , ] = iden_mat %*% X[i, ,]
    }
    
    result = list(
      X = X,
      Z = Z,
      Y_cat = Y_cat,
      X_iden = X_iden
    )
  }
  
  # Original generate_identified_choice_data function (slow version)
  generate_identified_choice_data_original <<- function(
      n_obs,
      x_range,
      coef_true,
      Sigma_iden,
      seed = 1
  ) {
    p = nrow(coef_true)
    n_choices = nrow(Sigma_iden)+1
    
    set.seed(seed)
    covariates = array(
      data = runif(n_obs * (n_choices-1) * p,
                   min = x_range[1], max = x_range[2]),
      dim = c(n_obs, n_choices-1, p),
      dimnames = list("obs" = paste0("obs_", 1:n_obs),
                      "choice" = paste0("choice_", 1:(n_choices-1)),
                      "covariate" = paste0("covariate_", 1:p))
    )
    
    relative_utilities = t(apply(covariates, 1, function(x) rmvnorm(1, as.matrix(x) %*% coef_true, Sigma_iden)))
    Y = apply(relative_utilities, 1, function(z) {
      if (all(z < 0)) {
        return(1)
      } else {
        return(which.max(z)+1)
      }
    })
    Y_cat = factor(Y, levels = 1:n_choices)
    
    result = list(
      X = covariates,
      Z_rel = relative_utilities,
      Y_cat = Y_cat
    )
  }
}

# Test function to verify 3D array handling
test_3d_array_handling <- function() {
  cat("Testing 3D array handling...\n")
  
  # Test parameters
  n_obs = 50
  n_choices = 3
  p = 2
  x_range = c(-1, 1)
  coef_true = matrix(c(0.5, -0.3), nrow = p)
  Sigma_true = matrix(c(1, 0.3, 0.3, 1), nrow = n_choices-1)
  
  # Test generate_choice_data
  result = generate_choice_data(n_obs, n_choices, x_range, coef_true, 
                               diag(n_choices), seed = 123)
  
  # Check dimensions
  stopifnot(dim(result$X) == c(n_obs, n_choices, p))
  stopifnot(dim(result$Z) == c(n_obs, n_choices))
  stopifnot(dim(result$X_iden) == c(n_obs, n_choices-1, p))
  stopifnot(length(result$Y_cat) == n_obs)
  
  cat("✓ generate_choice_data 3D array handling correct\n")
  
  # Test generate_identified_choice_data
  result2 = generate_identified_choice_data(n_obs, x_range, coef_true, 
                                           Sigma_true, seed = 123)
  
  # Check dimensions
  stopifnot(dim(result2$X) == c(n_obs, n_choices-1, p))
  stopifnot(dim(result2$Z_rel) == c(n_obs, n_choices-1))
  stopifnot(length(result2$Y_cat) == n_obs)
  
  cat("✓ generate_identified_choice_data 3D array handling correct\n")
  
  # Test custom sampler
  custom_sampler <- function(n) rnorm(n, 0, 0.5)
  result3 = generate_custom_identified_choice_data(n_obs, custom_sampler, 
                                                  coef_true, Sigma_true, seed = 123)
  
  stopifnot(dim(result3$X) == c(n_obs, n_choices-1, p))
  stopifnot(dim(result3$Z_rel) == c(n_obs, n_choices-1))
  
  cat("✓ generate_custom_identified_choice_data 3D array handling correct\n")
  
  cat("All 3D array handling tests passed!\n\n")
}

# Performance comparison function
compare_performance <- function() {
  cat("Performance Comparison Tests\n")
  cat("============================\n")
  
  # Create original functions for comparison
  create_original_functions()
  
  # Test parameters
  n_obs_small = 100
  n_obs_large = 1000
  n_choices = 4
  p = 3
  x_range = c(-2, 2)
  coef_true = matrix(rnorm(p), nrow = p)
  Sigma_true = diag(n_choices) + 0.2
  Sigma_iden = Sigma_true[2:n_choices, 2:n_choices] - 
               Sigma_true[2:n_choices, 1] %*% t(Sigma_true[1, 2:n_choices]) / Sigma_true[1,1]
  
  cat("Test parameters:\n")
  cat(sprintf("Small dataset: n_obs=%d, n_choices=%d, p=%d\n", n_obs_small, n_choices, p))
  cat(sprintf("Large dataset: n_obs=%d, n_choices=%d, p=%d\n", n_obs_large, n_choices, p))
  cat("\n")
  
  # Test 1: generate_choice_data
  cat("Testing generate_choice_data...\n")
  
  # Small dataset
  time_orig_small = system.time({
    result_orig_small = generate_choice_data_original(n_obs_small, n_choices, x_range, 
                                                     coef_true, Sigma_true, seed = 42)
  })
  
  time_opt_small = system.time({
    result_opt_small = generate_choice_data(n_obs_small, n_choices, x_range, 
                                           coef_true, Sigma_true, seed = 42)
  })
  
  # Large dataset
  time_orig_large = system.time({
    result_orig_large = generate_choice_data_original(n_obs_large, n_choices, x_range, 
                                                     coef_true, Sigma_true, seed = 42)
  })
  
  time_opt_large = system.time({
    result_opt_large = generate_choice_data(n_obs_large, n_choices, x_range, 
                                           coef_true, Sigma_true, seed = 42)
  })
  
  cat(sprintf("Small dataset (n=%d):\n", n_obs_small))
  cat(sprintf("  Original: %.4f seconds\n", time_orig_small[3]))
  cat(sprintf("  Optimized: %.4f seconds\n", time_opt_small[3]))
  cat(sprintf("  Speedup: %.1fx\n", time_orig_small[3] / time_opt_small[3]))
  
  cat(sprintf("Large dataset (n=%d):\n", n_obs_large))
  cat(sprintf("  Original: %.4f seconds\n", time_orig_large[3]))
  cat(sprintf("  Optimized: %.4f seconds\n", time_opt_large[3]))
  cat(sprintf("  Speedup: %.1fx\n", time_orig_large[3] / time_opt_large[3]))
  cat("\n")
  
  # Test 2: generate_identified_choice_data
  cat("Testing generate_identified_choice_data...\n")
  
  time_orig_small_id = system.time({
    result_orig_small_id = generate_identified_choice_data_original(n_obs_small, x_range, 
                                                                   coef_true, Sigma_iden, seed = 42)
  })
  
  time_opt_small_id = system.time({
    result_opt_small_id = generate_identified_choice_data(n_obs_small, x_range, 
                                                         coef_true, Sigma_iden, seed = 42)
  })
  
  time_orig_large_id = system.time({
    result_orig_large_id = generate_identified_choice_data_original(n_obs_large, x_range, 
                                                                   coef_true, Sigma_iden, seed = 42)
  })
  
  time_opt_large_id = system.time({
    result_opt_large_id = generate_identified_choice_data(n_obs_large, x_range, 
                                                         coef_true, Sigma_iden, seed = 42)
  })
  
  cat(sprintf("Small dataset (n=%d):\n", n_obs_small))
  cat(sprintf("  Original: %.4f seconds\n", time_orig_small_id[3]))
  cat(sprintf("  Optimized: %.4f seconds\n", time_opt_small_id[3]))
  cat(sprintf("  Speedup: %.1fx\n", time_orig_small_id[3] / time_opt_small_id[3]))
  
  cat(sprintf("Large dataset (n=%d):\n", n_obs_large))
  cat(sprintf("  Original: %.4f seconds\n", time_orig_large_id[3]))
  cat(sprintf("  Optimized: %.4f seconds\n", time_opt_large_id[3]))
  cat(sprintf("  Speedup: %.1fx\n", time_orig_large_id[3] / time_opt_large_id[3]))
  cat("\n")
}

# Statistical equivalence test
test_statistical_equivalence <- function() {
  cat("Statistical Equivalence Tests\n")
  cat("=============================\n")
  
  create_original_functions()
  
  # Test parameters
  n_obs = 500
  n_choices = 3
  p = 2
  x_range = c(-1, 1)
  coef_true = matrix(c(0.5, -0.3), nrow = p)
  Sigma_true = matrix(c(1, 0.2, 0.2, 1), nrow = n_choices-1)
  
  # Generate multiple datasets with same seed
  seeds = c(123, 456, 789)
  
  for (seed in seeds) {
    # Original
    result_orig = generate_identified_choice_data_original(n_obs, x_range, 
                                                          coef_true, Sigma_true, seed = seed)
    # Optimized
    result_opt = generate_identified_choice_data(n_obs, x_range, 
                                                coef_true, Sigma_true, seed = seed)
    
    # Test if utilities are statistically similar (should be identical with same seed)
    max_diff = max(abs(result_orig$Z_rel - result_opt$Z_rel))
    cat(sprintf("Seed %d: Max utility difference = %.2e\n", seed, max_diff))
    
    # Test if choices are identical
    choice_match = all(result_orig$Y_cat == result_opt$Y_cat)
    cat(sprintf("Seed %d: Choices identical = %s\n", seed, choice_match))
    
    # The utilities should be identical (or very close due to floating point)
    if (max_diff > 1e-10) {
      warning(sprintf("Large difference in utilities for seed %d: %.2e", seed, max_diff))
    }
  }
  
  cat("Statistical equivalence tests completed.\n\n")
}

# Test fast_mvrnorm function specifically
test_fast_mvrnorm <- function() {
  cat("Testing fast_mvrnorm function...\n")
  
  n_obs = 1000
  n_vars = 3
  
  # Create test means and covariance
  means = matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
  sigma = matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), nrow = n_vars)
  
  # Time the fast function
  time_fast = system.time({
    result_fast = fast_mvrnorm(means, sigma)
  })
  
  # Time using mvtnorm (row by row)
  time_slow = system.time({
    result_slow = t(apply(means, 1, function(mu) rmvnorm(1, mu, sigma)))
  })
  
  cat(sprintf("fast_mvrnorm: %.4f seconds\n", time_fast[3]))
  cat(sprintf("mvtnorm (row-wise): %.4f seconds\n", time_slow[3]))
  cat(sprintf("Speedup: %.1fx\n", time_slow[3] / time_fast[3]))
  
  # Check dimensions
  stopifnot(dim(result_fast) == c(n_obs, n_vars))
  
  cat("✓ fast_mvrnorm test passed\n\n")
}

# Run all tests
main <- function() {
  cat("=== Optimized simdata_tools.R Performance Tests ===\n\n")
  
  test_3d_array_handling()
  test_fast_mvrnorm()
  compare_performance()
  test_statistical_equivalence()
  
  cat("=== All tests completed ===\n")
}

# Run the tests
if (!interactive()) {
  main()
}
