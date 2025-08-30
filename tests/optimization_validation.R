# Validation Test Script for Optimized simdata_tools.R
# This script tests edge cases, parameter combinations, and correctness

library(mvtnorm)
library(CVXR)
source("R/simdata_tools.R")

# Test edge cases and parameter combinations
test_edge_cases <- function() {
  cat("Testing edge cases and parameter combinations...\n")
  
  # Test 1: Small datasets
  cat("✓ Testing small datasets...\n")
  n_obs = 5
  n_choices = 2
  p = 1
  x_range = c(0, 1)
  coef_true = matrix(0.5, nrow = p)
  Sigma_true = diag(n_choices)  # Fixed: proper square matrix
  
  result = generate_choice_data(n_obs, n_choices, x_range, coef_true, Sigma_true, seed = 1)
  stopifnot(dim(result$X) == c(n_obs, n_choices, p))
  stopifnot(dim(result$Z) == c(n_obs, n_choices))
  
  # Test 2: Larger number of choices
  cat("✓ Testing larger number of choices...\n")
  n_obs = 100
  n_choices = 6
  p = 2
  coef_true = matrix(c(0.3, -0.2), nrow = p)
  Sigma_true = diag(n_choices) + 0.1
  
  result = generate_choice_data(n_obs, n_choices, x_range, coef_true, Sigma_true, seed = 2)
  stopifnot(dim(result$X) == c(n_obs, n_choices, p))
  stopifnot(dim(result$Z) == c(n_obs, n_choices))
  
  # Test 3: Higher dimensional covariates
  cat("✓ Testing higher dimensional covariates...\n")
  n_obs = 50
  n_choices = 3
  p = 5
  coef_true = matrix(rnorm(p), nrow = p)
  Sigma_true = diag(n_choices)
  
  result = generate_choice_data(n_obs, n_choices, x_range, coef_true, Sigma_true, seed = 3)
  stopifnot(dim(result$X) == c(n_obs, n_choices, p))
  stopifnot(dim(result$Z) == c(n_obs, n_choices))
  
  # Test 4: Different covariance structures
  cat("✓ Testing different covariance structures...\n")
  n_obs = 100
  n_choices = 4
  p = 2
  coef_true = matrix(c(1, -1), nrow = p)
  
  # High correlation
  Sigma_high_corr = matrix(c(1, 0.8, 0.8, 0.8,
                            0.8, 1, 0.8, 0.8,
                            0.8, 0.8, 1, 0.8,
                            0.8, 0.8, 0.8, 1), nrow = n_choices)
  
  result = generate_choice_data(n_obs, n_choices, x_range, coef_true, Sigma_high_corr, seed = 4)
  stopifnot(all(is.finite(result$Z)))
  
  cat("All edge case tests passed!\n\n")
}

# Test random seed consistency
test_seed_consistency <- function() {
  cat("Testing random seed consistency...\n")
  
  n_obs = 100
  n_choices = 3
  p = 2
  x_range = c(-1, 1)
  coef_true = matrix(c(0.5, -0.3), nrow = p)
  Sigma_true = diag(n_choices)
  
  # Generate same data multiple times with same seed
  seeds = c(123, 456, 789)
  
  for (seed in seeds) {
    result1 = generate_choice_data(n_obs, n_choices, x_range, coef_true, Sigma_true, seed = seed)
    result2 = generate_choice_data(n_obs, n_choices, x_range, coef_true, Sigma_true, seed = seed)
    
    # Check that results are identical
    stopifnot(all.equal(result1$X, result2$X))
    stopifnot(all.equal(result1$Z, result2$Z))
    stopifnot(all.equal(result1$Y_cat, result2$Y_cat))
    
    cat(sprintf("✓ Seed %d: Results are consistent\n", seed))
  }
  
  # Test different seeds produce different results
  result_seed1 = generate_choice_data(n_obs, n_choices, x_range, coef_true, Sigma_true, seed = 111)
  result_seed2 = generate_choice_data(n_obs, n_choices, x_range, coef_true, Sigma_true, seed = 222)
  
  # Results should be different
  are_equal = isTRUE(all.equal(result_seed1$Z, result_seed2$Z))
  stopifnot(!are_equal)
  cat("✓ Different seeds produce different results\n\n")
}

# Test custom sampler functionality
test_custom_sampler <- function() {
  cat("Testing custom sampler functionality...\n")
  
  n_obs = 100
  p = 2
  coef_true = matrix(c(0.5, -0.3), nrow = p)
  Sigma_iden = matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  
  # Test 1: Normal sampler
  normal_sampler <- function(n) rnorm(n, 0, 1)
  result1 = generate_custom_identified_choice_data(n_obs, normal_sampler, 
                                                  coef_true, Sigma_iden, seed = 1)
  
  # Test 2: Uniform sampler
  uniform_sampler <- function(n) runif(n, -2, 2)
  result2 = generate_custom_identified_choice_data(n_obs, uniform_sampler, 
                                                  coef_true, Sigma_iden, seed = 1)
  
  # Test 3: Exponential sampler
  exp_sampler <- function(n) rexp(n, 1)
  result3 = generate_custom_identified_choice_data(n_obs, exp_sampler, 
                                                  coef_true, Sigma_iden, seed = 1)
  
  # Check dimensions and finite values
  for (result in list(result1, result2, result3)) {
    stopifnot(dim(result$X) == c(n_obs, 2, p))
    stopifnot(dim(result$Z_rel) == c(n_obs, 2))
    stopifnot(all(is.finite(result$Z_rel)))
    stopifnot(length(result$Y_cat) == n_obs)
  }
  
  # Results should be different due to different samplers
  are_equal_12 = isTRUE(all.equal(result1$X, result2$X))
  are_equal_23 = isTRUE(all.equal(result2$X, result3$X))
  stopifnot(!are_equal_12)
  stopifnot(!are_equal_23)
  
  cat("✓ Custom sampler tests passed\n\n")
}

# Test target identified choice data generation
test_target_identified <- function() {
  cat("Testing target identified choice data generation...\n")
  
  n_obs = 50
  x_range = c(-1, 1)
  p = 2
  coef_true = matrix(c(0.5, -0.3), nrow = p)
  
  # Create a test precision matrix
  Prec_iden = matrix(c(2, 0.5, 0.5, 1.5), nrow = 2)
  
  # Test the function
  result = generate_target_identified_choice_data(n_obs, x_range, coef_true, 
                                                 Prec_iden, seed = 1)
  
  # Check dimensions
  n_choices = nrow(Prec_iden) + 1
  stopifnot(dim(result$X) == c(n_obs, n_choices, p))
  stopifnot(dim(result$Z) == c(n_obs, n_choices))
  stopifnot(dim(result$X_iden) == c(n_obs, n_choices-1, p))
  stopifnot(length(result$Y_cat) == n_obs)
  
  # Check that precision matrix is valid
  stopifnot(all(is.finite(result$Prec_not_iden)))
  stopifnot(dim(result$Prec_not_iden) == c(n_choices, n_choices))
  
  cat("✓ Target identified choice data generation test passed\n\n")
}

# Test rchoicemvn function
test_rchoicemvn <- function() {
  cat("Testing rchoicemvn function...\n")
  
  n_obs = 100
  n_choices = 3
  
  # Create test means and precision matrix
  means = matrix(rnorm(n_obs * n_choices), nrow = n_obs, ncol = n_choices)
  precision = solve(matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), nrow = n_choices))
  
  # Test the function
  result = rchoicemvn(means, precision)
  
  # Check dimensions
  stopifnot(dim(result) == c(n_obs, n_choices))
  stopifnot(all(is.finite(result)))
  
  # Test with different precision matrices
  # Identity precision (independent normals)
  precision_id = diag(n_choices)
  result_id = rchoicemvn(means, precision_id)
  stopifnot(dim(result_id) == c(n_obs, n_choices))
  
  cat("✓ rchoicemvn function test passed\n\n")
}

# Test choice probability distributions
test_choice_probabilities <- function() {
  cat("Testing choice probability distributions...\n")
  
  n_obs = 1000
  n_choices = 3
  p = 1
  x_range = c(0, 0)  # Constant covariates
  coef_true = matrix(c(1), nrow = p)
  
  # Test 1: Equal utility case (should have roughly equal choice probabilities)
  Sigma_equal = diag(n_choices)
  result_equal = generate_choice_data(n_obs, n_choices, x_range, 
                                     matrix(0, nrow = p), Sigma_equal, seed = 1)
  
  choice_props = table(result_equal$Y_cat) / n_obs
  cat(sprintf("Equal utility case - Choice proportions: %s\n", 
              paste(round(choice_props, 3), collapse = ", ")))
  
  # Test 2: Unequal utility case
  coef_unequal = matrix(c(2), nrow = p)
  x_range_unequal = c(-1, 1)
  result_unequal = generate_choice_data(n_obs, n_choices, x_range_unequal, 
                                       coef_unequal, Sigma_equal, seed = 1)
  
  choice_props_unequal = table(result_unequal$Y_cat) / n_obs
  cat(sprintf("Unequal utility case - Choice proportions: %s\n", 
              paste(round(choice_props_unequal, 3), collapse = ", ")))
  
  cat("✓ Choice probability distribution test completed\n\n")
}

# Comprehensive validation function
validate_optimizations <- function() {
  cat("=== Validation Tests for Optimized simdata_tools.R ===\n\n")
  
  test_edge_cases()
  test_seed_consistency()
  test_custom_sampler()
  test_target_identified()
  test_rchoicemvn()
  test_choice_probabilities()
  
  cat("=== All validation tests completed successfully ===\n")
}

# Run validation if script is executed directly
if (!interactive()) {
  validate_optimizations()
}
