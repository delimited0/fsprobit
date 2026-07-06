# Test script for mnp_ops.R functions
# Tests each function with one random test example per case
# Fixed seed = 1 for reproducibility

library(fsprobit)

# Source the mnp_ops.R file to load the functions
source("R/mnp_ops.R")

# Set fixed seed for reproducibility
set.seed(1)

cat("Testing mnp_ops.R functions\n")
cat("===========================\n\n")

# Test 1: vec_relative_choice
cat("Test 1: vec_relative_choice\n")
cat("---------------------------\n")

# Create test vector
m = 20
p = 2
v_test = 1:(m-1)
constraints = utility_shift_constraints(m)
choice_idx = 2

A = constraints[[choice_idx]]
A %*% v_test 
vec_relative_choice(choice_idx-1, v_test)

# Apply function
result_vec <- vec_relative_choice(choice_idx, v_test)

cov_A = matrix(rnorm((m-1)*(m-1)), nrow = m-1)
covmat = cov_A %*% t(cov_A)
M_test = covmat

A %*% M_test
rows_relative_choice(choice_idx-1, M_test)

A %*% M_test %*% t(A)
sandwich_choice(choice_idx-1+10, M_test)