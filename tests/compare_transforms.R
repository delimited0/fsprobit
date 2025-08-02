library(fsprobit)
library(mvtnorm)
library(ggplot2)
library(data.table)

# Simulation Parameters
n_choices <- 25
n_obs <- 2000
max_iter = 500
tol = 1e-3

newton_tol = 1e-2
max_newton_iter = 50

# True Covariance/Precision Matrix (Compound Symmetric)
# User requested compound symmetric covariance = 0.5
# This implies a covariance matrix with 1s on the diagonal and 0.5 on the off-diagonal
off_diag_val <- 0.5
Sigma_iden <- matrix(off_diag_val, nrow = n_choices - 1, ncol = n_choices - 1)
diag(Sigma_iden) <- 1
Prec_iden <- solve(Sigma_iden)
coef_true = as.matrix(c(1,-1))

# Generate Data
simdata <- generate_custom_identified_choice_data(
  n_obs = n_obs,
  sampler = function(n) rnorm(n, mean = 0, sd = 1),
  coef_true = coef_true,
  Sigma_iden = Sigma_iden,
  seed = 12356
)

# conv_metric = "beta_relative"
conv_metric = "precision"

# Initial parameters for the model
sigma_init <- diag(n_choices - 1)
coef_init <- as.matrix(c(0, 0))

cl <- makeCluster(8) # Example: use all but one core
registerDoParallel(cl)

# Fit Model 1: With transform
probit_ep_transformed <- mnp_probit(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = sigma_init,
  E_method = "EP",
  M_method = "Newton",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  max_iter = max_iter,
  newton_tol = newton_tol,
  max_newton_iter = max_newton_iter,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  verbose = 5,
  record_history = TRUE,
  conv_metric = conv_metric,
  transform = TRUE
)
stopCluster(cl)


cl <- makeCluster(8) # Example: use all but one core
registerDoParallel(cl)
probit_ep_cvx = mnp_probit(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = sigma_init,
  E_method = "EP",
  M_method = "CVX",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  max_iter = max_iter,
  newton_tol = newton_tol,
  max_newton_iter = max_newton_iter,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  verbose = 5,
  record_history = TRUE,
  conv_metric = conv_metric,
  transform = TRUE
)
stopCluster(cl)


cl <- makeCluster(8) # Example: use all but one core
registerDoParallel(cl)
# Fit Model 2: Without transform
probit_ep_untransformed <- mnp_probit(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = sigma_init,
  E_method = "EP",
  M_method = "Newton",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  max_iter = max_iter,
  newton_tol = newton_tol,
  max_newton_iter = max_newton_iter,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  verbose = 5,
  conv_metric = conv_metric,
  transform = FALSE
)

stopCluster(cl)

cl <- makeCluster(8) # Example: use all but one core
registerDoParallel(cl)
# Fit Model 2: Without transform
probit_hmc_untransformed <- mnp_probit(
  X = simdata$X, Y = simdata$Y,
  beta_init = coef_init,
  Sigma_init = sigma_init,
  E_method = "HMC",
  M_method = "Newton",
  n_choices = n_choices,
  true_trace = sum(diag(Prec_iden)),
  tol = tol,
  max_iter = max_iter,
  newton_tol = newton_tol,
  max_newton_iter = max_newton_iter,
  shift_iden_method = "ref",
  scale_iden_method = "trace",
  verbose = 5,
  n_mc = 1000,
  conv_metric = conv_metric,
  transform = FALSE
)

stopCluster(cl)

###### visualization ########

# Extract the estimated precision matrices
prec_transformed <- probit_ep_transformed$Precision
prec_untransformed <- probit_ep_untransformed$Precision

# Prepare data for ggplot by melting the matrices into a long format
melt_prec <- function(prec_matrix, model_name) {
  dt <- as.data.table(prec_matrix)
  dt[, Var1 := .I]
  melted_dt <- melt(dt, id.vars = "Var1", variable.name = "Var2", value.name = "value")
  melted_dt[, Var2 := as.integer(sub("V", "", Var2))]
  melted_dt[, model := model_name]
  melted_dt[, type := ifelse(Var1 == Var2, "Diagonal", "Off-Diagonal")]
  return(melted_dt)
}

df_transformed <- melt_prec(prec_transformed, "Transformed")
df_untransformed <- melt_prec(prec_untransformed, "Untransformed")

# Combine the data from both models
plot_data <- rbindlist(list(df_transformed, df_untransformed))

# Create a data frame for the true values to draw horizontal lines
true_vals_df <- data.frame(
    type = c("Diagonal", "Off-Diagonal"),
    true_value = c(diag(true_prec)[1], true_prec[1, 2])
)

# Create the ggplot visualization
comparison_plot <- ggplot(plot_data, aes(x = model, y = value, fill = model)) +
  geom_boxplot() +
  facet_wrap(~type, scales = "free_y") +
  geom_hline(data = true_vals_df, aes(yintercept = true_value), color = "black", linetype = "dashed", size = 1) +
  labs(
    title = "Comparison of Estimated Precision Matrix Elements",
    subtitle = "EP Method: With vs. Without Transformation",
    x = "Model Type",
    y = "Estimated Precision Value"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 12))

# Print the plot to the viewer
print(comparison_plot)

# Optionally, save the plot to a file
# ggsave("tests/precision_comparison.png", plot = comparison_plot, width = 8, height = 6)

#$####### training history.  #####
probit_ep_transformed$llik |> plot(main = "Log-Likelihood History")
probit_transformed$m_step_bound

probit_ep_transformed$beta_history[, 1] |> plot(type = "l", main = "Beta History")
probit_ep_transformed$beta_history[, 2] |> plot(type = "l", main = "Beta History")

# Create a faceted plot for the history of selected precision matrix elements
prec_history <- probit_ep_transformed$Prec_history
d <- dim(prec_history)[2]

# Get all off-diagonal indices (upper triangle)
off_diag_indices <- which(upper.tri(matrix(nrow = d, ncol = d)), arr.ind = TRUE)

# Randomly sample 5 of these indices
set.seed(42) # for reproducibility
h = 4
selected_indices <- off_diag_indices[sample(nrow(off_diag_indices), h), ]

# Prepare data for ggplot
history_list <- lapply(1:nrow(selected_indices), function(k) {
  i <- selected_indices[k, "row"]
  j <- selected_indices[k, "col"]
  data.table(
    iteration = 1:dim(prec_history)[1],
    value = prec_history[, i, j],
    element = paste0("Prec[", i, ",", j, "]"),
    true_value = true_prec[i, j]
  )
})

history_df <- rbindlist(history_list)

# Create the ggplot
precision_history_plot <- ggplot(history_df, aes(x = iteration, y = value)) +
  geom_line(aes(color = element)) +
  geom_hline(aes(yintercept = true_value), linetype = "dashed", color = "black") +
  facet_wrap(~element, scales = "free_y", ncol = 2) +
  labs(
    title = "Convergence of Randomly Selected Off-Diagonal Precision Elements",
    x = "Iteration",
    y = "Precision Value"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(precision_history_plot)


