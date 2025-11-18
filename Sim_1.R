library(devtools)
devtools::install_github('mprice747/scrdiff')
library(scrdiff)


# Simulations for 1 Data Point
set.seed(8392013)
# Sample Size
#n <- 50
n <- 100
#n <- 200
#n <- 300
num_sims <- 100

# Posterior HPD Intervals for Each Method
stat_points_1_90_diff_mcmc <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_1_95_diff_mcmc <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_1_99_diff_mcmc <- matrix(rep(0, num_sims * 2), ncol = 2)

stat_points_1_90_dgp  <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_1_95_dgp  <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_1_99_dgp  <- matrix(rep(0, num_sims * 2), ncol = 2)

# MAP of Stationary Point
diff_map <- rep(0, num_sims)
dgp_map <- rep(0, num_sims)

# Posterior Mean of Stationary Point
diff_mean <- rep(0, num_sims)
dgp_mean <- rep(0, num_sims )

# Posterior Variance of Stationary Point
post_var_1_diff_mcmc <- rep(0, num_sims)
post_var_1_dgp <- rep(0, num_sims)


# Parameters
# n = 50, num_betas = 4
# n = 100, num_betas = 5
# n = 200, num_betas = 6
# n = 300, num_betas = 7

# Diffeomorphism Parameters, see example.R to for explanation 
num_betas <- 5
num_stationary <- 1
first_direction <- 1
zero_is_zero <- FALSE
interpolation <- "cubic"
b_vec <- NULL
n_chains <- 100
n_per_chain <- 2500

# n = 50, cov_scale = 0.93
# n = 100, cov_scale = 0.88
# n = 200, cov_scale = 0.8
# n = 300, cov_scale = 0.72
cov_scale <- 0.88
apply_sm_correction <- TRUE
num_cores <- "ALL"

#X_matrix <- matrix(rep(0, num_sims * n), nrow = num_sims)
#Y_matrix <- matrix(rep(0, num_sims * n), nrow = num_sims)

# For DGP method
total_resampled <- 10000

X_matrix <- as.matrix(read.csv("Data/Sim_1/X_matrix_100.csv"))
Y_matrix <- as.matrix(read.csv("Data/Sim_1/Y_matrix_100.csv"))


for (i in 1:num_sims){
  
  print(i)
  
  #0.39973
  #X <- runif(n, 0, 1)
  #Y <- (-2*X^2) + (3*X) + sin(2*X) + cos(3 *X) + 1 + rnorm(n, sd = 0.25)
  
  #X_matrix[i, ] <- X
  #Y_matrix[i, ] <- Y
  
  X <- X_matrix[i, ]
  Y <- Y_matrix[i, ]
  
  prior_mean <- c(rep(0, num_betas), mean(Y),
                  rep(0, num_stationary + 1), 0)
  prior_sd <- c(rep(1, num_betas), 2 * sd(Y), 
                rep(max(Y) - min(Y), num_stationary + 1), 
                0.5)
  
  
  
  # HPD Intervals of Stationary Points using Diffeomorphism Method
  diff_method_mcmc <- posterior_stat_points_diff_mcmc(X, Y, num_betas, num_stationary, first_direction, 
                                                      zero_is_zero, interpolation, b_vec, 
                                                      prior_mean, prior_sd,
                                                      n_chains, n_per_chain, cov_scale,
                                                      apply_sm_correction, 
                                                      num_cores)
  
  # HPD Intervals of Stationary Points using DGP Method
  print("DGP Method")
  dgp_method <- posterior_stat_points_dgp(X, Y, num_stationary, total_resampled)
  
  diff_map[i] <- diff_method_mcmc$diff_map_stat_points
  dgp_map[i] <- dgp_method$dgp_map_stat_points
  
  diff_mean[i] <- mean(diff_method_mcmc$diff_stat_points[[1]])
  dgp_mean[i] <- mean(dgp_method$dgp_stat_points[[1]])
  
  # Get Stationary Point Intervals
  stat_points_1_90_diff_mcmc[i, ] <- diff_method_mcmc$intervals_list_diff$intervals_90[1, ]
  stat_points_1_95_diff_mcmc[i, ] <- diff_method_mcmc$intervals_list_diff$intervals_95[1, ]
  stat_points_1_99_diff_mcmc[i, ] <- diff_method_mcmc$intervals_list_diff$intervals_99[1, ]
  
  stat_points_1_90_dgp[i, ] <- dgp_method$intervals_list_dgp$intervals_90[1, ]
  stat_points_1_95_dgp[i, ] <- dgp_method$intervals_list_dgp$intervals_95[1, ]
  stat_points_1_99_dgp[i, ] <- dgp_method$intervals_list_dgp$intervals_99[1, ] 
  
  # Posterior Variance
  post_var_1_diff_mcmc[i] <- diff_method_mcmc$diff_post_cov[1, ]
  post_var_1_dgp[i] <- dgp_method$dgp_post_var[1]
}

# Report Coverage Probabilities for Each Method
sum((stat_points_1_90_diff_mcmc[, 1] <  0.39973 ) & (stat_points_1_90_diff_mcmc[, 2] >  0.39973))
sum((stat_points_1_95_diff_mcmc[, 1] <  0.39973 ) & (stat_points_1_95_diff_mcmc[, 2] >  0.39973))
sum((stat_points_1_99_diff_mcmc[, 1] <  0.39973 ) & (stat_points_1_99_diff_mcmc[, 2] >  0.39973))

sum((stat_points_1_90_dgp[, 1] <  0.39973 ) & (stat_points_1_90_dgp[, 2] >  0.39973))
sum((stat_points_1_95_dgp[, 1] <  0.39973 ) & (stat_points_1_95_dgp[, 2] >  0.39973))
sum((stat_points_1_99_dgp[, 1] <  0.39973 ) & (stat_points_1_99_dgp[, 2] >  0.39973))

# Report RMSE for MAP Estimate of Stationary Points
sqrt(mean((diff_map - 0.39973)^2))
sqrt(mean((dgp_map - 0.39973)^2))

# Report MSE for Posterior Mean of Stationary Points
sqrt(mean((diff_mean - 0.39973)^2))
sqrt(mean((dgp_mean - 0.39973)^2))

# Report Mean Bias for MAP 
(mean(diff_map - 0.39973))
(mean(dgp_map - 0.39973))

# Report Mean Bias for Posterior Mean 
(mean(diff_mean - 0.39973))
(mean(dgp_mean - 0.39973))

# Report Mean Posterior Standard Deviation
mean(sqrt(post_var_1_diff_mcmc))
mean(sqrt(post_var_1_dgp))

X_matrix_df <- as.data.frame(X_matrix)
Y_matrix_df <- as.data.frame(Y_matrix)

# Data Used (CHANGE)
write.csv(X_matrix_df, file = "Data/Sim_1/X_matrix_300.csv", row.names = FALSE)
write.csv(Y_matrix_df, file = "Data/Sim_1/Y_matrix_300.csv", row.names = FALSE)


# Record HDP Intervals
HDP_df <- as.data.frame(cbind(stat_points_1_90_diff_mcmc, stat_points_1_95_diff_mcmc, 
                              stat_points_1_99_diff_mcmc, stat_points_1_90_dgp, 
                              stat_points_1_95_dgp, stat_points_1_99_dgp ))

names(HDP_df) <- c("SP1_90_Lower_Diff", "SP1_90_Upper_Diff", "SP1_95_Lower_Diff",
                   "SP1_95_Upper_Diff", "SP1_99_Lower_Diff", "SP1_99_Upper_Diff", 
                   "SP1_90_Lower_DGP", "SP1_90_Upper_DGP", "SP1_95_Lower_DGP",
                   "SP1_95_Upper_DGP", "SP1_99_Lower_DGP", "SP1_99_Upper_DGP")

write.csv(HDP_df, file = "Sim_Results/Sim_1/HDP_Ints/HDPs_300.csv", row.names = FALSE)


# Record MAP, Posterior Mean and Posterior Standard Deviation of Stationary Point
metrics_df <- as.data.frame(cbind(diff_map, diff_mean, sqrt(post_var_1_diff_mcmc),
                                  dgp_map, dgp_mean, sqrt(post_var_1_dgp)))

names(metrics_df) <- c("SP1_MAP_Diff", "SP1_Post_Mean_Diff", "SP1_Post_SD_Diff", 
                       "SP1_MAP_DGP", "SP1_Post_Mean_DGP", "SP1_Post_SD_DGP")

write.csv(metrics_df, file = "Sim_Results/Sim_1/Metrics/Metrics_300.csv", row.names = FALSE)
