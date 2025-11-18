library(devtools)
devtools::install_github('mprice747/scrdiff')
library(scrdiff)


set.seed(5903891)
# Sample Size
#n <- 50
n <- 100
#n <- 200
#n <- 300
num_sims <- 100

#Posterior HPD Intervals for Each Stationary Point
stat_points_1_90_diff_mcmc <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_1_95_diff_mcmc <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_1_99_diff_mcmc <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_1_bonf_diff_mcmc <- matrix(rep(0, num_sims * 2), ncol = 2)

stat_points_2_90_diff_mcmc <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_2_95_diff_mcmc <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_2_99_diff_mcmc <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_2_bonf_diff_mcmc <- matrix(rep(0, num_sims * 2), ncol = 2)

# MAP and Posterior Mean for Each Stationary Point
diff_map_1 <- rep(0, num_sims)
diff_mean_1 <- rep(0, num_sims)

diff_map_2 <- rep(0, num_sims)
diff_mean_2 <- rep(0, num_sims)

# Posterior Variance for Each Stationary Point
post_var_1_diff_mcmc <- rep(0, num_sims)
post_var_2_diff_mcmc <- rep(0, num_sims)

# Parameters
# n = 50, num_betas = 7
# n = 100, num_betas = 8
# n = 200, num_betas = 9
# n = 300, num_betas = 10

# Diffeomorphism Parameters, see example.R to for explanation
num_betas <- 8
num_stationary <- 2
first_direction <- 1
zero_is_zero <- FALSE
interpolation <- "cubic"
b_vec <- NULL
n_chains <- 100
n_per_chain <- 2500

# n = 50, cov_scale  = 0.7
# n = 100, cov_scale = 0.63
# n = 200, cov_scale = 0.57
# n = 300, cov_scale = 0.5
cov_scale <- 0.63
apply_sm_correction <- TRUE
num_cores <- "ALL"

#X_matrix <- matrix(rep(0, num_sims * n), nrow = num_sims)
#Y_matrix <- matrix(rep(0, num_sims * n), nrow = num_sims)

# Load Data
X_matrix <- as.matrix(read.csv("Data/Sim_2/X_matrix_100.csv"))
Y_matrix <- as.matrix(read.csv("Data/Sim_2/Y_matrix_100.csv"))

for (i in 1:num_sims){
  
  print(i)
  
  #0.24204, 0.687345
  #X <- runif(n, 0, 1)
  #Y <- sin(1.1* 2 *X) + cos(4*X) + 0.7*2*X + rnorm(n, sd = 0.15)
  
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
  

  
  # Get Stationary Point Intervals
  stat_points_1_90_diff_mcmc[i, ] <- diff_method_mcmc$intervals_list_diff$intervals_90[1, ]
  stat_points_1_95_diff_mcmc[i, ] <- diff_method_mcmc$intervals_list_diff$intervals_95[1, ]
  stat_points_1_99_diff_mcmc[i, ] <- diff_method_mcmc$intervals_list_diff$intervals_99[1, ]
  stat_points_1_bonf_diff_mcmc[i, ] <- diff_method_mcmc$intervals_list_diff$intervals_bonf[1, ]
  
  stat_points_2_90_diff_mcmc[i, ] <- diff_method_mcmc$intervals_list_diff$intervals_90[2, ]
  stat_points_2_95_diff_mcmc[i, ] <- diff_method_mcmc$intervals_list_diff$intervals_95[2, ]
  stat_points_2_99_diff_mcmc[i, ] <- diff_method_mcmc$intervals_list_diff$intervals_99[2, ]
  stat_points_2_bonf_diff_mcmc[i, ] <- diff_method_mcmc$intervals_list_diff$intervals_bonf[2, ]
  
  # Posterior Means and MAPs for Stationary Points
  diff_map_1[i] <- diff_method_mcmc$diff_map_stat_points[1]
  diff_mean_1[i] <- mean(diff_method_mcmc$diff_stat_points[[1]])
  
  diff_map_2[i] <- diff_method_mcmc$diff_map_stat_points[2]
  diff_mean_2[i] <- mean(diff_method_mcmc$diff_stat_points[[2]])
  
  post_var_1_diff_mcmc[i] <- diff_method_mcmc$diff_post_cov[1, 1]
  post_var_2_diff_mcmc[i] <- diff_method_mcmc$diff_post_cov[2, 2]
}

# Report Coverage Probabilities for Each Stationary Point
sum((stat_points_1_90_diff_mcmc[, 1] < 0.24204 ) & (stat_points_1_90_diff_mcmc[, 2] > 0.24204))
sum((stat_points_1_95_diff_mcmc[, 1] < 0.24204) & (stat_points_1_95_diff_mcmc[, 2] > 0.24204))
sum((stat_points_1_99_diff_mcmc[, 1] < 0.24204 ) & (stat_points_1_99_diff_mcmc[, 2] > 0.24204))

sum((stat_points_2_90_diff_mcmc[, 1] < 0.687345) & (stat_points_2_90_diff_mcmc[, 2] > 0.687345))
sum((stat_points_2_95_diff_mcmc[, 1] < 0.687345) & (stat_points_2_95_diff_mcmc[, 2] > 0.687345))
sum((stat_points_2_99_diff_mcmc[, 1] < 0.687345) & (stat_points_2_99_diff_mcmc[, 2] >  0.687345))

# Bonferonni Coverage
sum((stat_points_1_bonf_diff_mcmc[, 1] < 0.24204 ) & (stat_points_1_bonf_diff_mcmc[, 2] > 0.24204) & (stat_points_2_bonf_diff_mcmc[, 1] < 0.687345) & (stat_points_2_bonf_diff_mcmc[, 2] > 0.687345))

# RMSE
sqrt(mean((diff_map_1 -0.24204)^2))
sqrt(mean((diff_mean_1 -0.24204)^2))

sqrt(mean((diff_map_2 - 0.687345)^2))
sqrt(mean((diff_mean_2 - 0.687345)^2))

# Bias
(mean(diff_map_1 -0.24204))
(mean(diff_mean_1 -0.24204))

(mean(diff_map_2 - 0.687345))
(mean(diff_mean_2 - 0.687345))

# Posterior Standard Deviation
mean(sqrt(post_var_1_diff_mcmc))
mean(sqrt(post_var_2_diff_mcmc))

X_matrix_df <- as.data.frame(X_matrix)
Y_matrix_df <- as.data.frame(Y_matrix)

# Data Used
write.csv(X_matrix_df, file = "Data/Sim_2/X_matrix_300.csv", row.names = FALSE)
write.csv(Y_matrix_df, file = "Data/Sim_2/Y_matrix_300.csv", row.names = FALSE)

# Record HDP Intervals
HDP_df <- as.data.frame(cbind(stat_points_1_90_diff_mcmc, stat_points_1_95_diff_mcmc, 
                              stat_points_1_99_diff_mcmc, stat_points_1_bonf_diff_mcmc, 
                              stat_points_2_90_diff_mcmc, stat_points_2_95_diff_mcmc, 
                              stat_points_2_99_diff_mcmc, stat_points_2_bonf_diff_mcmc))

names(HDP_df) <- c("SP1_90_Lower", "SP1_90_Upper", "SP1_95_Lower", "SP1_95_Upper", 
                   "SP1_99_Lower", "SP1_99_Upper", "SP1_Bonf_Lower", "SP1_Bonf_Upper", 
                   "SP2_90_Lower", "SP2_90_Upper", "SP2_95_Lower", "SP2_95_Upper", 
                   "SP2_99_Lower", "SP2_99_Upper", "SP2_Bonf_Lower", "SP2_Bonf_Upper")

write.csv(HDP_df, file = "Sim_Results/Sim_2/HDP_Ints/HDPs_300.csv", row.names = FALSE)

# Record MAP, Posterior Mean and Posterior Standard Deviation of Stationary Points
metrics_df <- as.data.frame(cbind(diff_map_1, diff_mean_1, sqrt(post_var_1_diff_mcmc),
                                  diff_map_2, diff_mean_2, sqrt(post_var_2_diff_mcmc)))

names(metrics_df) <- c("SP1_MAP", "SP1_Post_Mean", "SP1_Post_SD", 
                       "SP2_MAP", "SP2_Post_Mean", "SP2_Post_SD")

write.csv(metrics_df, file = "Sim_Results/Sim_2/Metrics/Metrics_300.csv", row.names = FALSE)

