Raw_ERP <- read.csv("Raw_ERP.csv")

all_signals <- Raw_ERP[, 7:256]
names(all_signals)

X <- seq(-100, 896, by = 4)
Y <- as.vector(colMeans(all_signals))

num_betas <- 15
num_stationary <- 10
first_direction <- 1
zero_is_zero <- FALSE
interpolation <- "cubic"
n_chains <- 100
n_per_chain <- 3000
cov_scale <- 1/6


prior_mean <- c(rep(0, num_betas), mean(Y),
                rep(0, num_stationary + 1), 0)
prior_sd <- c(rep(1, num_betas), 2 * sd(Y), 
              rep(max(Y) - min(Y), num_stationary + 1), 
              0.5)


new_b_vec <- seq(min(X), max(X), length.out =12)
new_b_vec[2:11] <- c(0.2, 0.27, 0.38, 0.5, 0.68, 0.75, 0.8, 0.85, 0.9, 0.95)

diff_method_erp <- posterior_stat_points_diff_mcmc(X, Y, num_betas, num_stationary,
                                                   first_direction, 
                                zero_is_zero, interpolation, new_b_vec, 
                                prior_mean, prior_sd,
                                n_chains, n_per_chain, cov_scale,
                                apply_sm_correction = TRUE)

dgp_method <- posterior_stat_points_dgp(X, Y, num_stationary, 100000)


bonf_intervals_dgp <- get_hpd_interval_from_den(dgp_method$post_prob_gp, 
                                                seq(min(X), max(X), length.out = 5000), 
                                                target_prob = 0.995)

map_post_process <- process_post_samples(diff_method_erp$diff_map_par, num_betas,
                                         first_direction)

X_for_pred <- seq(-100, 896, length.out = 10000)

map_predictions <- diffeo_predictions(X_for_pred, 
                                      map_post_process$betas, 
                                      new_b_vec,
                   map_post_process$lambdas)$predictions[, 1]

erp_s1 <- density(diff_method_erp$diff_stat_points[[1]])
erp_s2 <- density(diff_method_erp$diff_stat_points[[2]])
erp_s3 <- density(diff_method_erp$diff_stat_points[[3]])
erp_s4 <- density(diff_method_erp$diff_stat_points[[4]])
erp_s5 <- density(diff_method_erp$diff_stat_points[[5]])
erp_s6 <- density(diff_method_erp$diff_stat_points[[6]])
erp_s7 <- density(diff_method_erp$diff_stat_points[[7]])
erp_s8 <- density(diff_method_erp$diff_stat_points[[8]])
erp_s9 <- density(diff_method_erp$diff_stat_points[[9]])
erp_s10 <- density(diff_method_erp$diff_stat_points[[10]])

bonf_intervals <- diff_method_erp$intervals_list_diff$intervals_bonf

col_list <- c(
  "#1B4F72", # dark blue
  "#8B2500", # dark orange/red
  "#006400", # dark green
  "#A35C00", # dark goldenrod/orange
  "#7B1F5C", # dark magenta/purple
  "#000000", # black
  "#665C00", # dark gold/olive
  "#174E64", # dark teal/blue
  "#5B0A2D", # dark maroon/brown
  "#004D40"  # dark cyan/teal
)

par(mfrow = c(2, 1))

plot(X, Y,  pch = 19, cex = 0.5, xlab = "MSEC", ylab = "ERP Signal", 
     main = "ERP Data Fit with Stationary Point HPDs", ylim = c(-11, 9), col = "gray")
lines(X_for_pred, map_predictions,  col = "coral", lwd = 3)


for (j in 1:10){
  
  in_interval <- (X_for_pred >= bonf_intervals[j, 1]) & (X_for_pred <= bonf_intervals[j, 2])
  sub_X <- X_for_pred[in_interval]
  sub_predictions <- map_predictions[which(in_interval)]
  
  lines(sub_X, sub_predictions, col = col_list[j], lwd = 5)
}

plot(erp_s1$x, erp_s1$y, type = 'l', col = col_list[1], lwd = 3, xlim = c(-100, 896), 
     ylim = c(0, 0.33), main = "Stationary Points Posteriors", xlab = "MSEC", 
     ylab = "Posterior Likelihood")
lines(erp_s2$x, erp_s2$y, type = 'l', col = col_list[2], lwd = 3)
lines(erp_s3$x, erp_s3$y, type = 'l', col = col_list[3], lwd = 3)
lines(erp_s4$x, erp_s4$y, type = 'l', col = col_list[4], lwd = 3)
lines(erp_s5$x, erp_s5$y, type = 'l', col = col_list[5], lwd = 3)
lines(erp_s6$x, erp_s6$y, type = 'l', col = col_list[6], lwd = 3)
lines(erp_s7$x, erp_s7$y, type = 'l', col = col_list[7], lwd = 3)
lines(erp_s8$x, erp_s8$y, type = 'l', col = col_list[8], lwd = 3)
lines(erp_s9$x, erp_s9$y, type = 'l', col = col_list[9], lwd = 3)
lines(erp_s10$x, erp_s10$y, type = 'l', col = col_list[10], lwd = 3)
lines(c(-100, 896), c(0, 0), col = col_list[1], lwd = 3)



par(mfrow = c(1, 1))



map_predictions_1 <- diffeo_predictions(sort(diff_method_erp$diff_stat_points[[1]]), 
                                      map_post_process$betas, 
                                      new_b_vec,
                                      map_post_process$lambdas)

lines(sort(diff_method_erp$diff_stat_points[[1]]), map_predictions_1$predictions[, 1],
      col = 'blue', lwd =2)


segments(bonf_intervals[1, 1], -12, bonf_intervals[1, 1], 10, col = "blue", lwd = 2.5)
segments(bonf_intervals[1, 2], -12, bonf_intervals[1, 2], 10, col = "blue", lwd = 2.5)

segments(bonf_intervals[2, 1], -12, bonf_intervals[2, 1], 10, col = "red", lwd = 2.5)
segments(bonf_intervals[2, 2], -12, bonf_intervals[2, 2], 10, col = "red", lwd = 2.5)

segments(bonf_intervals[3, 1], -12, bonf_intervals[3, 1], 10, col = "green", lwd = 2.5)
segments(bonf_intervals[3, 2], -12, bonf_intervals[3, 2], 10, col = "green", lwd = 2.5)

segments(bonf_intervals[4, 1], -12, bonf_intervals[4, 1], 10, col = "orange", lwd = 2.5)
segments(bonf_intervals[4, 2], -12, bonf_intervals[4, 2], 10, col = "orange", lwd = 2.5)

segments(bonf_intervals[5, 1], -12, bonf_intervals[5, 1], 10, col = "purple", lwd = 2.5)
segments(bonf_intervals[5, 2], -12, bonf_intervals[5, 2], 10, col = "purple", lwd = 2.5)

segments(bonf_intervals[6, 1], -12, bonf_intervals[6, 1], 10, col = "black", lwd = 2.5)
segments(bonf_intervals[6, 2], -12, bonf_intervals[6, 2], 10, col = "black", lwd = 2.5)

segments(bonf_intervals[7, 1], -12, bonf_intervals[7, 1], 10, col = "gold3", lwd = 2.5)
segments(bonf_intervals[7, 2], -12, bonf_intervals[7, 2], 10, col = "gold3", lwd = 2.5)

segments(bonf_intervals[8, 1], -12, bonf_intervals[8, 1], 10, col = "brown", lwd = 2.5)
segments(bonf_intervals[8, 2], -12, bonf_intervals[8, 2], 10, col = "brown", lwd = 2.5)

segments(bonf_intervals[9, 1], -12, bonf_intervals[9, 1], 10, col = "pink", lwd = 2.5)
segments(bonf_intervals[9, 2], -12, bonf_intervals[9, 2], 10, col = "pink", lwd = 2.5)

segments(bonf_intervals[10, 1], -12, bonf_intervals[10, 1], 10, col = "green4", lwd = 2.5)
segments(bonf_intervals[10, 2], -12, bonf_intervals[10, 2], 10, col = "green4", lwd = 2.5)





segments(bonf_intervals_dgp$ci_lower[1], -12, bonf_intervals_dgp$ci_lower[1], 
         10, col = "blue", lwd = 2.5)
segments(bonf_intervals_dgp$ci_upper[1], -12, bonf_intervals_dgp$ci_upper[1], 
         10, col = "blue", lwd = 2.5)

segments(bonf_intervals_dgp$ci_lower[2], -12, bonf_intervals_dgp$ci_lower[2], 
         10, col = "red", lwd = 2.5)
segments(bonf_intervals_dgp$ci_upper[2], -12, bonf_intervals_dgp$ci_upper[2], 
         10, col = "red", lwd = 2.5)

segments(bonf_intervals_dgp$ci_lower[3], -12, bonf_intervals_dgp$ci_lower[3], 
         10, col = "green", lwd = 2.5)
segments(bonf_intervals_dgp$ci_upper[3], -12, bonf_intervals_dgp$ci_upper[3], 
         10, col = "green", lwd = 2.5)

segments(bonf_intervals_dgp$ci_lower[4], -12, bonf_intervals_dgp$ci_lower[4], 
         10, col = "orange", lwd = 2.5)
segments(bonf_intervals_dgp$ci_upper[4], -12, bonf_intervals_dgp$ci_upper[4], 
         10, col = "orange", lwd = 2.5)

segments(bonf_intervals_dgp$ci_lower[5], -12, bonf_intervals_dgp$ci_lower[5], 
         10, col = "purple", lwd = 2.5)
segments(bonf_intervals_dgp$ci_upper[5], -12, bonf_intervals_dgp$ci_upper[5], 
         10, col = "purple", lwd = 2.5)

segments(bonf_intervals_dgp$ci_lower[6], -12, bonf_intervals_dgp$ci_lower[6], 
         10, col = "black", lwd = 2.5)
segments(bonf_intervals_dgp$ci_upper[6], -12, bonf_intervals_dgp$ci_upper[6], 
         10, col = "black", lwd = 2.5)

segments(bonf_intervals_dgp$ci_lower[7], -12, bonf_intervals_dgp$ci_lower[7], 
         10, col = "gold3", lwd = 2.5)
segments(bonf_intervals_dgp$ci_upper[7], -12, bonf_intervals_dgp$ci_upper[7], 
         10, col = "gold3", lwd = 2.5)

segments(bonf_intervals_dgp$ci_lower[8], -12, bonf_intervals_dgp$ci_lower[8], 
         10, col = "brown", lwd = 2.5)
segments(bonf_intervals_dgp$ci_upper[8], -12, bonf_intervals_dgp$ci_upper[8], 
         10, col = "brown", lwd = 2.5)

segments(bonf_intervals_dgp$ci_lower[9], -12, bonf_intervals_dgp$ci_lower[9], 
         10, col = "pink", lwd = 2.5)
segments(bonf_intervals_dgp$ci_upper[9], -12, bonf_intervals_dgp$ci_upper[9], 
         10, col = "pink", lwd = 2.5)

segments(bonf_intervals_dgp$ci_lower[10], -12, bonf_intervals_dgp$ci_lower[10], 
         10, col = "green4", lwd = 2.5)
segments(bonf_intervals_dgp$ci_upper[10], -12, bonf_intervals_dgp$ci_upper[10], 
         10, col = "green4", lwd = 2.5)
