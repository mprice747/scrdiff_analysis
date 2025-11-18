library(devtools)
devtools::install_github('mprice747/scrdiff')
library(scrdiff)

# Data from http://dsenturk.bol.ucla.edu/supplements.html
Raw_ERP <- read.csv("Data/Applied_Data_Analysis/Raw_ERP.csv")

# Get Means for ERP
all_signals <- Raw_ERP[, 7:256]
X <- seq(-100, 896, by = 4)
Y <- as.vector(colMeans(all_signals))

plot(X, Y)

# Diffeomorphism Model Parameters
num_betas <- 15
num_stationary <- 10
first_direction <- 1
zero_is_zero <- FALSE
interpolation <- "cubic"
n_chains <- 100
n_per_chain <- 5000
cov_scale <- 1/5
num_cores <- "ALL"


prior_mean <- c(rep(0, num_betas), mean(Y),
                rep(0, num_stationary + 1), 0)
prior_sd <- c(rep(1, num_betas), 2 * sd(Y), 
              rep(max(Y) - min(Y), num_stationary + 1), 
              0.5)


# Guess on stationary points
new_b_vec <- seq(min(X), max(X), length.out =12)
new_b_vec[2:11] <- c(100, 180, 300, 400, 550, 660, 700, 740, 790, 840)

# Will take around an hour on MacBook Air
diff_method_erp <- posterior_stat_points_diff_mcmc(X, Y, num_betas, num_stationary,
                                                   first_direction, 
                                                   zero_is_zero, interpolation, new_b_vec, 
                                                   prior_mean, prior_sd,
                                                   n_chains, n_per_chain, cov_scale,
                                                   apply_sm_correction, num_cores)

 # Density of Marginal Posterior
erp_s1 <- density(diff_method_erp$diff_stat_points[[1]], bw = "nrd")
erp_s2 <- density(diff_method_erp$diff_stat_points[[2]], bw = "nrd")
erp_s3 <- density(diff_method_erp$diff_stat_points[[3]], bw = "nrd")
erp_s4 <- density(diff_method_erp$diff_stat_points[[4]], bw = "nrd")
erp_s5 <- density(diff_method_erp$diff_stat_points[[5]], bw = "nrd")
erp_s6 <- density(diff_method_erp$diff_stat_points[[6]], bw = "nrd")
erp_s7 <- density(diff_method_erp$diff_stat_points[[7]])
erp_s8 <- density(diff_method_erp$diff_stat_points[[8]])
erp_s9 <- density(diff_method_erp$diff_stat_points[[9]])
erp_s10 <- density(diff_method_erp$diff_stat_points[[10]])

hpd_intervals <- diff_method_erp$intervals_list_diff$intervals_bonf
# Color List
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

# Plot with MAP estimate
plot(X, Y,  pch = 19, cex = 0.5, xlab = "MSEC", ylab = "ERP Signal", 
     main = "ERP Data Fit with Stationary Point HPDs", ylim = c(-11, 9), col = "gray")
lines(diff_method_erp$map_predictions[1, ], diff_method_erp$map_predictions[2, ],
      col = "coral", lwd = 3)

X_for_pred <- diff_method_erp$map_predictions[1, ]

# Add Bonferonni Intervals
for (j in 1:10){
  in_interval <- (X_for_pred >= hpd_intervals[j, 1]) & (X_for_pred <= hpd_intervals[j, 2])
  sub_X <- X_for_pred[in_interval]
  sub_predictions <- diff_method_erp$map_predictions[2, ][which(in_interval)]
  
  lines(sub_X, sub_predictions, col = col_list[j], lwd = 5)
}

# Stationary Point Posteriors
plot(erp_s1$x, erp_s1$y, type = 'l', col = col_list[1], lwd = 3, xlim = c(-100, 896), 
     ylim = c(0, 0.30), main = "Stationary Points Posteriors", xlab = "MSEC", 
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
