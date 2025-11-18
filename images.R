pkgs <- c("pracma", "mvtnorm", "splines2", "devtools")
# Install necessary packages
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

devtools::install_github('mprice747/scrdiff')
library(scrdiff)



################################################
# Figure 1, Template and Composition Visualization 
################################################
X_ex <- seq(0, 1, length.out = 1000)
diff_ex <- (tanh(3*X_ex))/tanh(3)

# Template and Composition
template_func <- -5* (X_ex - 0.5)^2 + 3
comp_func <-  -5* (diff_ex - 0.5)^2 + 3

new_stat_point <- interp1(diff_ex, X_ex, 0.5)

# Plot Template and Composition
plot(X_ex, comp_func, type = 'l', col = "#009E73", lwd = 3, 
     xlab = "x", ylab = "y", main = expression("Template and "*g(gamma(x))))
lines(X_ex, template_func, col = "#E69F00", lwd = 3)

# Plot Stationary Points
segments(new_stat_point, 0, new_stat_point, 3, 
         col = "#56B4E9", lty = "dashed" , lwd = 3)
segments(0.5, 0, 0.5, 3,  col = "#56B4E9", lty = "dashed" , lwd = 3)
points(new_stat_point, 3, col = "#009E73",  lwd = 10)
points(0.5,  3, col = "#E69F00", lwd = 10)
legend("topright", c("Template", expression(g(gamma(x)))), col = c("#E69F00", "#009E73"), 
       lwd = c(3, 3))

# Plot Diffeomorphism with Stationary Points of Template and Composition
plot(X_ex, diff_ex, type = "l", col = "#D55E00", lwd = 3, 
     xlab = "x", ylab = expression(gamma(X)),
     main = expression(gamma(x)*" Visualization") )
points(new_stat_point, 0.5, col = "#56B4E9", lwd = 5)
segments(new_stat_point, 0, new_stat_point, 0.5, 
         col = "#56B4E9", lty = "dashed" , lwd = 3)
segments(0, 0.5, new_stat_point, 0.5, 
         col = "#56B4E9", lty = "dashed" , lwd = 3)


################################################
# Figure 2, Hermite and Spline Template Strategies 
################################################

# Solve Constrained Least Squares  
# minimize ||X\beta - y||^2 s.t. A\beta = b
constrained_least_squares <- function(X, y, A, b, P_mat = 0){
  
  t_X <- t(X)
  X_t_X <- 2 *  t_X %*% X
  X_t_y <- 2 * t_X %*% y
  
  k <- nrow(A)
  
  block_matrix <- rbind(cbind(X_t_X + P_mat, t(A)), 
                        cbind(A, zeros(k, k)))
  
  solve_all <- solve(block_matrix, c(X_t_y, b))
  
  constrained_beta <- solve_all[1:ncol(X)]
  
  return(constrained_beta)
}

# Create points num_in_between evenly spaced points in between each x_point
create_points_in_between <- function(x_points, num_in_between){
  
  length_x <- length(x_points)
  
  matrix_diff <- sapply(1:(length_x - 1), function(x){return(x_points[x:(x + 1)])})
  
  matrix_ib <- apply(matrix_diff, 2, 
                     function(x){x[1] +  seq(0, num_in_between)* (x[2] - x[1])/(num_in_between + 1)})
  
  x_ib <- c(as.vector(matrix_ib), x_points[length_x])
  
  return(x_ib)
  
  
}

# Example b_vec and lambda_vec
b_vec <- c(0, 1/3, 2/3, 1)
lambda_vec <- c(0, 2, -3, 1)
# Num knots in between points in b_vec
num_knots_ib <- 2
# Data Points per segment
data_points_per_seg <- 10
# Degree of B-Spline
degree <- 3

# Create Knots
knots <- create_points_in_between(b_vec, num_knots_ib)

# Create Generate Data
x_data_points <- create_points_in_between(knots, data_points_per_seg)
y_data_points <- interp1(b_vec, lambda_vec, x_data_points, "linear")

knots_spline <- knots[2:(length(knots) - 1)]
knots_deriv <- b_vec[2:(length(b_vec)- 1)]

# Create B-spline
X_mat <- bSpline(x_data_points, knots = knots_spline, degree = degree, 
                 intercept = TRUE)

# Create Constrains (B-spline must equal lambda_vec at b_vec, and derivatives need to be 0)
A_mat <- rbind(bSpline(b_vec , knots = knots_spline, degree = degree, 
                       intercept = TRUE), 
               bSpline(b_vec , knots = knots_spline, degree = degree, 
                       intercept = TRUE, derivs = 1 )[2:(length(b_vec) -1), ])

b_vec_cons <- c(lambda_vec, rep(0, length(b_vec) -2 ))

# P-Spline Penalty
D_mat <- diff(diag(ncol(X_mat)), differences = 2)
P_mat <- 10^6 * t(D_mat) %*% D_mat

# Construct B-Spline Template Function
beta_fit <- constrained_least_squares(X_mat, y_data_points, A_mat, b_vec_cons, P_mat)

# For Visualization for Cubic Spline
test_values <- seq(0, 1, length.out = 1000)
test_mat <- bSpline(test_values, knots = knots_spline, degree = degree, 
                    intercept = TRUE)

test_y <- test_mat %*% beta_fit

# Hermite Cubic Interpolation
cubic_interp <- interp1(b_vec, lambda_vec, x_data_points, "cubic")
plot(x_data_points, cubic_interp, type = 'l', lwd = 3, xlab = "X", ylab = "Y", 
     main = "Hermite Cubic Template Function", xlim = c(-0.05, 1.05), ylim = c(-3.6, 3), 
     col = "#56B4E9")
points(b_vec, lambda_vec, pch = 19, lwd = 3)
text(b_vec[1], lambda_vec[1] - 0.4,
     bquote( ( .(bquote(b[1])) * ", " * .(bquote(lambda[1])) ) ), cex = 1.1)
text(b_vec[2], lambda_vec[2] + 0.4,
     bquote( ( .(bquote(b[2])) * ", " * .(bquote(lambda[2])) ) ), cex = 1.1)
text(b_vec[3], lambda_vec[3] - 0.4,
     bquote( ( .(bquote(b[3])) * ", " * .(bquote(lambda[3])) ) ), cex = 1.1)
text(b_vec[4], lambda_vec[4] + 0.4,
     bquote( ( .(bquote(b[4])) * ", " * .(bquote(lambda[4])) ) ), cex = 1.1)


# B-Spline Cubic Interpolation
plot(test_values, test_y, type = 'l', lwd = 3, xlab = "X", ylab = "Y", 
     main = "B-Spline Template Function", xlim = c(-0.05, 1.05), ylim = c(-3.6, 3), 
     col = "#56B4E9")
lines(x_data_points, y_data_points, col = '#E69F00', lwd = 3)
points(b_vec, lambda_vec, pch = 19, lwd = 4)
text(b_vec[1], lambda_vec[1] - 0.4,
     bquote( ( .(bquote(b[1])) * ", " * .(bquote(lambda[1])) ) ), cex = 1.1)
text(b_vec[2], lambda_vec[2] + 0.4,
     bquote( ( .(bquote(b[2])) * ", " * .(bquote(lambda[2])) ) ), cex = 1.1)
text(b_vec[3], lambda_vec[3] - 0.4,
     bquote( ( .(bquote(b[3])) * ", " * .(bquote(lambda[3])) ) ), cex = 1.1)
text(b_vec[4], lambda_vec[4] + 0.4,
     bquote( ( .(bquote(b[4])) * ", " * .(bquote(lambda[4])) ) ), cex = 1.1)
legend("topright", legend = c(expression(g[lambda]), "Generated Data"), 
       col = c("#56B4E9", "#E69F00"), lwd = c(3, 3))


################################################
# Figure 3, Diffeomorphism and DGP Method on Flat Function
################################################

# Generate Data
n <- 250
set.seed(120)
X <- runif(n, 2, 8) 
Y <- dgamma(X, 16, 4) + rnorm(n, sd = 0.05)
X <- (X - 2)/6

plot(X, Y)

X_ref <- seq(2, 8, length.out = 1000)
X_ref_norm <- (X_ref - 2)/6
true_sp <- (15/4 - 2)/6

# Diffeomorphism Parameters
num_betas <- 8
num_stationary <- 1
first_direction <- 1
zero_is_zero <- FALSE
interpolation <- "cubic"
b_vec <- NULL
n_chains <- 200
n_per_chain <- 1000
cov_scale <- 0.8
apply_sm_correction <- TRUE
num_cores <- "ALL"

prior_mean <- c(rep(0, num_betas), mean(Y),
                rep(0, num_stationary + 1), 0)
prior_sd <- c(rep(1, num_betas), 2 * sd(Y), 
              rep(max(Y) - min(Y), num_stationary + 1), 
              0.5)

# Run Inference on Flat Function
flat_func_mcmc <- posterior_stat_points_diff_mcmc(X, Y, num_betas, num_stationary,
                                first_direction, 
                                zero_is_zero, interpolation, b_vec, 
                                prior_mean, prior_sd,
                                n_chains, n_per_chain, cov_scale,
                                apply_sm_correction, num_cores)



# HPD Intervals of Stationary Points using DGP Method
flat_func_dgp <- posterior_stat_points_dgp(X, Y, num_stationary, total_resampled)

# Diffeomorphism Stationary Point Posterior
diff_dens <- density(flat_func_mcmc$diff_stat_points[[1]], bw = "nrd", 
                     n = 1000, from = 0, to =1)
# DGP Stationary Point Posterior
dgp_dens <- density(flat_func_dgp$dgp_stat_points[[1]], bw = "nrd", 
                    n = 1000, from = 0, to = 1 )

# Plot True Function and Different Posteriors
plot(X, Y, main = "Flat Function Example", xlab = "Simulated X", 
     ylab = "Simulated Y")
lines(X_ref_norm, dgamma(X_ref, 16, 4), 
      col = '#D55E00', lwd = 3)
segments(true_sp, -0.2, true_sp, dgamma(15/4, 16, 4), col = '#0072B2', lwd = 5)
legend("topright", legend = c("True Function", "True Stat. Point" ), 
       col = c("#D55E00", "#0072B2"), lwd = c(3, 3))


plot(diff_dens$x, diff_dens$y, type = 'l', col = '#009E73', lwd = 3, 
     xlab = "Stationary Point", ylab = "Posterior Density",
     main = "Stationary Point Posterior for Flat Function")

lines(dgp_dens$x, dgp_dens$y, type = 'l', col = '#E69F00', lwd = 3)
segments(true_sp, 0, true_sp, 25, col = '#0072B2', lwd = 5)
legend("topright", legend = c("Diff. Posterior", "DGP Posterior", "True Stat. Point"), 
       col = c("#009E73", "#E69F00", "#0072B2"), lwd = c(3, 3, 3))


################################################
# Figure 4, HDP Intervals 
################################################


HDP_100_df <- read.csv("Sim_Results/Sim_1/HDP_Ints/HDPs_100.csv")

diff_method_90 <- as.matrix(HDP_100_df[, c("SP1_90_Lower_Diff", "SP1_90_Upper_Diff")])
dgp_method_90 <- as.matrix(HDP_100_df[, c("SP1_90_Lower_DGP", "SP1_90_Upper_DGP")])


# Plot Confidence Intervals 
plot_conf_intervals <- function(ci_matrix, true_value, main_title,
                                col_intervals = "#444444", 
                                col_points = "black",
                                col_truth = "#0072B2", 
                                lwd_intervals = 2, 
                                lwd_truth = 5, 
                                xlim = NULL )
{
  n <- nrow(ci_matrix)
  
  # Basic plot setup
  if(is.null(xlim)){
    plot(NA, xlim = range(ci_matrix), ylim = c(0.5, n + 0.5),
         xlab = "Stationary Point", ylab = "", yaxt = "n",
         main = main_title)
  } else{
    plot(NA, xlim = xlim, ylim = c(0.5, n + 0.5),
         xlab = "Stationary Point", ylab = "", yaxt = "n",
         main = main_title)
  }
  
  
  # Draw confidence intervals
  for (i in 1:n) {
    segments(ci_matrix[i, 1], i, ci_matrix[i, 2], i, 
             col = col_intervals, lwd = lwd_intervals)
  }
  
  # Add true parameter line
  abline(v = true_value, col = col_truth, lwd = lwd_truth)
  legend("topright", legend = "True Stat. Point", col = "#0072B2", lwd = 3)
}


plot_conf_intervals(diff_method_90[1:50, ], 0.39973, 
                    main_title = "90% HDP Intervals (Diffeomorphism, n = 100)",
                    lwd_truth = 2.5, 
                    xlim = c(0, 1))
plot_conf_intervals(dgp_method_90[1:50, ], 0.39973, 
                    main_title = "90% HDP Intervals (DGP, n = 100)",
                    lwd_truth = 2.5, 
                    xlim = c(0, 1))

################################################
# Figure 5, Boxplots of Bias and Posteior Standard Deviation
################################################

Metrics_100_df <- read.csv("Sim_Results/Sim_1/Metrics/Metrics_100.csv")

diff_map <- Metrics_100_df$SP1_MAP_Diff
diff_mean <- Metrics_100_df$SP1_Post_Mean_Diff
diff_sd <- Metrics_100_df$SP1_Post_SD_Diff

dgp_map <- Metrics_100_df$SP1_MAP_DGP
dgp_mean <- Metrics_100_df$SP1_Post_Mean_DGP
dgp_sd <- Metrics_100_df$SP1_Post_SD_DGP

plot_two_boxplots <- function(data1, data2, main_title, 
                              ylab, y_line = NULL,
                              col1 = "#009E73",  # blue
                              col2 = "#E69F00",  # orange
                              line_lwd = 4) {
  
  # Combine data
  combined_data <- list(data1, data2)
  
  # Make boxplot
  boxplot(combined_data, names = c("Diffeomorphism", "DGP"),
          col = c(col1, col2), border = "black",
          ylab = ylab, 
          main = main_title,
          outline = TRUE)
  
  # Add user-specified horizontal line
  if(!is.null(y_line)) {
    abline(h = y_line, col = line_col, lwd = line_lwd, lty = 1)
    legend("topright", legend = "True Stat. Point", col = "#0072B2", lwd = 3)
  }
  
}

plot_two_boxplots(diff_mean - 0.39973, dgp_mean - 0.39973, 
                  "Posterior Mean Bias Boxplot (n = 100)", "Mean Error")
plot_two_boxplots(diff_sd, dgp_sd, 
                  "Posterior SD Boxplot (n = 100)", 
                  "Posteror SD")

