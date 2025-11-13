library(devtools)
devtools::install_github('mprice747/scrdiff')

library(scrdiff)

set.seed(49284)
X <- runif(100, 0, 1)
Y <- sin(1.1* 2 *X) + cos(4*X) + 0.7*2*X + rnorm(length(X), sd = 0.15)
plot(X, Y)

# Length of diffeomorphism weight vector (beta)
num_betas <- 8
# Number of stationary points
num_stationary <- 2
# Whether the function is increasing from min(X) to the first stationary point (1), or
# decreasing from min(X) to the first stationary point (-1)
first_direction <- 1
# Whether f(min(X)) = 0
zero_is_zero <- FALSE
# Interpolation type for the template
interpolation <- "cubic"
# Location of template's stationary points
b_vec <- NULL
# Number of MCMC chains
n_chains <- 100
# Number of samples per MCMC chain
n_per_chain <- 2500
# Multiplied to Covariance proposal matrices, which are the inverse
# Hessians at the posterior modes
cov_scale <- 0.63
# Whether to filter out posterior modes with low likelihood 
apply_sm_correction <- TRUE
# Number of cores for parallelization
num_cores <- "ALL"

# Prior Mean and Standard Deviation
prior_mean <- c(rep(0, num_betas), mean(Y),
                rep(0, num_stationary + 1), 0)
prior_sd <- c(rep(1, num_betas), 2 * sd(Y), 
              rep(max(Y) - min(Y), num_stationary + 1), 
              0.5)

# Run MCMC
diff_method_results <- posterior_stat_points_diff_mcmc(X, Y, num_betas, num_stationary,
                                                       first_direction, zero_is_zero, 
                                                       interpolation, b_vec, prior_mean,
                                                       prior_sd, n_chains, n_per_chain, 
                                                       cov_scale, apply_sm_correction, 
                                                       num_cores)
# Histogram of Posteior of Stationary Points
hist(diff_method_results$diff_stat_points[[1]], xlim = c(0, 1), 
     col = "#0072B2", main = "Posterior Samples of Stationary Points", 
     xlab = "X", prob = TRUE, ylim = c(0, 20))
hist(diff_method_results$diff_stat_points[[2]], add = TRUE, breaks = 20, 
     col = "#E69F00", prob = TRUE)
segments(0.2420, 0, 0.2420, 15, col = "#009E73", lwd = 5)
segments(0.6873, 0, 0.6873, 15, col = "#009E73", lwd = 5)
legend("topleft", legend = c("SP 1", "SP 2", "True SPs"), 
       col = c("#0072B2", "#E69F00", "#009E73"), pch = 19)

