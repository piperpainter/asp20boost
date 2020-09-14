library(asp20boost)

# set true model parameters
beta_0  <- 2
beta_1  <- 2
gamma_0 <- 2
gamma_1 <- 2



# simulate data -------------------------
set.seed(1337)
n <- 500
x <- runif(n)
y <- beta_0 + beta_1 * x + rnorm(
  n,
  sd = exp(gamma_0 + gamma_1 * x)
)

# estimate data by boosting ------------
model <- LocationScaleRegressionBoost$new(y ~ x, ~ x)
gradient_boost(model, stepsize = 0.001, maxit = 5000, abstol = 0.00000000000000000000000000000000000000000001, verbose = TRUE)


# plot results -------------------------
x_grid <- seq(0,1, length.out = 100)
y_grid <- beta_0 + beta_1 * x_grid

# limits for plot
ylim_upper <- model$beta[1] + model$beta[2] * 1 + 1.96 * exp(model$gamma[1] + model$gamma[2] * 1)
ylim_lower <- model$beta[1] + model$beta[2] * 1 - 1.96 * exp(model$gamma[1] + model$gamma[2] * 1)


# plot data points
plot(x, y, col = 2,
     xlim = c(-0.3,1.3),
     ylim = c(ylim_upper * 1.05, ylim_lower * 1.05))

# plot true model
lines(x_grid, y_grid)
lines(x_grid, y_grid + 1.96 * exp(gamma_0 + gamma_1 * x_grid))
lines(x_grid, y_grid - 1.96 * exp(gamma_0 + gamma_1 * x_grid))

#plot estimated effect
location_estimate <- model$beta[1] + model$beta[2] * x_grid
scale_estimate_sd <- exp(model$gamma[1] + model$gamma[2] * x_grid)
lines(x_grid, location_estimate, col = 3)
lines(x_grid, location_estimate + 1.96 * scale_estimate_sd, col = 3)
lines(x_grid, location_estimate - 1.96 * scale_estimate_sd, col = 3)

model$beta
model$gamma
