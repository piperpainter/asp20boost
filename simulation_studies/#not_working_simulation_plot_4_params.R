# Since the alogrithm does not yet work properly for all possible models
# this script analyses its functionality for 2-dim beta and gamma
# also scale is regressed on the same predictor as location
library(R6)
library(asp20model)
library(asp20boost)


# in order to make this script work an extra step needs to be done:
# the error is as follow: model creation inside the "boost_true_model"-function
# below does not work, UNLESS a LocationScaleRegressionBoost-model class has been created before
# for this reason a simple empty model call and estimation is performed here
{ x <- runif(100)
  y <- y <- x + rnorm(100, sd = exp(-3 + 2 * x))
  model <- LocationScaleRegressionBoost$new(y ~ x, ~ x)
  gradient_boost(model, stepsize = 0.001, maxit = 1000, abstol = 0.0001, verbose = FALSE)
}



boost_true_model <- function(beta_0, beta_1, gamma_0, gamma_1){

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
  gradient_boost(model, stepsize = 0.001, maxit = 5000, abstol = 0.0001, verbose = FALSE)


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



  return(list(
    beta = c(model$beta[1], model$beta[2]),
    gamma = c(model$gamma[1], model$gamma[2])
  ))
}

#windows()

boost_true_model(0,1,-3,2)

#boost_true_model(0,0,0,0)
#boost_true_model(0.5,0.5,0.5,0.5)
#boost_true_model(1,1,1,1)
#boost_true_model(2,2,2,2)
#boost_true_model(-1,-1,-1,-1)
#boost_true_model(-2,-2,-2,-2)

#boost_true_model(-2,-2,1,1)
#boost_true_model(2,2,0,0.5)
