set.seed(1337)


# create a model class - simplest case
n <- 500
x <- runif(n)
z <- runif(n)
y <- x + rnorm(n, sd = exp(z))
model_simple <- LocationScaleRegressionBoost$new(y ~ x, ~ z)


# create a model class - large model
n <- 500
x <- matrix(runif(5 * n), ncol = 5) # 5 location-covariates x_1 to x_10
z <- matrix(runif(5 * n), ncol = 5) # 5 scale-covariates z_1 to z_10
# formulate the true model - dependent on only some covariates
y <- 21 + 22 * x[, 1] + 22 * x[, 3] + 23 * x[, 5] + rnorm(n,
                                                          sd = 31 + 32 * z[, 2] + 33 * z[, 3] + 34 * z[, 4])
model_large <- LocationScaleRegressionBoost$new(
  y ~  x[, 1] + x[, 2] + x[, 3] + x[, 4] + x[, 5],
  ~ z[, 1] + z[, 2] + z[, 3] + z[, 4] + z[, 5])
