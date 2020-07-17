library(R6)
library(asp20model)
library(asp20boost)

set.seed(1337)
n <- 500
x <- matrix(runif(10 * n), ncol = 10) #10 covariates x_1 to x_10

y <- 0 + x[, 1] + x[, 3] + x[, 5] + rnorm(n, sd = 1)
# only half of them influence response

start_time <- Sys.time()
model <- asp20boost::LocationScaleRegressionBoost$new(
  y ~ x[, 1] + x[, 2] + x[, 3] + x[, 4] + + x[, 5],
  ~ 1 #homogenous scale
  )
asp20boost::gradient_boost(model, stepsize = 0.0001, maxit = 5000, abstol = 0.0001, verbose = FALSE)
end_time <- Sys.time()


# True Coefficients Location: 0, 1, 0, 1, 0, 1
# Identified Coefficients Location:
model$beta

#Time Taken
end_time - start_time
