library(R6)
library(asp20model)
library(asp20boost)

set.seed(1337)
n <- 500
x <- runif(n)
y <- x + rnorm(n, sd = exp(-3 + 2 * x))

start_time <- Sys.time()
model <- asp20boost::LocationScaleRegressionBoost$new(y ~ x, ~ x)
asp20boost::gradient_boost(model, stepsize = 0.001, maxit = 2000, abstol = 0.0001, verbose = FALSE)
end_time <- Sys.time()

# True Coefficients Location: 0, 1
# Identified Coefficients Location:
model$beta[1]
model$beta[2]

# True Coefficients Scale: -3, 2
# Identified Coefficients Scale:
model$gamma[1]
model$gamma[2]

#Time Taken
end_time - start_time
