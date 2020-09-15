# This script analyses if variable selection in the case of multidemnsional
# beta and gamma vectors works
# For this reason some coefficients are 0.
library(asp20boost)

beta_vec <- c(1,0,10,0,10,0)
gamma_vec <- c(0,0.1,0,0.1,0,0.1)


# beta_2 <- c(0,1,0,1,0,1)
# gamma_2 <- c(0,0,1,0,1,0)

# beta_3 <- c(1,1,0,1,0,1)
# gamma_3 <- c(1,0,1,0,1,0)

# homoskedastic model
# beta_4 <- c(1,1,0,1,0,1)
# gamma_4 <- c(1,0,0,0,0,0)

# intercept model (heteroskedastic)
# beta_5 <- c(1,0,0,0,0,0)
# gamma_5 <- c(1,0,1,0,1,0)






# simulate data -------------------------
set.seed(1337)
n <- 500
X <- matrix(runif(n*5), ncol = 5)
y <- beta_vec[1] + beta_vec[2] * X[, 1] + beta_vec[3] * X[, 2] +  beta_vec[4] * X[, 3] + beta_vec[5] * X[, 4] + beta_vec[6] * X[, 5] + rnorm(
  n,
  sd = exp(gamma_vec[1] + gamma_vec[2] * X[, 1] + gamma_vec[3] * X[, 2] +  gamma_vec[4] * X[, 3] + gamma_vec[5] * X[, 4] + gamma_vec[6] * X[, 5])
)

# estimate data by boosting ------------
model <- LocationScaleRegressionBoost$new(
  y ~ X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5],
  ~ X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5])
gradient_boost(model, stepsize = c(0.1, 0.01), maxit = 1000, verbose = TRUE, componentwise=T)


model$beta
model$gamma

beta_vec
gamma_vec
