# This script analyses if variable selection in the case of multidemnsional
# beta and gamma vectors works
# For this reason some coefficients are 0.
library(R6)
library(asp20model)
library(asp20boost)


# in order to make this script work an extra step needs to be done:
# the error is as follow: model creation inside the "boost_true_model"-function
# below does not work, UNLESS a LocationScaleRegressionBoost-model class has been created before
# for this reason asimple model call and estimation is performed here
{ x <- runif(100)
  y <- runif(100)
  model <- LocationScaleRegressionBoost$new(y ~ x, ~ x)
  gradient_boost(model, stepsize = 0.001, maxit = 1000, abstol = 0.0001, verbose = FALSE)
}




boost_true_model <- function(beta_vec, gamma_vec){

  # simulate data -------------------------
  set.seed(1337)
  n <- 500
  X <- matrix(runif(n*5), ncol = 5)
  y <- beta_vec[1] + beta_vec[2] * X[, 1] + beta_vec[3] * X[, 2] +  beta_vec[4] * X[, 3] + beta_vec[5] * X[, 4] + beta_vec[6] * X[, 5] + rnorm(
    n,
    sd = exp(gamma_vec[1] + gamma_vec[2] * X[, 1] + gamma_vec[3] * X[, 2] +  gamma_vec[4] * X[, 3] + gamma_vec[5] * X[, 4] + gamma_vec[6] * X[, 5])
  )

  # estimate data by boosting ------------
  model <- asp20boost::LocationScaleRegressionBoost$new(
    y ~ X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5],
    ~ X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5])
  asp20boost::gradient_boost(model, stepsize = 0.001, maxit = 5000, abstol = 0.0001, verbose = FALSE)


  return(list(
    beta = model$beta,
    gamma = model$gamma)
  )
}









# beta_1 <- c(1,1,0,0,0,0)
# gamma_1 <- c(1,1,0,0,0,0)
# boost_true_model(beta_1, gamma_1)

# beta_2 <- c(0,1,0,1,0,1)
# gamma_2 <- c(0,0,1,0,1,0)
# boost_true_model(beta_2, gamma_2)

# beta_3 <- c(1,1,0,1,0,1)
# gamma_3 <- c(1,0,1,0,1,0)
# boost_true_model(beta_3, gamma_3)

# homoskedastic model
# beta_4 <- c(1,1,0,1,0,1)
# gamma_4 <- c(1,0,0,0,0,0)
# boost_true_model(beta_4, gamma_4)

# intercept model
beta_5 <- c(1,0,0,0,0,0)
gamma_5 <- c(1,0,1,0,1,0)
boost_true_model(beta_5, gamma_5)






