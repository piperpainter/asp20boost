# set up main plotting function ----------------------------------
init_model <- function(beta, gamma, n) {

  # set up linear predictor functions -------------------------------
  eta_mu <- function(x) beta[1] + beta[2] * x
  eta_sigma <- function(x) gamma[1] + gamma[2] * x

  # set up random variables -----------------------------------------
  x <- runif(n)
  y <- eta_mu(x) + rnorm(n, sd = exp(eta_sigma(x)))

  # create model object ---------------------------------------------
  LocationScaleRegressionBoost$new(y ~ x, ~ x)

}
