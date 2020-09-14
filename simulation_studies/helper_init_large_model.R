init_large_model <- function(beta_range = c(10, 20),
                             gamma_range = c(0.1, 0.3)
){


  # sample coefficients ------------------------------------------------------------
  beta <- round(runif(11, min = 10, max = 20), digits = 0)
  gamma <- round(runif(11, min = 0.1, max = 0.3), digits = 1)

  beta[seq(3,11,3)] <- 0 # set every third beta coef to 0
  gamma[seq(2,10,2)] <- 0 # set every second gamma coef to 0
  beta[1] <- 0 # set intercept to 0 for convenience
  gamma[1] <- -3 # set gamma intercept to keep sds in an estimable range
  #beta[2] <- 280 # set one very high effect for checking]


  # sample covariates --------------------------------------------------------------
  # we sample 100 covariates. intercept is automatically added by model class
  covariates_X <- runif(100 * 11) # generate 100 obs each for 100 covariates
  covariates_Z <- runif(100 * 11)

  X <- matrix(covariates_X, nrow = 100, byrow = TRUE)
  X[, 1] <- 1
  Z <- matrix(covariates_Z, nrow = 100, byrow = TRUE)
  Z[, 1] <- 1


  # sample response ----------------------------------------------------------------
  y <- rep(NA, times = 100)
  eta_mu <- rep(NA, 100)
  eta_sigma <- rep(NA, 100)

  for(i in 1:100) {
    x_i <- X[i, ]
    eta_mu[i] <- x_i %*% beta
    eta_sigma[i] <- x_i %*% gamma
    y[i] <- rnorm(n = 1,
                  mean = x_i %*% beta,
                  sd = exp(x_i %*% gamma)
                  #sd = exp(z_i) %*% gamma
                  #sd = 0.01
    )
  }


  # create model call --------------------------------------------------------------
  xnam <- paste("x", 1:10, sep="")
  for(i in 1:10) assign(xnam[i], X[, i+1]) # i+1 since first term is no covariate, but intercept-1-vector
  fmla_loc <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  fmla_scale <- as.formula(paste("~ ", paste(xnam, collapse= "+")))
  mod_1 <- LocationScaleRegressionBoost$new(fmla_loc, fmla_scale)


  # return model -------------------------------------------------------------------
  result_list <- list()
  result_list$model <- mod_1
  result_list$design <- X
  result_list$response <- y
  result_list$fitted_response <- eta_mu
  result_list$fitted_sd <- exp(eta_sigma)
  return(result_list)


}
