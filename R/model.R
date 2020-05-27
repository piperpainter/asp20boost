#' R6 class for location-scale regression models
#'
#' This model class assumes a normally distributed response variable with
#' one linear predictor for the location (i.e. the mean) and one for the scale
#' (i.e. the standard deviation). The linear predictors for the location and
#' the scale are called \eqn{X\beta} and \eqn{Z\gamma} respectively. The scale
#' uses a log link.
#'
#' @field beta A numeric vector with the `beta` parameters.
#' @field gamma A numeric vector with the `gamma` parameters.
#' @field fitted_location A numeric vector with the fitted values
#'                        for the location.
#' @field fitted_scale A numeric vector with the fitted values
#'                     for the scale.
#'
#' @import R6
#' @import asp20model
#' @export

#Uncomment for testing
# library(R6)

LocationScaleRegressionBoost <- R6Class("LocationScaleRegression",
                                        inherit = LocationScaleRegression,
                                        public = list(

                                          X_tr_X_inv = numeric(),
                                          Proj_M_X = numeric(),
                                          Z_tr_Z_inv = numeric(),
                                          Proj_M_Z = numeric(),

                                          initialize = function(location,
                                                                scale = ~1,
                                                                data = environment(location),
                                                                ...)
                                          {
                                            super$initialize(location, scale,data,...)

                                            #𝑃𝑟𝑜𝑗𝑀𝑎𝑡𝑟𝑖𝑥=(𝑋𝑇𝑋)−1𝑋𝑇
                                            #https://stats.stackexchange.com/questions/154485/least-squares-regression-step-by-step-linear-algebra-computation?noredirect=1&lq=1
                                            #Converts X and Z to Projection Matrix for further calculations
                                            self$X_tr_X_inv <- solve(t(private$X) %*% private$X)
                                            self$Proj_M_X <- self$X_tr_X_inv %*% t(private$X)

                                            self$Z_tr_Z_inv <- solve(t(private$Z) %*% private$Z)
                                            self$Proj_M_Z <- self$Z_tr_Z_inv %*% t(private$Z)

                                          }
                                        ),
                                        active = list(
                                             bestFittingVariableGamma = function()
                                            {
                                            #Fit separate linear models for all covariates of gamma
                                            loss <- numeric()
                                            ui <- self$resid("deviance")
                                            for(n in 1:dim(private$Z)[2]) {
                                               zn <- private$Z[,n]

                                               u <- ((self$resid()^2)*private$Z[,n]*exp(-2*(drop(private$Z %*% self$gamma)))-private$Z[,n])

                                               bjhat <- solve(t(private$Z[,n])%*%private$Z[,n])%*%t(private$Z[,n])%*% u
                                               loss[n]<-sum((ui-(self$resid()/zn%*%bjhat))^2)
                                               }
                                            #determine index of the best-fitting variable
                                            indexOfGammaUpdate = which.min(loss)
                                            #build update vector
                                            updateGamma <- replicate(length(model$gamma), 0)
                                            #Calculates current first derivates for Gamma (index of the best-fitting variable) of Log Likelehood
                                            updateGamma[indexOfGammaUpdate] <- sum((self$resid()^2)*private$Z[,indexOfGammaUpdate]*exp(-2*(drop(private$Z %*% self$gamma)))-private$Z[,indexOfGammaUpdate])
                                            updateGamma
                                          },
                                          bestFittingVariableBeta = function()
                                          {
                                            #Fit separate linear models for all covariates of beta
                                            loss <- numeric()
                                            ui <- (self$resid())
                                            for(n in 1:dim(private$X)[2]) {
                                              xn <- private$X[,n]
                                              bjhat <- solve(t(private$X[,n])%*%private$X[,n])%*%t(private$X[,n])%*%ui
                                              loss[n]<-sum((ui-xn%*%bjhat)^2)
                                            }
                                            #determine index of the best-fitting variable
                                            indexOfBetaUpdate = which.min(loss)

                                            #build update vector
                                            updateBeta <- replicate(length(model$beta), 0)
                                            #Calculates current first derivates for Beta - the residuals
                                            updateBeta[indexOfBetaUpdate] <- solve(t(private$X[,indexOfBetaUpdate])%*%private$X[,indexOfBetaUpdate])%*% t(private$X[,indexOfBetaUpdate]) %*% model$resid()
                                            updateBeta

                                          }
                                        )


)

#' Gradient bost for the `LocationScaleRegressionBoost` model class
#'
#' This function optimizes the log-likelihood of the given location-scale
#' regression model by gradient descent. It has a side effect on the `model`
#' object.
#'
#' @param model A [`LocationScaleRegressionBoost`] object.
#' @param stepsize The scaling factor for the gradient.
#' @param maxit The maximum number of iterations.
#' @param abstol The absolute convergence tolerance. The algorithm stops if the
#'               absolute value of the gradient drops below this value.
#' @param verbose Whether to print the progress of the algorithm.
#'
#' @return
#' The updated model, invisibly.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionBoost$new(y ~ 1)
#' gradient_boost(model)
#'
#' @export



gradient_boost = function(model,
                          stepsize = 0.001,
                          maxit = 1000,
                          abstol = 0.001,
                          verbose = TRUE) {
  grad_beta <- model$grad_beta()
  grad_gamma <- model$grad_gamma()
  v <- stepsize
  #Page 226 init b0 wiht mean of y?
  tmp_DerviatesB <- c(0,0)
  tmp_DerviatesZ <- c(0,0)

  #check if location scale model
  for (i in seq_len(maxit)) {

    #Set values for checking progress -> abstol
    old_grad_beta <-grad_beta
    old_grad_gamma <-grad_gamma

    #needs custom step size for Beta, otherwise it would take to long ? increase with this factor
    stepsize <-10

    model$beta<-model$beta + stepsize*v*model$bestFittingVariableBeta
    grad_beta <- model$grad_beta()


    model$gamma<-model$gamma + v*model$bestFittingVariableGamma
    grad_gamma <- model$grad_gamma()

    if (verbose) {
      par_msg <- c(model$beta, model$gamma)
      par_msg <- format(par_msg, trim = TRUE, digits = 3)
      par_msg <- paste(par_msg, collapse = " ")

      grad_msg <- c(grad_beta, grad_gamma)
      grad_msg <- format(grad_msg, trim = TRUE, digits = 3)
      grad_msg <- paste(grad_msg, collapse = " ")

      loglik_msg <- format(model$loglik(), digits = 3)

      abs_msg <- format((abs(c(grad_beta, grad_gamma))), digits = 3)

      message(
        "Iteration:      ", i, "\n",
        "Parameters:     ", par_msg, "\n",
        "Gradient:       ", grad_msg, "\n",
        "Squared Resid:       ", sum(model$resid()^2), "\n",
        "Log-likelihood: ", loglik_msg, "\n",
        "------------Beta---------------\n",
        #"Beta Update Coeff: ", which.min(tmp_DerviatesB)," ", tmp_DerviatesB[indexOfBetaUpdate], "\n",
        "------------Gamma---------------\n",
        #"Gamma Update Coeff: ", which.min(tmp_DerviatesZ)," ", tmp_DerviatesZ[indexOfGammaUpdate], "\n",
        #"ABS:",all(abs(c(grad_beta-old_grad_beta, grad_gamma-old_grad_gamma)))  , "\n",
        "==============="
      )
    }


    if (i>1&& all(abs(c(grad_beta-old_grad_beta, grad_gamma-old_grad_gamma)) <= abstol)) {
       message("break")
       break()
    }


  }
  invisible(model)
}

## for testing p for testing purpose
# library(asp20model)
# set.seed(1337)
# n <- 500
# x <- runif(n)
# y <- x + rnorm(n, sd = exp(-3 + 2 * x))
# model <- LocationScaleRegressionBoost$new(y ~ x, ~ x)
# model,stepsize = 0.001, maxit = 20000, abstol = 0.0001, verbose = TRUE)
