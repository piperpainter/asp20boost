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


library(R6)


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
                                          #' @details
                                          #' Least Square Residual
                                          #'
                                          #' @return
                                          #' Current gradient of the coefficents
                                          #'
                                          #' @examples
                                          #' model$<-lstSqrResid
                                          lstSqrResid = function()
                                          {
                                            self$Proj_M_X %*% self$resid()
                                          },
                                          #' @details
                                          #' First derivatives of Log Likelehood
                                          #'
                                          #' @return
                                          #' Current first derivates of each Coeffizent
                                          #'
                                          #' @examples
                                          #' model$<-derivativesZ
                                          derivativesZ = function()
                                          {
                                            #Calculates current first derivates for Gamma n of Log Likelehood
                                            #Could be further improved by moving the matrix multiplaction to the update gamma function to reduce the calculation time
                                            Z_countCoefficient <- dim(model$getZ)[2]
                                            derivatives <- numeric()
                                            for(n in 1:dim(private$Z)[2]) {
                                              derivatives[n]<-sum((self$resid()^2)*private$Z[,n]*exp(-2*(drop(private$Z %*% self$gamma)))-private$Z[,n])
                                            }
                                            derivatives
                                          },
                                          derivativesB = function()
                                          {
                                            #Calculates current first derivates for Beta n of Squared Residuals ()
                                            #Could be further improved by moving the matrix multiplaction to the update beta function to reduce the calculation time
                                            derivatives<-numeric()
                                            for(n in 1:dim(private$X)[2]) {
                                              ui <- (self$resid())
                                              xn <- private$X[,n]
                                              bjhat <- solve(t(xn)%*%xn)%*%t(xn)%*%ui
                                              derivatives[n]<-sum((ui-xn%*%bjhat)^2)
                                            }
                                            derivatives
                                          },
                                          getX = function()
                                          {
                                            private$X
                                          },
                                          getY = function()
                                          {
                                            private$y
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
  #Page 226 init b0 wiht mean of y
  tmp_DerviatesB <- c(0,0)
  tmp_DerviatesZ <- c(0,0)

  #check if location scale model
  for (i in seq_len(maxit)) {

    #Set values for checking progress -> abstol
    old_grad_beta <-grad_beta
    old_grad_gamma <-grad_gamma


    #Update Beta with greatest gradient
    #should be evaluated to check Least Squares Criterion to identify component to update?
    # indexOfBetaUpdate = which.max(abs(tmp_DerviatesB))
    # tmp_newbeta<-model$beta
    # tmp_newbeta[indexOfBetaUpdate] <- model$beta[indexOfBetaUpdate] + v*tmp_DerviatesB[indexOfBetaUpdate]
    # model$beta<-tmp_newbeta
    # tmp_DerviatesB <- model$derivativesB
    # grad_beta <- model$grad_beta()
    #END OF OLD

    #needs custom step size for Beta, otherwise it would take to long
    stepsize <- 0.1
    indexOfBetaUpdate <-which.min(model$derivativesB)

    newbeta <- model$beta
    xn <- model$getX[,indexOfBetaUpdate]
    newbeta[indexOfBetaUpdate] <- newbeta[indexOfBetaUpdate] + stepsize*solve(t(xn)%*%xn)%*% t(xn) %*% model$resid()
    model$beta <- newbeta



    #Update Gamma with greatest gradient
    #needs to be changes to fit sperate linear models for all coveraites and use min function to determin min loss
    indexOfGammaUpdate = which.max(abs(tmp_DerviatesZ))
    tmp_newgamma <-model$gamma
    tmp_newgamma[indexOfGammaUpdate] <- model$gamma[indexOfGammaUpdate] + v*tmp_DerviatesZ[indexOfGammaUpdate]
    model$gamma<-tmp_newgamma
    tmp_DerviatesZ <- model$derivativesZ
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
        "Componentwise boosting: ", loglik_msg, "\n",
        "------------Beta---------------\n",
        "Beta Update Coeff: ", which.min(tmp_DerviatesB)," ", tmp_DerviatesB[indexOfBetaUpdate], "\n",
        "------------Gamma---------------\n",
        "Gamma Update Coeff: ", which.min(tmp_DerviatesZ)," ", tmp_DerviatesZ[indexOfGammaUpdate], "\n",
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

## for testing purpose
library(asp20model)
set.seed(1337)
n <- 500
x <- runif(n)
y <- x + rnorm(n, sd = exp(-3 + 2 * x))
model <- LocationScaleRegressionBoost$new(y ~ x, ~ x)
gradient_boost(model,stepsize = 0.001, maxit = 10000, abstol = 0.0001, verbose = TRUE)

