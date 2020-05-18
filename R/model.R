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


#


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
                                            #Calculates current first derivates for Zn of Log Likelehood
                                            #Could be further improved by moving the matrix multiplaction to the update gamma function to reduce the calculation time
                                            Z_countCoefficient <- dim(model$getZ)[2]
                                            coefficients <- numeric()
                                            for(n in 1:dim(private$Z)[2]) {
                                              coefficients[n]<-sum((self$resid()^2)*private$Z[,n]*exp(-2*(drop(private$Z %*% self$gamma)))-private$Z[,n])
                                            }
                                            coefficients

                                          },

                                          lstSqrCovariates = function()
                                          {
                                            #Page 219

                                            lstSqrsCovariate <- numeric()

                                            for(j in 1:dim(private$X)[2]) {
                                              xj_mean <- mean(private$X[,j])

                                              #Compute the resdiuals for the jth covariate
                                              ui <- y - self$beta[j]*private$X[,j]

                                              #Create Prjection Matrix for jth covariate
                                              X_tr_X_inv <- solve(t(private$X[,j]) %*% private$X[,j])
                                              Proj_M_X <- X_tr_X_inv %*% t(private$X[,j])

                                              #Fit least-squares base-learning procedures for all the parameters yielding
                                              bjhat <- Proj_M_X * self$resid()

                                              lstSqrs <- sum((ui - private$X[,j]*bjhat)^2)
                                              lstSqrsCovariate[j] <- lstSqrs

                                            }
                                            #which.min(lstSqrsCovariate)
                                            lstSqrsCovariate
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

  message("boost js")
  v <- stepsize

  tmp_DerviatesB <- c(0,0)
  tmp_DerviatesZ <- c(0,0)

  #check if location scale model
  for (i in seq_len(maxit)) {

    #Compontentwise Boosting for Beta
    #
    tmp_lstSqrsCovariate <- model$lstSqrCovariates2
    indexOfCovariateUpdate = which.min(tmp_lstSqrsCovariate)

    # betanew <- model$beta
    #betanew[indexOfCovariateUpdate] <- model$beta[indexOfCovariateUpdate] + v*tmp_lstSqrsCovariate[indexOfCovariateUpdate]

    #model$beta<-betanew

    #old
    #


    model$beta <- model$beta + v*model$lstSqrResid
    grad_beta <- model$grad_beta()

    #new
    #tmp_DerviatesB <- model$lstSqrCovariates
    #indexOfBetaUpdate = which.max((tmp_DerviatesB))
    #tmp_newbeta <-model$beta
    #tmp_newbeta[indexOfBetaUpdate] <- model$beta[indexOfBetaUpdate] + v*tmp_DerviatesB[indexOfBetaUpdate]
    #model$beta<-tmp_newbeta
    #grad_beta <- model$grad_beta()

    #Update Gamma with greatest gradient

    indexOfGammaUpdate = which.max(abs(tmp_DerviatesZ))
    tmp_newgamma <-model$gamma
    tmp_newgamma[indexOfGammaUpdate] <- model$gamma[indexOfGammaUpdate] + v*tmp_DerviatesZ[indexOfGammaUpdate]
    model$gamma<-tmp_newgamma
    tmp_DerviatesZ <- model$derivativesZ
    grad_gamma <- model$grad_gamma()
    #old
    #model$gamma<-model$gamma+(v*model$derivativesZ)
    #grad_gamma <- model$grad_gamma()
    #lstsqrcov_msg <- format((abs(c( model$lstSqrCovariates))), digits = 3)
    #message("show covariats:", lstsqrcov_msg)

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
        "Log-likelihood: ", loglik_msg, "\n",
        "Componentwise boosting: ", loglik_msg, "\n",
        "------------Beta---------------\n",
        "Lst Sqr Covariate min: ", which.min(tmp_DerviatesB)," ", tmp_DerviatesB[which.min(tmp_DerviatesB)], "\n",
        "Lst Sqr Covariate: ", which.max(tmp_DerviatesB)," ", tmp_DerviatesB[which.max(tmp_DerviatesB)], "\n",
        "------------Gamma---------------\n",
        "Lst Sqr Covariate: ", which.min(tmp_DerviatesZ)," ", tmp_DerviatesZ[indexOfGammaUpdate], "\n",
        "Lst Sqr Covariate: ", which.max(tmp_DerviatesZ)," ", tmp_DerviatesZ[which.max(tmp_DerviatesZ)], "\n",
        "ABS:",all(abs(c(grad_beta, grad_gamma)) <= abstol), abs_msg, "\n",
        "==============="
      )
    }

    #early stopping needs to be checked
    if (all(abs(c(grad_beta, grad_gamma)) <= abstol))
    {
      message("abs stop at")
      message("Iteration:      ", i, "\n")
      break
    }


  }
  invisible(model)
}

## for testing purpose
#library(asp20model)
#set.seed(1337)
#
#n <- 500
#x <- runif(n)
#y <- x + rnorm(n, sd = exp(-3 + 2 * x))
#model <- LocationScaleRegressionBoost$new(y ~ x, ~ x)

#gradient_boost(model,stepsize = 0.0001, maxit = 10000, abstol = 0.001, verbose = TRUE)
