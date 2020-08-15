#' @title asp20boost
#'
#' This model class [...]
#'
#' @field componentwiseLossGamma Calculates loss function for
#'   every component in scale-dimension/gamma
#' @field componentwiseLossBeta Calculates loss function for
#'   every component in location-dimension/beta
#' @field bestFittingVariableGamma Calculates the actual update for gamma
#' @field bestFittingVariableBeta Calculates the actual update for beta
#'
#' @import R6
#' @import asp20model
#' @export

LocationScaleRegressionBoost <- R6Class(
  "LocationScaleRegressionBoost",
  inherit = LocationScaleRegression,




  private = list(

    componentwise = TRUE,

    X_HAT = numeric(),
    X_HAT_componentwise = numeric(),
    X_number_of_columns = numeric(),

    Z_HAT = numeric(),
    Z_HAT_componentwise = numeric(),
    Z_number_of_columns = numeric(),

    #Extend existing update beta to calculate least square estimates for current parameters
    update_beta = function(value) {
      super$update_beta(value)
      #update least squares estimate only when model is initialized
      if(length(private$X_number_of_columns)>0){
        if(private$componentwise)
        {
          self$lstSqrEstBetaLossLoss_componentwise <- private$componentwiseLossBeta()
        }
        else
        {

          self$lstSqrEstBetaLoss <- private$X_HAT %*% crossprod(private$X, self$resid())
        }
      }
    },

    #Extend existing update gamma to calculate least square estimates for current parameters
    update_gamma = function(value) {
      super$update_gamma(value)

      #update least squares estimate only when model is initialized
      if(length(private$X_number_of_columns)>0){
        if(private$componentwise){
          self$lstSqrEstGammaLossLoss_componentwise <- private$componentwiseLossGamma()
        }
        else{
          self$lstSqrEstGammaLoss <- private$Z_HAT %*% crossprod(private$Z, self$resid("deviance"))
        }
      }
    },

    #' @details
    #' Returns the loss - squared error of the partial deviance for each covariate of Beta
    #' @return
    #' Vector of sum of squared errors for each covariate

    componentwiseLossBeta = function(){

      #Fit separate linear models for all covariates of beta
      loss <- numeric()
      resid_working <- (self$resid()) #  Our u_i for beta in page 226 step 2.

      for(n in 1:private$X_number_of_columns) {
        predictor_colwise <- private$X[,n]

        # old way to do it: estimate residuals with predictor componentwise (columnwise)
        # XX_inv_colwise <- t(predictor_colwise) %*% predictor_colwise)
        # bjhat <- solve(XX_inv_colwise %*% t(predictor_colwise) %*% resid_working

        # new way including cholesky calculations
        bjhat <- private$X_HAT_componentwise[n] %*% crossprod(predictor_colwise, resid_working)

        loss[n] <- sum((resid_working - predictor_colwise %*% bjhat) ^ 2)
      }

      return(loss)

    },



    #' @details
    #' Returns the loss - squared error of the partial deviance for each covariate of Gamma
    #' @return
    #' Vector of sum of squared errors for each covariate

    componentwiseLossGamma = function(){

      #Fit separate linear models for all covariates of gamma
      loss <- numeric()
      resid_deviance <- self$resid("deviance") # Our u_i for gamma in page 226 step 2
      resid_working <- self$resid()
      number_of_columns <- dim(private$Z)[2]

      for(n in 1:number_of_columns) {
        predictor_colwise <- private$Z[,n]

        # u -> first derivation of Log Likelihood - see documentation for further reading
        grad_loss_function <- ((resid_working ^ 2) * predictor_colwise * exp(-2 * (drop(private$Z %*% self$gamma))) - predictor_colwise)
        # for location-estimation the gradient of the loss function results in the residuals
        # for scale-estimation another term results

        #Kneib book Page 226

        #Can be further improvied by moving ProjectionMatrix Calculation to initializalisation

        # old way to calculate "residual-estimator" bjhat
        # ZZ_inv_colwise <- solve(t(predictor_colwise) %*% predictor_colwise)
        # bjhat <- ZZ_inv_colwise %*% t(predictor_colwise) %*% grad_loss_function

        # new way including cholesky calculations
        bjhat <- private$Z_HAT_componentwise[n] %*% crossprod(predictor_colwise, grad_loss_function)


        # Calculate new squared error by subtracting the partial deviance of covariate with new gamma (first derivation u) of the total deviance
        partial_deviance <- (resid_working / predictor_colwise %*% bjhat)
        loss[n] <- sum((resid_deviance - partial_deviance) ^ 2) # Punkt 3 Formel
      }

      return(loss)
    }

  ),

  public = list(
    lstSqrEstBetaLoss = numeric(),
    lstSqrEstBetaLossLoss_componentwise = numeric(),
    lstSqrEstGammaLoss = numeric(),
    lstSqrEstGammaLossLoss_componentwise = numeric(),

    par_log = list(),

    initialize = function(location,
                          scale = ~1,
                          data = environment(location),
                          ...) {
      #first init super class LocationScaleRegression
      super$initialize(location, scale,data)

      private$X_number_of_columns <- dim(private$X)[2]
      private$Z_number_of_columns <- dim(private$Z)[2]


      private$X_HAT <- chol2inv(chol(crossprod(private$X, private$X)))
      private$Z_HAT <- chol2inv(chol(crossprod(private$Z, private$Z)))


      #init componentwise Hat Matrix to improve computation time for X
      for(n in 1:private$X_number_of_columns) {
        predictor_colwise <- private$X[,n]
        private$X_HAT_componentwise[n] <- chol2inv(chol(crossprod(predictor_colwise, predictor_colwise)))
      }
      #init componentwise Hat Matrix to improve computation time for Z
      for(n in 1:private$Z_number_of_columns) {
        predictor_colwise <- private$Z[,n]
        private$Z_HAT_componentwise[n] <- chol2inv(chol(crossprod(predictor_colwise, predictor_colwise)))
      }

      #Init least square estimates for components of beta and gamma
      self$lstSqrEstBetaLossLoss_componentwise <- private$componentwiseLossBeta()
      self$lstSqrEstGammaLossLoss_componentwise <- private$componentwiseLossGamma()
    },

    #enable or disable componentwise boosting
    setComponentwiseBoosting = function(value)
    {
      private$componentwise = value
    },



    #' @details
    #' Determines the best variable for componentwise boosting
    #'
    #' @return
    #' New Beta Vector with updated values on the best fitted covariate
    #' @export

    coefficentUpdateBeta = function()
    {
      #build update vector
      updateBeta <- replicate(length(model$beta), 0)
      if(private$componentwise){
      #determine index of the best-fitting variable
      indexOfBetaUpdate = which.min(model$lstSqrEstBetaLossLoss_componentwise)



      #Calculates current first derivates for Beta - the residuals
      #updateBeta[indexOfBetaUpdate] <- solve(t(private$X[,indexOfBetaUpdate])%*%private$X[,indexOfBetaUpdate])%*% t(private$X[,indexOfBetaUpdate]) %*% model$resid()
      updateBeta[indexOfBetaUpdate] <- private$X_HAT_componentwise[indexOfBetaUpdate] %*% crossprod(private$X[,indexOfBetaUpdate], model$resid())
      }
      else{

        if(length(self$lstSqrEstBetaLoss)>0)
        {
          updateBeta <- updateBeta+self$lstSqrEstBetaLoss
        }


      }
      updateBeta

    },



    #' @details
    #' Determines the best variable for componentwise boosting
    #' @return
    #' New Gamma Vector with updated values on the best fitted covariate
    #' @export

    coefficentUpdateGamma = function()
    {
      #build update vector
      updateGamma <- replicate(length(model$gamma), 0)
      if(private$componentwise){
      #determine index of the best-fitting variable, the covariate with the greates influance on the deviance will decrease the loss at most
      indexOfGammaUpdate = which.min(model$lstSqrEstGammaLossLoss_componentwise)

      #Calculates current first derivates for Gamma (index of the best-fitting variable) again like u, but as sum, since Gamma covariate is a single
      updateGamma[indexOfGammaUpdate] <- sum((self$resid()^2)*private$Z[,indexOfGammaUpdate]*exp(-2*(drop(private$Z %*% self$gamma)))-private$Z[,indexOfGammaUpdate])
      }
      else{
        if(length(self$lstSqrEstGammaLoss)>0)
        {
          updateGamma<- updateGamma+self$lstSqrEstGammaLoss
        }
      }


      updateGamma
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
                          abstol = 0.001, componentwise = TRUE,
                          verbose = TRUE, plot = TRUE) {

  model$setComponentwiseBoosting(componentwise)
  grad_beta <- model$grad_beta()
  grad_gamma <- model$grad_gamma()
  v <- stepsize
  #Page 226 init b0 wiht mean of y?
  tmp_DerviatesB <- c(0,0)
  tmp_DerviatesZ <- c(0,0)


  model$par_log <- list()

  #check if location scale model
  for (i in seq_len(maxit)) {

    #Set values for checking progress -> abstol
    old_grad_beta <-grad_beta
    old_grad_gamma <-grad_gamma

    #needs custom step size for Beta, otherwise it would take to long ? increase with this factor
    stepsizemultiplierbeta <-10

    model$beta<-model$beta + stepsizemultiplierbeta*v*model$coefficentUpdateBeta()
    model$gamma<-model$gamma + v*model$coefficentUpdateGamma()


    grad_beta <- model$grad_beta()
    grad_gamma <- model$grad_gamma()


    model$par_log[[i]]<-c(model$beta, model$gamma)
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


    if (i>1&& all(abs(c(grad_beta-old_grad_beta, grad_gamma-old_grad_gamma)) <= abstol)) {
      message("early stopping at iteration:",i)
      #break()
    }


  }

  }

  if(plot)
  {
    df <- data.frame(matrix(unlist(model$par_log), nrow=length(model$par_log), byrow=T))

    plot(df$X1, type = "l",
         ylim = c(min(df),max(df)),
         xaxt = "n",
         xlab = "No. of iterations",
         ylab = "Parameter")
    axis(1, at = c(1:dim(df)[1]), c(1:dim(df)[1]))
    for(n in 2:dim(df)[2])
    {
      points(df[n], type = "l", col = "black")
    }
  }

  invisible(model)
}

# library(R6)
# library(asp20model)
# start_time <- Sys.time()
# library(asp20model)
# set.seed(1337)
# n <- 500
# x <- runif(n)
# y <- x + rnorm(n, sd = exp(-3 + 2 * x))
# model <- LocationScaleRegressionBoost$new(y ~ x, ~ x)
# gradient_boost(model,stepsize = 0.001, maxit = 10000, abstol = 0.0001,componentwise = TRUE, verbose = TRUE, plot=TRUE)
# end_time <- Sys.time()
# end_time - start_time



