#' R6 Class for Specifying Location-Scale Regression Models
#'
#' @description
#' Model obejcts created by this class may be estimated by the `gradient_boost` function.
#' This model class builds upon and thus inherits the basic structure
#' from the `LocationScaleRegression` model class of the `asp20model` package.
#' Consequently, it shares similar syntax and structure. (Please  see the help
#' page of that class *Note: It would be nice to just link the vignette*)
#'
#' @usage
#' example_model <- LocationScaleRegressionBoost$new(formula, ...)
#'
#' @details
#' It assumes a normally distributed response variable with
#' one linear predictor eta_mu for the location (i.e. the mean) and one linear predictor
#' eta_sigma for scale (i.e. the standard deviation). Scale is modeled by a log-link.
#'
#'
#' @field gradients_loglik_mu Caclulates the unit-wise gradients of the loss-functions,
#' which is the negative loglikelihood.
#' @field gradients_loglik_sigma Caclulates the unit-wise gradients of the loss-functions,
#' which is the negative loglikelihood.
 #' @field update_parameters_conventional ccc
#' @field update_parameters_compwise ddd
#' @field par_log eee
#' @field plot fff
#'
# #' @inheritSection asp20model::LocationScaleRegression Section description
#'
#' @export

LocationScaleRegressionBoost <- R6Class(
  "LocationScaleRegressionBoost",
  inherit = LocationScaleRegression,

  private = list(

    #helper variables to improve computation time
    X_HAT = numeric(),
    X_HAT_componentwise = numeric(),
    X_number_of_columns = numeric(),

    Z_HAT = numeric(),
    Z_HAT_componentwise = numeric(),
    Z_number_of_columns = numeric(),


    # estimators (bj_hat) for location and scale,
    # conventional boosting --------------------------------------------------------
    gradient_estimators_mu = function() {
      private$X_HAT %*% crossprod(private$X, self$gradients_loglik_mu())
    },
    gradient_estimators_sigma = function() {
      private$Z_HAT %*% crossprod(private$Z, self$gradients_loglik_sigma())
    },
    # componentwise boosting --------------------------------------------------------
    gradient_estimators_mu_compwise = function() {
      number_of_cols <- ncol(private$X)
      result_list <- rep(NA, number_of_cols)
      for(i in 1:number_of_cols){
        result_list[i] <- private$X_HAT_componentwise[i] %*% crossprod(private$X[, i], self$gradients_loglik_mu())
       }
      return(result_list)
    },
    gradient_estimators_sigma_compwise = function() {
      number_of_cols <- ncol(private$Z)
      result_list <- rep(NA, number_of_cols)
      for(i in 1:number_of_cols){
        result_list[i] <- private$Z_HAT_componentwise[i] %*% crossprod(private$Z[, i], self$gradients_loglik_sigma())
      }
      return(result_list)
    }
),

  public = list(

    #stepsize for boosting
    stepsize_beta = numeric(),
    stepsize_gamma = numeric(),

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

    },




    #for unit testing purposes
    eta_mu = function(){
      private$X %*% self$beta
    },
    eta_sigma = function(){
      private$Z %*% self$gamma
    },

    # gradients of the log-likelihood wrt eta_mu
    # these gradients are estimated and then stepwisely added
    #' @details Returns the unit-wise gradients of the loss function, which is the negative loglikelihood.
    #' @return a numeric vector nx1
    #' @examples
    #' y <- rnorm(30)
    #' model <- LocationScaleRegressionBoost$new(y ~ 1)
    #' model$gradients_loglik_mu()
    gradients_loglik_mu = function() {
      resid <- self$resid("response")
      scale <- self$fitted_scale
      return(resid / (scale^2) )
    },

    # gradients of the log-likelihood eta_sigma.
    # these gradients are estimated and then stepwisely added
    #' @details Returns the unit-wise gradients of the loss function, which is the negative loglikelihood.
    #' @return a numeric vector nx1
    #' @examples
    #' y <- rnorm(30)
    #' model <- LocationScaleRegressionBoost$new(y ~ 1)
    #' model$gradients_loglik_sigma()
    gradients_loglik_sigma = function() {
      resid <- self$resid("response")
      scale <- self$fitted_scale
      return(resid^2 / scale^2 - 1)
    },




    # helper functions for boosting conventional <-> compwise ----------------------------------
    # a function using the gradient estimators to update parameters beta and gamma
    #' @return a numeric vector nx1
    update_parameters_conventional = function() {
      self$beta <- self$beta + self$stepsize_beta * private$gradient_estimators_mu()
      self$gamma <- self$gamma + self$stepsize_gamma * private$gradient_estimators_sigma()
    },

    # a function using the gradient estimators to update parameters beta and gamma
    #' @return a numeric vector nx1
    update_parameters_compwise = function() {

      number_of_cols_X <- dim(private$X)[2]
      number_of_cols_Z <- dim(private$Z)[2]

      losses_mu <- rep(NA, number_of_cols_X)
      losses_sigma <- rep(NA, number_of_cols_Z)

      #check if least squares is appropriate for calculating loss
      for(j in 1:number_of_cols_X) {
        losses_mu[j] <- sum((self$gradients_loglik_mu() - private$X[, j] * private$gradient_estimators_mu_compwise()[j]) ^2)
      }
      for(j in 1:number_of_cols_Z) {
        losses_sigma[j] <- sum((self$gradients_loglik_sigma() - private$Z[, j] * private$gradient_estimators_sigma_compwise()[j]) ^2)
      }

      least_loss_index_mu <- which.min(losses_mu)
      least_loss_index_sigma <- which.min(losses_sigma)

      update_vector_beta <- rep(0, number_of_cols_X)
      update_vector_beta[least_loss_index_mu] <- private$gradient_estimators_mu_compwise()[least_loss_index_mu]
      update_vector_gamma <- rep(0, number_of_cols_Z)
      update_vector_gamma[least_loss_index_sigma] <-  private$gradient_estimators_sigma_compwise()[least_loss_index_sigma]

      self$beta <- self$beta + self$stepsize_beta * update_vector_beta
      self$gamma <- self$gamma + self$stepsize_gamma * update_vector_gamma
    },


    # the following list enables plotting by logging all parameters at each iteration -------------------------------------------------
    #' enables plotting
    par_log = list(),

    #' Plots the logged model parameters for each iteration
    #' @examples
    #' y <- rnorm(30)
    #' model <- LocationScaleRegressionBoost$new(y ~ 1)
    #' gradient_boost(model, plot=TRUE)
    #' model$plot()
    plot = function() {
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


  )

)

#' Gradient Boost Algorithm for Location-Scale Regression
#'
#' This function optimizes the log-likelihood of the given location-scale
#' regression model by gradient boosting. It is designed to work specifically
#' together with the LocationScaleRegressionBoost model class.
#'
#' @param model A [LocationScaleRegressionBoost] object.
#' @param stepsize The learning rate for the parameter updates.
#' @param maxit The maximum number of iterations.
#' @param abstol The absolute convergence tolerance. The algorithm stops if the
#'               absolute value of the unit-wise gradients all drop below this value.
#' @param verbose Whether to print the progress of the algorithm.
#' @param plot Whether to plot the result of the algorithm.
#'
#' @return The updated model, invisibly.
#'
#' @examples
#' y <- rnorm(30)
#' model <- LocationScaleRegressionBoost$new(y ~ 1)
#' gradient_boost(model)
#'
#' @export

gradient_boost = function(model,
                          stepsize = c(0.01, 0.1),
                          maxit = 1000,
                          abstol = 0.001,
                          componentwise = FALSE,
                          verbose = FALSE,
                          plot = FALSE) {

  model$par_log <- list()

  #Init step size of model
  model$stepsize_beta <- stepsize[1]
  if(length(stepsize)>1)
    model$stepsize_gamma <- stepsize[2]
  else
    model$stepsize_gamma <- stepsize[1]

  for(iter in 1:maxit) {

    if(componentwise == T) model$update_parameters_compwise()
    if(componentwise == F) model$update_parameters_conventional()


    if(plot == T)
      {
      #Log parameters for plotting
      model$par_log[[iter]]<-c(model$beta, model$gamma)
      }

    if(verbose == T) {
      par_msg <- c(model$beta, model$gamma)
      par_msg <- format(par_msg, trim = TRUE, digits = 3)
      par_msg <- paste(par_msg, collapse = " ")

      loglik_msg <- format(model$loglik(), digits = 6)

      message(
        "Iteration:      ", iter, "\n",
        "Parameters:     ", par_msg, "\n",
        "Squared Resid:  ", sum(model$resid()^2), "\n",
        "Log-likelihood: ", loglik_msg, "\n",
        "==============="
      )

    }

  }
  if(plot == T)
  {
    #plot logged parameters
    model$plot()
  }

  invisible(model)

}


