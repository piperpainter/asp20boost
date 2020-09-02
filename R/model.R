#' @title R6 class for performing either simple and generic componentwise
#'  boosting on location-scale regression models.
#'
#'
#' This model class builts upon and thus inherits the basic structure
#' from the `LocationScaleRegression` model class of the `asp20model` package.
#' Consequently, it share similar syntax and structure. (Please  see the help
#' page of that class *Note: It would be nice to just link the vignette*)
#'
#' It assumes a normally distributed response variable with
#' one linear predictor for the location (i.e. the mean) and one for the
#' scale (i.e. the standard deviation). The linear predictors for the
#' location and the scale are called Xβ and Zγ respectively. The scale
#' uses a log link.
#'
#' @field componentwiseLossGamma Calculates loss function for
#'   every component in scale-dimension `gamma`
#' @field componentwiseLossBeta Calculates loss function for
#'   every component in location-dimension `beta`
#' @field bestFittingVariableGamma Calculates the actual update for `gamma`
#' @field bestFittingVariableBeta Calculates the actual update for `beta`
#'
#' @import R6
#' @import asp20model
#' @export

LocationScaleRegressionBoost <- R6Class(
  "LocationScaleRegressionBoost",
  inherit = LocationScaleRegression,

  public = list(

    # gradients of the log-likelihood wrt eta_mu / eta_sigma.
    # these gradients are estimated and then stepwisely added
    gradients_loglik_mu = function() {
      resid <- self$resid("response")
      scale <- self$fitted_scale
      return(resid / (scale^2) )
    },
    gradients_loglik_sigma = function() {
      resid <- self$resid("response")
      scale <- self$fitted_scale
      return(resid^2 / scale^2 - 1)
    },


    # estimators (bj_hat) for location and scale,
    # conventional <-> componentwise --------------------------------------------------------

    gradient_estimators_mu = function() {
      chol2inv(chol(crossprod(private$X, private$X))) %*% crossprod(private$X, self$gradients_loglik_mu())
    },
    gradient_estimators_sigma = function() {
      chol2inv(chol(crossprod(private$Z, private$Z))) %*% crossprod(private$Z, self$gradients_loglik_sigma())
    },

    gradient_estimators_mu_compwise = function() {
      number_of_cols <- ncol(private$X)
      result_list <- rep(NA, number_of_cols)
      for(i in 1:number_of_cols){
        result_list[i] <- chol2inv(chol(crossprod(private$X[, i], private$X[, i]))) %*% crossprod(private$X[, i], self$gradients_loglik_mu())
      }
      return(result_list)
    },
    gradient_estimators_sigma_compwise = function() {
      number_of_cols <- ncol(private$Z)
      result_list <- rep(NA, number_of_cols)
      for(i in 1:number_of_cols){
        result_list[i] <- chol2inv(chol(crossprod(private$Z[, i], private$Z[, i]))) %*% crossprod(private$Z[, i], self$gradients_loglik_sigma())
      }
      return(result_list)
    },


    # helper functions for boosting conventional <-> compwise ----------------------------------
    update_parameters_conventional = function() {
      self$beta <- self$beta + 0.01 * self$gradient_estimators_mu()
      self$gamma <- self$gamma + 0.1 * self$gradient_estimators_sigma()
    }, #possible naming: update parameters conventional
    update_parameters_compwise = function() {

      number_of_cols_X <- dim(private$X)[2]
      number_of_cols_Z <- dim(private$Z)[2]

      losses_mu <- rep(NA, number_of_cols_X)
      losses_sigma <- rep(NA, number_of_cols_Z)

      #check if least squares is appropriate for calculating loss
      for(j in 1:number_of_cols_X) {
        losses_mu[j] <- sum((self$gradients_loglik_mu() - private$X[, j] * self$gradient_estimators_mu_compwise()[j]) ^2)
      }
      for(j in 1:number_of_cols_Z) {
        losses_sigma[j] <- sum((self$gradients_loglik_sigma() - private$Z[, j] * self$gradient_estimators_sigma_compwise()[j]) ^2)
      }

      least_loss_index_mu <- which.min(losses_mu)
      least_loss_index_sigma <- which.min(losses_sigma)

      update_vector_beta <- rep(0, number_of_cols_X)
      update_vector_beta[least_loss_index_mu] <- self$gradient_estimators_mu_compwise()[least_loss_index_mu]
      update_vector_gamma <- rep(0, number_of_cols_Z)
      update_vector_gamma[least_loss_index_sigma] <-  self$gradient_estimators_sigma_compwise()[least_loss_index_sigma]

      self$beta <- self$beta + 0.01 * update_vector_beta
      self$gamma <- self$gamma + 0.1 * update_vector_gamma
    },


    # the following publics enable plotting -------------------------------------------------
    par_log = list(),

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
                          stepsize = c(0.01, 0.1),
                          maxit = 1000,
                          abstol = 0.001,
                          componentwise = FALSE,
                          verbose = FALSE,
                          plot = FALSE) {

  # store gradients from last iteration step to compare them with newly calculated ones
  grad_old_mu <- model$gradients_loglik_mu()
  grad_old_sigma <- model$gradients_loglik_sigma()

  model$par_log <- list()
  for(iter in 1:maxit) {

    if(componentwise == T) model$update_parameters_compwise()
    if(componentwise == F) model$update_parameters_conventional()

    # check for substantial change in gradients of loglikelihood (wrt eta_mu, eta_sigma)
    {
    grad_change_mu <- abs(grad_old_mu - model$gradients_loglik_mu())
    grad_change_sigma <- abs(grad_old_sigma - model$gradients_loglik_sigma())
    grad_changes <- c(grad_change_mu, grad_change_sigma)
    }

    if(iter > 1 && all(grad_changes <= abstol)){
      message("early stopping at iteration: ",iter)
      break()
    }

    model$par_log[[iter]]<-c(model$beta, model$gamma)
    if(plot == T) model$plot()

    if(verbose == T) {
      par_msg <- c(model$beta, model$gamma)
      par_msg <- format(par_msg, trim = TRUE, digits = 3)
      par_msg <- paste(par_msg, collapse = " ")

      loglik_msg <- format(model$loglik(), digits = 3)

      message(
        "Iteration:      ", iter, "\n",
        "Parameters:     ", par_msg, "\n",
        "Squared Resid:  ", sum(model$resid()^2), "\n",
        "Log-likelihood: ", loglik_msg, "\n",
        "==============="
      )

    }

  }

  invisible(model)

}






