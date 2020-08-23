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

    # gradients of the log-likelihood wrt eta_mu / eta_sigma. These gradients are estimated
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


    # gradient estimators, also "bj_hat", which are added iteratively to the final estimates beta and gamma
    gradient_estimators_mu = function() {
      chol2inv(chol(crossprod(private$X, private$X))) %*% crossprod(private$X, self$gradients_loglik_mu())
    },
    gradient_estimators_sigma = function() {
      chol2inv(chol(crossprod(private$Z, private$Z))) %*% crossprod(private$Z, self$gradients_loglik_sigma())
    },


    # functions useful for caching component-wise boosting stuff
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

    # current version of algorithm needs access to full design matrices X and Z
    # impossible if private
    # --> later on cache results after checking correctness of loss calculation and get rid of public X, Z
    X_pub = function(){
      private$X
    },

    Z_pub = function(){
      private$X
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
                          componentwise = TRUE,
                          verbose = TRUE,
                          plot = TRUE) {

  # store gradients from last iteration step to compare them with newly calculated ones
  grad_old_mu <- model$gradients_loglik_mu()
  grad_old_sigma <- model$gradients_loglik_sigma()

  for(iter in 1:maxit) {

    if(componentwise == T) helper_boost_compwise(model, stepsize)
    if(componentwise == F) helper_boost_conventional(model, stepsize)

    # check for substantial change in gradients of loglikelihood (wrt eta_mu, eta_sigma)
    {
    grad_change_mu <- abs(grad_old_mu - model$gradients_loglik_mu())
    grad_change_sigma <- abs(grad_old_sigma - model$gradients_loglik_sigma())
    grad_changes <- c(grad_change_mu, grad_change_sigma)
    }

    if(iter > 1 && all(grad_changes <= abstol)){
      message("early stopping at iteration: ",i)
      break()
    }

    if (verbose == T) {
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

helper_boost_conventional <- function(model,
                                      stepsize = c(0.01, 0.1)) {

  model$beta <- model$beta + stepsize[1] * model$gradient_estimators_mu()
  model$gamma <- model$gamma + stepsize[2] * model$gradient_estimators_sigma()

  invisible(model)
}



helper_boost_compwise <- function(model,
                                  stepsize = c(0.01, 0.1)
) {

  number_of_cols_X <- dim(model$X_pub())[2]
  number_of_cols_Z <- dim(model$Z_pub())[2]

  losses_mu <- rep(NA, number_of_cols_X)
  losses_sigma <- rep(NA, number_of_cols_Z)


  #check if least squares is appropriate for calculating loss
  for(j in 1:number_of_cols_X) {
    losses_mu[j] <- sum((model$gradients_loglik_mu() - model$X_pub()[, j] * model$gradient_estimators_mu_compwise()[j]) ^2)
  }

  #check if least squares is appropriate for calculating loss
    for(j in 1:number_of_cols_Z) {
    losses_sigma[j] <- sum((model$gradients_loglik_sigma() - model$Z_pub()[, j] * model$gradient_estimators_sigma_compwise()[j]) ^2)
  }


  least_loss_index_mu <- which.min(losses_mu)
  least_loss_index_sigma <- which.min(losses_sigma)


  update_vector_beta <- rep(0, number_of_cols_X)
  update_vector_beta[least_loss_index_mu] <- model$gradient_estimators_mu_compwise()[least_loss_index_mu]

  update_vector_gamma <- rep(0, number_of_cols_X)
  update_vector_gamma[least_loss_index_sigma] <-  model$gradient_estimators_sigma_compwise()[least_loss_index_sigma]


  model$beta <- model$beta + stepsize[1] * update_vector_beta
  model$gamma <- model$gamma + stepsize[2] * update_vector_gamma

  invisible(model)
}
