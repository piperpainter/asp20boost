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

    grad_loglik_mu = function() {
      resid <- self$resid("response")
      scale <- self$fitted_scale
      return(resid / (scale^2) )
    },

    grad_loglik_sigma = function() {
      resid <- self$resid("response")
      scale <- self$fitted_scale
      return(resid^2 / scale^2 - 1)
    },

    X_pub = function() {
      private$X
    },

    Z_pub = function() {
      private$Z
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
  grad_old_mu <- model$grad_loglik_mu()
  grad_old_sigma <- model$grad_loglik_sigma()

  for(iter in 1:maxit) {

    if(componentwise == T) helper_boost_compwise(model, stepsize)
    if(componentwise == F) helper_boost_conventional(model, stepsize)

    # check for substantial change in gradients
    grad_change_mu <- abs(grad_old_mu - model$grad_loglik_mu())
    grad_change_sigma <- abs(grad_old_sigma - model$grad_loglik_sigma())
    grad_changes <- c(grad_change_mu, grad_change_sigma)

    if(iter > 1 && all(grad_changes) <= abstol){
      message("early stopping at iteration: ",i)
      break()
    }

  #   if (verbose == T) {
  #     par_msg <- c(model$beta, model$gamma)
  #     par_msg <- format(par_msg, trim = TRUE, digits = 3)
  #     par_msg <- paste(par_msg, collapse = " ")
  #
  #     grad_msg <- c(grad_beta, grad_gamma)
  #     grad_msg <- format(grad_msg, trim = TRUE, digits = 3)
  #     grad_msg <- paste(grad_msg, collapse = " ")
  #
  #     loglik_msg <- format(model$loglik(), digits = 3)
  #
  #     abs_msg <- format((abs(c(grad_beta, grad_gamma))), digits = 3)
  #
  #     message(
  #       "Iteration:      ", i, "\n",
  #       "Parameters:     ", par_msg, "\n",
  #       "Gradient:       ", grad_msg, "\n",
  #       "Squared Resid:       ", sum(model$resid()^2), "\n",
  #       "Log-likelihood: ", loglik_msg, "\n",
  #       "------------Beta---------------\n",
  #       #"Beta Update Coeff: ", which.min(tmp_DerviatesB)," ", tmp_DerviatesB[indexOfBetaUpdate], "\n",
  #       "------------Gamma---------------\n",
  #       #"Gamma Update Coeff: ", which.min(tmp_DerviatesZ)," ", tmp_DerviatesZ[indexOfGammaUpdate], "\n",
  #       #"ABS:",all(abs(c(grad_beta-old_grad_beta, grad_gamma-old_grad_gamma)))  , "\n",
  #       "==============="
  #     )
  #
  # }
  #
  # }
  # if(plot == T) {
  #   df <- data.frame(matrix(unlist(model$par_log), nrow=length(model$par_log), byrow=T))
  #   plot(df$X1, type = "l",
  #        ylim = c(min(df),max(df)),
  #        xaxt = "n",
  #        xlab = "No. of iterations",
  #        ylab = "Parameter")
  #   axis(1, at = c(1:dim(df)[1]), c(1:dim(df)[1]))
  #   for(n in 2:dim(df)[2])
  #   {
  #     points(df[n], type = "l", col = "black")
  #   }
  # }


  invisible(model)

}

helper_boost_conventional <- function(model,
                                      stepsize = c(0.01, 0.1)) {


  X <- model$X_pub()
  Z <- model$Z_pub()

  number_of_cols_X <- dim(X)[2]
  number_of_cols_Z <- dim(Z)[2]

  gradients_mu <- model$grad_loglik_mu()
  gradients_sigma <- model$grad_loglik_sigma()

  gradient_estimators_mu <-  solve(t(X) %*% X) %*% t(X) %*% gradients_mu
  gradient_estimators_sigma <-  solve(t(Z) %*% Z) %*% t(Z) %*% gradients_sigma

  model$beta <- model$beta + stepsize[1] * gradient_estimators_mu
  model$gamma <- model$gamma + stepsize[2] * gradient_estimators_sigma

  invisible(model)
}



helper_boost_compwise <- function(model,
                         stepsize = c(0.01, 0.1)
                         ) {


  X <- model$X_pub()
  Z <- model$Z_pub()

  number_of_cols_X <- dim(X)[2]
  number_of_cols_Z <- dim(Z)[2]

  gradients_mu <- model$grad_loglik_mu()
  gradients_sigma <- model$grad_loglik_sigma()


  gradient_estimators_mu <- rep(NA, number_of_cols_X)
  gradient_estimators_sigma <- rep(NA, number_of_cols_Z)
  losses_mu <- rep(NA, number_of_cols_X)
  losses_sigma <- rep(NA, number_of_cols_Z)


  for(j in 1:number_of_cols_X) {
    gradient_estimators_mu[j] <-  solve(t(X[, j]) %*% X[, j]) %*% t(X[, j]) %*% gradients_mu
    losses_mu[j] <- sum((gradients_mu - X[, j] * gradient_estimators_mu[j]) ^2)
  }

  for(j in 1:number_of_cols_Z) {
    gradient_estimators_sigma[j] <-  solve(t(Z[, j]) %*% Z[, j]) %*% t(Z[, j]) %*% gradients_sigma
    losses_sigma[j] <- sum((gradients_sigma - Z[, j] * gradient_estimators_sigma[j]) ^2)
  }


  least_loss_index_mu <- which.min(losses_mu)
  least_loss_index_sigma <- which.min(losses_sigma)


  update_vector_beta <- rep(0, number_of_cols_X)
  update_vector_beta[least_loss_index_mu] <- gradient_estimators_mu[least_loss_index_mu]

  update_vector_gamma <- rep(0, number_of_cols_X)
  update_vector_gamma[least_loss_index_sigma] <- gradient_estimators_sigma[least_loss_index_sigma]


  model$beta <- model$beta + stepsize[1] * update_vector_beta
  model$gamma <- model$gamma + stepsize[2] * update_vector_gamma

  invisible(model)
}
