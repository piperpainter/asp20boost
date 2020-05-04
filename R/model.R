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
#'
#' @importFrom R6 R6Class
#' @export

LocationScaleRegression <- R6Class(
  classname = "LocationScaleRegression",
  public = list(
    #' @details
    #' Create a new `LocationScaleRegression` object.
    #'
    #' @param location A two-sided formula with the response variable on the
    #'                 LHS and the predictor for the location (i.e. the mean)
    #'                 on the RHS.
    #' @param scale A one-sided formula with the predictor for the scale
    #'              (i.e. the standard deviation) on the RHS.
    #' @param data A data frame (or list or environment) in which to evaluate
    #'             the `location` and `scale` formulas.
    #' @param ... Passed on to [stats::model.matrix()].
    #'
    #' @return
    #' A `LocationScaleRegression` object.
    #'
    #' @examples
    #' y <- rnorm(30)
    #' LocationScaleRegression$new(y ~ 1)
    #'
    #' @importFrom stats model.matrix

    initialize = function(location,
                          scale = ~1,
                          data = environment(location),
                          ...) {
      scale <- update(scale, paste(location[[2]], "~ ."))
      private$y <- eval(location[[2]], data, environment(location))
      private$X <- model.matrix(location, data, ...)
      private$Z <- model.matrix(scale, data, ...)

      self$beta <- rep.int(0, ncol(private$X))
      self$gamma <- rep.int(0, ncol(private$Z))

      invisible(self)
    },

    #' @details
    #' Returns the log-likelihood of a `LocationScaleRegression` object
    #' at the current parameter values.
    #'
    #' @return
    #' A single number.
    #'
    #' @examples
    #' y <- rnorm(30)
    #' model <- LocationScaleRegression$new(y ~ 1)
    #' model$loglik()
    #'
    #' @importFrom stats dnorm

    loglik = function() {
      location <- private$fitted$location
      scale <- private$fitted$scale

      sum(dnorm(private$y, location, scale, log = TRUE))
    },

    #' @details
    #' Returns the gradient of the log-likelihood of a
    #' `LocationScaleRegression` object with respect to \eqn{\beta}
    #' at the current parameter values.
    #'
    #' @return
    #' A numeric vector.
    #'
    #' @examples
    #' y <- rnorm(30)
    #' model <- LocationScaleRegression$new(y ~ 1)
    #' model$grad_beta()

    grad_beta = function() {
      location <- private$fitted$location
      scale <- private$fitted$scale
      resid <- private$resid

      drop((resid / scale^2) %*% private$X)
    },

    #' @details
    #' Returns the gradient of the log-likelihood of a
    #' `LocationScaleRegression` object with respect to \eqn{\gamma}
    #' at the current parameter values.
    #'
    #' @return
    #' A numeric vector.
    #'
    #' @examples
    #' y <- rnorm(30)
    #' model <- LocationScaleRegression$new(y ~ 1)
    #' model$grad_gamma()

    grad_gamma = function() {
      location <- private$fitted$location
      scale <- private$fitted$scale
      resid <- private$resid

      drop(((resid / scale)^2 - 1) %*% private$Z)
    }
  ),
  private = list(
    y = numeric(),
    X = numeric(),
    Z = numeric(),
    .beta = numeric(),
    .gamma = numeric(),
    fitted = list(location = numeric(), scale = numeric()),
    resid = numeric(),
    update_beta = function(value) {
      private$.beta <- value
      private$fitted$location <- drop(private$X %*% self$beta)
      private$resid <- private$y - private$fitted$location
      invisible(self)
    },
    update_gamma = function(value) {
      private$.gamma <- value
      private$fitted$scale <- exp(drop(private$Z %*% self$gamma))
      invisible(self)
    }
  ),
  active = list(
    beta = function(value) {
      if (missing(value)) private$.beta else private$update_beta(value)
    },
    gamma = function(value) {
      if (missing(value)) private$.gamma else private$update_gamma(value)
    }
  )
)


#' Gradient descent for the `LocationScaleRegression` model class
#'
#' This function optimizes the log-likelihood of the given location-scale
#' regression model by gradient descent. It has a side effect on the `model`
#' object.
#'
#' @param model A [`LocationScaleRegression`] object.
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
#' model <- LocationScaleRegression$new(y ~ 1)
#' gradient_descent(model)
#'
#' @export

gradient_descent <- function(model,
                             stepsize = 0.001,
                             maxit = 1000,
                             abstol = 0.001,
                             verbose = FALSE) {
  grad_beta <- model$grad_beta()
  grad_gamma <- model$grad_gamma()

  for (i in seq_len(maxit)) {
    model$beta <- model$beta + stepsize * grad_beta
    model$gamma <- model$gamma + stepsize * grad_gamma

    grad_beta <- model$grad_beta()
    grad_gamma <- model$grad_gamma()

    if (verbose) {
      par_msg <- c(model$beta, model$gamma)
      par_msg <- format(par_msg, trim = TRUE, digits = 3)
      par_msg <- paste(par_msg, collapse = " ")

      grad_msg <- c(grad_beta, grad_gamma)
      grad_msg <- format(grad_msg, trim = TRUE, digits = 3)
      grad_msg <- paste(grad_msg, collapse = " ")

      loglik_msg <- format(model$loglik(), digits = 3)

      message(
        "Iteration:      ", i, "\n",
        "Parameters:     ", par_msg, "\n",
        "Gradient:       ", grad_msg, "\n",
        "Log-likelihood: ", loglik_msg, "\n",
        "==============="
      )
    }

    if (all(abs(c(grad_beta, grad_gamma)) <= abstol)) break
  }

  message("Finishing after ", i, " iterations")
  invisible(model)
}





# boost-function first try -----------------------------------------------------------
boost_levin <- function(model,
                        weight_learn = 0.1,
                        m_stop = 1000,
                        verbose = FALSE) {
  print("Boost-Function called. Nothing implemented yet. Sorry, Dude.")
  invisible(model)
}

# Boost Sebastian:

boost_sebastian <- function(model,
                        weight_learn = 0.001,
                        m_stop = 500,
                        verbose = FALSE) {
  print("nothing yet")
  invisible(model)
}

