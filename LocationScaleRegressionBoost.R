#devtools::install_gitlab("asp20/asp20model", host = "gitlab.gwdg.de", force=TRUE)
library(asp20model)
library(tidyverse)
library(R6)

LocationScaleRegressionBoost <- R6Class("LocationScaleRegression",
                                        inherit = LocationScaleRegression,
                                        private = list(
                                          update_beta = function(value) {
                                          private$.beta <- value
                                          self$fitted_location <- drop(private$X %*% self$beta)
                                          private$.resid <- private$y - self$fitted_location
                                          #message("update beta \n")
                                          #message( self$beta)
                                          invisible(self)
                                        }
                                        ),
                                   active = list(
                                          lstSqrResid = function()
                                          {
                                            
                                            X_tr_X_inv <- solve(t(private$X) %*% private$X)    
                                            Proj_M <- X_tr_X_inv %*% t(private$X)
                                            Proj_M %*% self$resid()
                                          },
                                          lstSqrResidZ = function()
                                          {
                                            
                                            X_tr_X_inv <- solve(t(private$X) %*% private$X)    
                                            Proj_M <- X_tr_X_inv %*% t(private$X)
                                            
                                            Proj_M %*% (self$resid("deviance"))
                                          }
                                          )
                                        
                                        
)




boostJohannes = function(model,
                         stepsize = 0.001,
                         maxit = 1000,
                         abstol = 0.001,
                         verbose = TRUE) {
  grad_beta <- model$grad_beta()
  grad_gamma <- model$grad_gamma()
  
  message("boost js")
  v <- 0.01
  
  
  #check if location scale model
  for (i in seq_len(maxit)) {
    model$beta <- model$beta + v*model$lstSqrResid
    grad_beta <- model$grad_beta()
    
    model$gamma <- model$gamma + v*model$lstSqrResidZ
   
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
  
  invisible(model)
}

##
##Test
##
##
set.seed(1337)
n <- 500
x <- runif(n)
y <- x + rnorm(n, sd = exp(-3 + 2 * x))
model <- LocationScaleRegressionBoost$new(y ~ x, ~ x)
model$gamma<-c(0,0)
model$beta<-c(0,0)
model$loglik()
boostJohannes(model, 
              stepsize = 0.0001, maxit = 40000,
              abstol = 0.001,
              verbose = TRUE)


#
####
#





df = data.frame(x,y,
                "resid"=model$resid(),"fitted_location" = model$fitted_location,
                "man_fitted_location" = model$beta[2]*x+model$beta[1] , 
                "y_hat"=model$fitted_location+0.01045657,
                "manResid" = y-model$beta[2]*x+model$beta[1] , 
                "fitted_scale" = model$fitted_scale,
                "upperTube" = model$gamma[1]+  model$fitted_location + model$fitted_location*model$fitted_scale,
                "downTube" =-1*model$gamma[1] + model$fitted_location - model$fitted_location*model$fitted_scale)

#"upperTube" = model$gamma[1]+  model$fitted_location + model$fitted_location*model$fitted_scale,
#"downTube" =-1*model$gamma[1] + model$fitted_location - model$fitted_location*model$fitted_scale)
#sd = exp(-3 + 2 * x)

plot(x, y, ylim = c(-2,2))
lines(x,df$fitted_location,type="l")
lines(x,df$upperTube,type="p",col="blue")
lines(x,df$downTube,type="p",col="green")


