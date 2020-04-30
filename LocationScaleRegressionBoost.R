library(asp20model)
library(tidyverse)


LocationScaleRegressionBoost <- R6Class("LocationScaleRegression",
                                        inherit = LocationScaleRegression,
                                        
                                        private = list(
                                          residScale = numeric(),
                                          update_gamma = function(value) {
                                            private$.gamma <- value
                                            private$fitted$scale <- exp(drop(private$Z %*% self$gamma))
                                            private$residScale <- private$y - private$fitted$scale
                                            invisible(self)
                                          }
                                        ),
                                        
                                        
                                        active = list(
                                          lstSqrResid = function()
                                          {
                                            
                                            X_tr_X_inv <- solve(t(private$X) %*% private$X)    
                                            Proj_M <- X_tr_X_inv %*% t(private$X)
                                            Proj_M %*% private$resid
                                          },
                                          
                                          lstSqrResidScale = function()
                                          {
                                            #
                                            # How to calculate the factor to change for gamma?
                                            #
                                            X_tr_X_inv2 <- solve(t(private$X) %*% private$X)    
                                            Proj_M2 <- X_tr_X_inv2 %*% t(private$X)
                                            #Proj_M2 %*% private$residScale
                                            #(resid / scale)^2 - 1
                                          }
                                          
                                        )
                                        
                                        
)




boostJohannes = function(model,
                         stepsize = 0.001,
                         maxit = 10000,
                         abstol = 0.1,
                         verbose = FALSE) {
  grad_beta <- model$grad_beta()
  grad_gamma <- model$grad_gamma()
  
  message("boost js")
  v <- 0.1
  
  for (i in seq_len(maxit)) {
    model$beta <- model$beta + v*model$lstSqrResid
    #message(model$beta)
    grad_beta <- model$grad_beta()
    ##todo gamma
    model$gamma <- model$gamma + v*model$lstSqrResid
    #message(model$loglik())
    
  }
  
  invisible(model)
}

##
##Test
##
n <- 500
x <- runif(n)
y <- x + rnorm(n, sd = exp(-3 + 2 * x))
model <- LocationScaleRegressionBoost$new(y ~ x, ~ x)
model$gamma<-c(0,0)
model$beta<-c(0,0)
model$loglik()
boostJohannes(model)
model$beta
model$gamma
model$grad_beta()

model$lstSqrResidScale
model$lstSqrResid


plot(x = x , y = y, type = 'p')

abline(model$beta)
abline(model$gamma)


model$grad_beta()
model$grad_gamma()
