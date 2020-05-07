#devtools::install_gitlab("asp20/asp20model", host = "gitlab.gwdg.de", force=TRUE)
library(asp20model)
library(tidyverse)
library(R6)

LocationScaleRegressionBoost <- R6Class("LocationScaleRegression",
                                        inherit = LocationScaleRegression,
                                        public = list(
                                          
                                          X_tr_X_inv = numeric(),
                                          Proj_M_X = numeric(),
                                          Z_tr_Z_inv = numeric(),
                                          Proj_M_Z = numeric(),
                                          y_upperTube = numeric(),
                                          y_downTube = numeric(),
                                          
                                          initProj_M = function() {
                                            #𝑃𝑟𝑜𝑗𝑀𝑎𝑡𝑟𝑖𝑥=(𝑋𝑇𝑋)−1𝑋𝑇
                                            #https://stats.stackexchange.com/questions/154485/least-squares-regression-step-by-step-linear-algebra-computation?noredirect=1&lq=1
                                            self$X_tr_X_inv <- solve(t(private$X) %*% private$X)    
                                            self$Proj_M_X <- self$X_tr_X_inv %*% t(private$X)
                                            
                                            self$Z_tr_Z_inv <- solve(t(private$Z) %*% private$Z)    
                                            self$Proj_M_Z <- self$Z_tr_Z_inv %*% t(private$Z)
                                            
                                          },
                                          negloglik = function() {
                                            location <- self$fitted_location
                                            scale <- self$fitted_scale
                                            
                                            sum(dnorm(private$y, location, scale, log = TRUE))
                                          }
                                          
                                          
                                          
                                          
                                        ),
                                        active = list(
                                          lstSqrResid = function()
                                          {
                                            
                                            self$Proj_M_X %*% self$resid()
                                          },
                                          
                                          derivative1 = function()
                                          {
                                            
                                            sum(self$resid()^2*private$Z[,1]*exp(-2*(private$Z[,1]*self$gamma[1]+self$gamma[2]*private$Z[,2]))-private$Z[,1])
                                            
                                            
                                            
                                          },
                                          derivative2 = function()
                                          {
                                            sum((self$resid()^2)*private$Z[,2]*exp(-2*(private$Z[,2]*self$gamma[2]+self$gamma[1]*private$Z[,1]))-private$Z[,2])
                                          },
                                          
                                          
                                          
                                          getX = function()
                                          {
                                            private$X
                                          },
                                          getZ = function()
                                          {
                                            private$Z
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
  v <- stepsize
  
  
  #check if location scale model
  for (i in seq_len(maxit)) {
    model$beta <- model$beta + v*model$lstSqrResid
    

    #message(par_msg2)
    
    grad_beta <- model$grad_beta()
    
    newgamma<-c(model$gamma[1] + v*model$derivative1,model$gamma[2] + v*model$derivative2)
    model$gamma<-newgamma
    
    #model$gamma <- c(model$gamma[1] + v*model$derivative1,model$gamma[2] + v*model$derivative2)
    
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
    
    
    
    if (all(abs(c(grad_beta, grad_gamma)) <= abstol))
    {
      message("abs stop at")
      message("Iteration:      ", i, "\n")
      break 
    }
    
    
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
model$initProj_M()
model$gamma<-c(10,20)
model$beta<-c(-20,50)

boostJohannes(model, 
              stepsize = 0.001, maxit = 5000,
              abstol = 0.001,
              verbose = TRUE)


#
####
#

model$loglik()

dftest = data.frame(x,y,
                    "resid"=model$resid(),"fitted_location" = model$fitted_location,
                    "man_fitted_location" = model$beta[2]*x+model$beta[1] , 
                    "y_hat"=model$fitted_location+0.01045657,
                    "manResid" = y-model$beta[2]*x+model$beta[1] , 
                    "fitted_scale" = model$fitted_scale,
                    "upperTube" = model$fitted_location + model$fitted_location*model$fitted_scale,
                    "downTube" =  model$fitted_location - model$fitted_location*model$fitted_scale)



#"upperTube" = model$gamma[1]+  model$fitted_location + model$fitted_location*model$fitted_scale,
#"downTube" =-1*model$gamma[1] + model$fitted_location - model$fitted_location*model$fitted_scale)
#sd = exp(-3 + 2 * x)



plot(x,y)
lines(dftest$x,dftest$upperTube, type = "p", col="blue")
lines(dftest$x,dftest$downTube, type = "p", col="green")
lines(x,model$fitted_location)
abline(0, 0) 

model$fitted_location + model$fitted_location*model$fitted_scale


model$gamma

model$grad_beta()
model$grad_gamma()
