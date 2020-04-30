
#First run only
#install.packages("devtools")
#devtools::install_gitlab("asp20/asp20model", host = "gitlab.gwdg.de")


library(asp20model)
library(tidyverse)

### Example Data

n <- 500
x <- runif(n)
y <- x + rnorm(n, sd = exp(-3 + 2 * x))
plot(x, y)
abline(0, 1, lwd = 2)
curve(x + 1.96 * exp(-3 + 2 * x), -0.1, 1.1, add = TRUE)
curve(x - 1.96 * exp(-3 + 2 * x), -0.1, 1.1, add = TRUE)

###

### Build simple linear model

simpleModel <- lm(x ~ y)
summary(simpleModel)
plot(x = x , y = y, type = 'p')
abline(simpleModel)


### Plot Resdiuals
simpleModelRes = resid(simpleModel)
plot(y, simpleModelRes, ylab="Residuals", xlab="X",  main="Residual Plot") 
abline(0, 0) 
#It shows Heteroscedasticity
#Note: Residuals should be standardized / studentized, P.184
#
LocationScaleRegression <- R6Class("LocationScaleRegression")

#
#Extend LocationScaleRegression Model 
#
gradient_descentBoost <- function(model,
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

gradient_descentBoost(model)
model$grad_gamma()
model <- LocationScaleRegression$new(y ~ x, ~ x)

