
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