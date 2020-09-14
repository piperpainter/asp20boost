rm(list = ls())
library(asp20boost)
options(scipen=999)
source("simulation_studies/helper_init_large_model.R")
windows()


# initialize model
result <- init_large_model()
mod_1 <- result$model
X <- result$design
y <- result$response
fitted_loc <- result$fitted_response
fitted_scale <- result$fitted_sd

# check model via plotting ------------------------------------------------------------------
par(mfrow = c(4,3))
for(i in 1:11) plot(X[,i], y, ylim = c(0,200))
# only one effect is clearly indicated, the others are too small compared to the very high response
# the effect is also indicated if heteroskedastic sds are added. Hence, signal-noise ratio is OK.
# ALSO FOR OTHER effects? or only for the exceptionally high one?



# check model via linear predictors ---------------------------------------------------------
fitted_loc
fitted_scale

# estimate huge model -----------------------------------------------------------------------
#gradient_boost(mod_1, maxit = 1000, componentwise = TRUE)



# ausprobieren ------------------------------------------------------------------------------
mod_1$beta
mod_1$gamma
