source("simulation_studies/helper_init_large_model.R")


# init_1 <- init_large_model(c(5,10), c(0.5,1), 5, -0.2) #works quite
# init_2 <- init_large_model(c(5,10), c(0.5,1), 10, -1) #bad performance
# init_3 <- init_large_model(c(3,6), c(0.2,0.8), 10, -0.2) #works quite
init <- init_large_model(c(5,10), c(0.01,0.08), 50, 0)

mod_1 <- init$model
init$fitted_response
init$fitted_sd

 init$beta
init$gamma



gradient_boost(mod_1, stepsize = c(0.1,0.01), maxit = 1000, componentwise = TRUE, verbose = TRUE)

init$beta
mod_1$beta

init$gamma
mod_1$gamma
