rm(list=ls())
source("simulation_studies/helper_plot_model.R")
source("simulation_studies/helper_init_model.R")
source("simulation_studies/helper_calc_mse.R")
library(asp20boost)
library(gamboostLSS)

beta_1 <- c(0,1)
gamma_1 <- c(-3,2)
beta_2 <- c(0,1)
gamma_2 <- c(-0.5,0)
beta_3 <- c(0,1)
gamma_3 <- c(-3,3)
beta_4 <- c(0,1)
gamma_4 <- c(3,-3)



mod_1 <- init_model(beta_1, gamma_1, 1000)
plot_model(beta_1, gamma_1, 1000)
gradient_boost(mod_1)
calc_MSE(beta_1, gamma_1, mod_1$beta, mod_1$gamma)


# mod_2 <- init_model(beta_2, gamma_2, 1000)
# plot_model(beta_2, gamma_2, 1000)
# gradient_boost(mod_2)
# calc_MSE(beta_2, gamma_2, mod_2$beta, mod_2$gamma)


# mod_3 <- init_model(beta_3, gamma_3, 1000)
# plot_model(beta_3, gamma_3, 1000)
# gradient_boost(mod_3)
# mod_3$beta
# mod_3$gamma


# mod_4 <- init_model(beta_4, gamma_4, 1000)
# plot_model(beta_4, gamma_4, 1000)
# gradient_boost(mod_4)
# mod_4$beta
# mod_4$gamma



# ausprobieren --------------------------------------












