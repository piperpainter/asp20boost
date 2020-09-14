helper_plot_model_high_p <- function(beta, gamma, X, index){
  
  # set up linear predictor functions -------------------------------
  eta_mu <- function(x) beta[1] + beta[index] * x
  eta_sigma <- function(x) gamma[1] + gamma[index] * x
  
  x <- X[ , index]
  
  plot(x, y) # plot y against column vectors of X
  curve(mean(y) + eta_mu(x) + 1 * exp(eta_sigma(x)), -0.1, 1.1, add = TRUE, col = 1, lty = 3, lwd = 2)
  curve(mean(y) + eta_mu(x) - 1 * exp(eta_sigma(x)), -0.1, 1.1, add = TRUE, col = 1, lty = 3, lwd = 2)
  
}
