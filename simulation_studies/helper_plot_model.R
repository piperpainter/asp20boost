# set up main plotting function ----------------------------------
plot_model <- function(beta, gamma, n, plot_title) {

  # set up linear predictor functions -------------------------------
  eta_mu <- function(x) beta[1] + beta[2] * x
  eta_sigma <- function(x) gamma[1] + gamma[2] * x

  # set up random variables -----------------------------------------
  x <- runif(n)
  y <- eta_mu(x) + rnorm(n, sd = exp(eta_sigma(x)))

  # configure plot --------------------------------------------------
  plot(x, y, col = "grey", cex = 0.5, pch = 16, main = plot_title)
  abline(beta[1], beta[2], lwd = 2, col = 1)
  curve(eta_mu(x) + 1 * exp(eta_sigma(x)), -0.1, 1.1, add = TRUE, col = 1, lty = 3, lwd = 2)
  curve(eta_mu(x) - 1 * exp(eta_sigma(x)), -0.1, 1.1, add = TRUE, col = 1, lty = 3, lwd = 2)
  abline(h = beta[1], lty = 2, col = 4)

  arrows(x0 = 1, x1 = 1,
         y0 = eta_mu(1),
         y1 = eta_mu(1) + 1 * exp(eta_sigma(1)),
         length = 0.05, code = 3, col = 2, lty = 1, lwd = 2)
  arrows(x0 = 0.985, x1 = 0.985,
         y0 = eta_mu(1),
         y1 = beta[1],
         length = 0.05, code = 3, col = 4, lty = 1, lwd = 2)

  #legend("topleft", legend=c("location slope", "standard deviation at x=1"),
  #       col=c(4, 2), lty=c(1,1), cex=0.8)

}
