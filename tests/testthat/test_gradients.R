context("Eta-Gradients")

library(numDeriv)

set.seed(1337)
n <- 500
x <- runif(n)
y <- rnorm(n, mean = x, sd = exp(-3 + 2* x))
mod_test <- LocationScaleRegressionBoost$new(y ~ x, ~ x)
gradient_boost(mod_test)



loglik <- function(a) {
  sum(dnorm(a, mod_test$fitted_location, mod_test$fitted_scale, log = TRUE))
}

# test_that("mu gradients work", {
#   model$clone()
#   expect_equivalent(
#     model$gradients_loglik_mu(),
#     grad(func = loglik, x = 1/2 * mod_test$eta_mu(), method = "simple")
#   )
# })
#
# test_that("sigma gradients work", {
#   model$clone()
#   expect_equivalent(
#     model$gradients_loglik_sigma(),
#     grad(func = loglik, x = mod_test$eta_sigma(), method = "simple")
#   )
# })
#
#

# test_that("sigma gradients work", {
#   expect_equivalent(
#     model$gradients_loglik_sigma(),
#     grad(func = model$loglik(), x = model$eta_sigma())
#   )
# })

# f <- function(x) {
#   model <- model$clone()
#   model$gamma <- x
#   model$loglik()
# }
#
# test_that("gamma gradient works", {
#   expect_equivalent(model$grad_gamma(), grad(f, model$gamma))
# })
