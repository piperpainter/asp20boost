library(numDeriv)

set.seed(1337)

n <- 500
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- rnorm(n, x1 + x3, exp(-3 + x2 + x3))

model <- LocationScaleRegressionBoost$new(y ~ x1 + x3, ~ x2 + x3)

f <- function(x) {
  model <- model$clone()
  model$beta <- x
  model$loglik()
}

test_that("beta gradient works", {
  expect_equivalent(model$grad_beta(), grad(f, model$beta))
})

f <- function(x) {
  model <- model$clone()
  model$gamma <- x
  model$loglik()
}

test_that("gamma gradient works", {
  expect_equivalent(model$grad_gamma(), grad(f, model$gamma))
})


#test_that("gradient boost does not throw an error", {
#  expect_error(gradient_boost(model, verbose = TRUE), NA)
#})

