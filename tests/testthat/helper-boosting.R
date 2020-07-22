set.seed(1337)

n <- 500
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- rnorm(n, x1 + x3, exp(-3 + x2 + x3))
model <- LocationScaleRegressionBoost$new(y ~ x1 + x3, ~ x2 + x3)
gradient_boost(model, stepsize = 0.001, maxit = 10000, abstol = 0.0001, verbose = FALSE)
