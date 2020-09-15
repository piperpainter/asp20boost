set.seed(1337)

n <- 500
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y <- rnorm(n, x1 + x3, exp(-3 + x2 + x3))
model <- LocationScaleRegressionBoost$new(y ~ x1 + x3, ~ x2 + x3)
gradient_boost(model,stepsize = c(0.01, 0.1), maxit = 1000, abstol = 0.0001,componentwise = TRUE, verbose = FALSE, plot=FALSE)
