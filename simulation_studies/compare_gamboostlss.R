library(asp20boost)
library(gamboostLSS)

# set true model parameters
beta_0  <- 0
beta_1  <- 1
gamma_0 <- -3
gamma_1 <- 2



# simulate data -------------------------
set.seed(1337)
n <- 1000
x <- runif(n)
y <- beta_0 + beta_1 * x + rnorm(
  n,
  sd = exp(gamma_0 + gamma_1 * x)
)

xy_df <- cbind(x,y)
head(xy_df)

# estimate data by boosting ------------
mod_1 <- LocationScaleRegressionBoost$new(y ~ x, ~ x)
gradient_boost(mod_1, maxit = 1000)


# estimate with external package --------
ctrl <- boost_control(trace = FALSE,
                      mstop = c(mu = 1000, sigma = 1000),
                      nu = c(mu = 0.01, sigma = 0.001)
                      )
mod_2 <- gamboostLSS(formula = y ~ x,
                     control = ctrl)
eta_mu <- fitted(mod_2)$mu
eta_sigma <- fitted(mod_2)$sigma


lm(eta_mu ~ x)
lm(eta_sigma ~ x)

