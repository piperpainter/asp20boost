# context("functionality boosting algorithm")
#
# source("helper_model_class_creation.R")
# # creates one simple and one large model
#
# # for ad-hoc loading use: source("tests/testthat/helper_model_class_creation.R")
#
#
# # estimate models
# gradient_boost(model_simple, stepsize = 0.01, maxit = 10000, abstol = 0.0001, verbose = TRUE)
# gradient_boost(model_large, stepsize = 0.001, maxit = 2000, abstol = 0.0001, verbose = FALSE)
#
#
# test_that("beta estimates are correct in simple model",{
#   model_1 <- model_simple$clone()
#   expect_equal(model_1$beta, c(0, 1))
# })
#
# test_that("gamma estimates are correct in simple model",{
#   model_1 <- model_simple$clone()
#   expect_equal(model_1$gamma, c(0, 1))
# })
