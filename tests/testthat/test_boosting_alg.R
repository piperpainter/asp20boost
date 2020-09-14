context("functionality boosting algorithm")

source("helper_test_boosting.R")
# creates one simple boosting model

# for ad-hoc loading use: source("tests/testthat/helper_model_class_creation.R")


test_that("length of covariates equals beta parameters",{
  model_test <- model$clone()
  expect_equal(length(model_test$beta), 3)
})

test_that("length of covariates equals gamma parameters",{
  model_test <- model$clone()
  expect_equal(length(model_test$gamma), 3)
})

test_that("model beta step size equals boosting parameter",{
  model_test <- model$clone()
  expect_equal(model_test$stepsize_beta, 0.01)
})

test_that("model gamma step size equals boosting parameter",{
  model_test <- model$clone()
  expect_equal(model_test$stepsize_gamma, 0.1)
})

test_that("beta values are finite values",{
  model_test <- model$clone()
  expect_equal(all(is.finite(model_test$beta)==TRUE), TRUE)
})
test_that("gama values are finite values",{
  model_test <- model$clone()
  expect_equal(all(is.finite(model_test$gamma)==TRUE), TRUE)
})


