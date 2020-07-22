context("model class creation")

source("helper_model_class_creation.R")
# for ad-hoc loading use: source("tests/testthat/helper_model_class_creation.R")
#creates one simple model and one large model


test_that("inheritance works", {
  model_1 <- model_simple$clone()
  model_2 <- model_large$clone()
  expect_is(model_1, "LocationScaleRegression")
  expect_is(model_1, "LocationScaleRegressionBoost")
  expect_is(model_2, "LocationScaleRegression")
  expect_is(model_2, "LocationScaleRegressionBoost")
})


test_that("beta vector initialized correctly", {
  model_1 <- model_simple$clone()
  model_2 <- model_large$clone()
  expect_equal(length(model_1$beta), 2)
  expect_equal(model_1$beta, c(0,0 ))
  expect_equal(length(model_2$beta), 6)
  expect_equal(model_2$beta, c(0,0,0,0,0,0))
})


test_that("gamma vector initialized correctly", {
  model_1 <- model_simple$clone()
  model_2 <- model_large$clone()
  expect_equal(length(model_1$gamma), 2)
  expect_equal(model_1$gamma, c(0,0 ))
  expect_equal(length(model_2$gamma), 6)
  expect_equal(model_2$gamma, c(0,0,0,0,0,0))
})
