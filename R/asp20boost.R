#' Gradient Boosting for Location-Scale Regression Models
#'
#' The package provides an R6-class `LocationScaleRegressionBoost` for model
#' specification, and a function `gradient_boost` for boosting the specified model object.
#'
#' @usage
#' example_model <- LocationScaleRegressionBoost(formula, ...)
#' gradient_boost(example_model, ...)
#'
#' @details
#' This package resulted from a series of student projects. Thus, it is part of a
#' bundle of packages all starting with `asp20...`. The overall goal of these packages
#' is to estimate location-scale regression models, while each project/package employs
#' different methods. All packages, including `asp20boost` build on the central
#' package `asp20model`. Here, an R6-class is specified, which is inherited and expanded
#' in the `asp20boost`-package.
#'
#' @author Johannes Strauß, Levin Wiebelt, Sebastiàn Aristizabal
#'
#' @import R6
#' @import asp20model
#'
#' @docType package
#' @name asp20boost
NULL
