#' Gradient Boosting for Location-Scale Regression Models HAI HAI
#'
#' The package provides an R6-class `LocationScaleRegressionBoost` for model
#' specification, and a function `gradient_boost` for boosting the specified
#' model object.
#'
#' @usage
#' model <- LocationScaleRegressionBoost(formula, ...)
#' gradient_boost(model, ...)
#'
#' @details
#' This package resulted from a series of student projects in the context of the
#' seminar *Advanced Statistical Programming 2020* of the University of
#' Göttingen. Thus, it is part of a bundle of the `asp20__` packages. The main
#' goal of these projects is implement a method to estimate location-scale
#' regression models. Each project's package implements a different method. All
#' packages, including `asp20boost`, build upon the `asp20model` package.
#' Specifically, an R6-class with inherited properties from `asp20model` is
#' `asp20model` expanded in the `asp20boost`-package.
#'
#'
#' @author Johannes Strauß, Levin Wiebelt, Sebastián Aristizábal
#'
#' @import R6
#' @import asp20model
#'
#' @docType package
#' @name asp20boost
NULL
