% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{gradient_boost}
\alias{gradient_boost}
\title{Gradient Boost Algorithm for Location-Scale Regression}
\usage{
gradient_boost(
  model,
  stepsize = c(0.01, 0.1),
  maxit = 1000,
  abstol = 0.001,
  componentwise = FALSE,
  verbose = FALSE,
  plot = FALSE
)
}
\arguments{
\item{model}{A \link{LocationScaleRegressionBoost} object.}

\item{stepsize}{The learning rate for the parameter updates.}

\item{maxit}{The maximum number of iterations.}

\item{abstol}{The absolute convergence tolerance. The algorithm stops if the
absolute value of the unit-wise gradients all drop below this value.}

\item{verbose}{Whether to print the progress of the algorithm.}

\item{plot}{Whether to plot the result of the algorithm.}
}
\value{
The updated model, invisibly.
}
\description{
This function optimizes the log-likelihood of the given location-scale
regression model by gradient boosting. It is designed to work specifically
together with the LocationScaleRegressionBoost model class.
}
\examples{
y <- rnorm(30)
model <- LocationScaleRegressionBoost$new(y ~ 1)
gradient_boost(model)

}
