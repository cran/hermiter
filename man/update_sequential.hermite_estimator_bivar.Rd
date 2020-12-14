% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hermite_estimator_bivar.R
\name{update_sequential.hermite_estimator_bivar}
\alias{update_sequential.hermite_estimator_bivar}
\title{Updates the Hermite series based estimator sequentially}
\usage{
\method{update_sequential}{hermite_estimator_bivar}(this, x)
}
\arguments{
\item{this}{A hermite_estimator_bivar object.}

\item{x}{A numeric vector of length 2. A bivariate observation to be 
incorporated into the estimator.}
}
\value{
An object of class hermite_estimator_bivar.
}
\description{
This method can be applied in sequential estimation settings.
}
\examples{
hermite_estimator <- hermite_estimator_bivar(N = 10, standardize = TRUE)
hermite_estimator <- update_sequential(hermite_estimator, x = c(1,2))
}