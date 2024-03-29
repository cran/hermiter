% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hermite_estimator_univar.R
\name{summary.hermite_estimator_univar}
\alias{summary.hermite_estimator_univar}
\title{Summarizes univariate hermite_estimator object.}
\usage{
\method{summary}{hermite_estimator_univar}(object, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
\item{object}{A hermite_estimator_univar object.}

\item{digits}{A numeric value. Number of digits to round to.}

\item{...}{Other arguments passed on to methods used in summary.}
}
\description{
Outputs key parameters of a univariate hermite_estimator object along with
estimates of the mean, standard deviation and deciles of the data that
the object has been updated with.
}
